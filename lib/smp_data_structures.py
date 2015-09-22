##
# This module conatins classes and functions used to 
# support the simulation of a semi-Markov process.
# 
# The underlying failure/repair distributions are assumed
# to be Weibull, which is Exponential when shape=1.
#
# Kevin Greenan (kmgreen@cs.ucsc.edu)
#
#

import mpmath
from mpmath import mpf
from mpmath import ln
from mpmath import findroot
import random


##
# Set precision used by the MP math lib
#
mpmath.mp.prec += 100
mpmath.mp.dps = 100

##
# Contains parameters, distribution functions and hazard rate function
# for a 3-parameter Weibull distribution based on shape, scale and location.
#
class Weibull:
    ##
    # Construct a Weibull object by specifying shape, scale and location.
    #
    # Note: when shape == 1, this is an Exponential distribution
    #
    def __init__(self, shape=1, scale=1, location=0):
        self.shape = mpf(shape)
        self.scale = mpf(scale)
        self.location = mpf(location)

    ##
    # Get the probability density of Weibull(shape, scale, location) at x 
    # 
    # @param x: random variable, most likely a time
    # @return density of Weibull(shape, scale, location) at x
    #
    def pdf_eval(self, x):
        if x < 0:
            return 0
        elif x < self.location:
            return 0
        else:
            x = mpf(x)
            a = self.shape/self.scale
            b = (x-self.location)/self.scale
            b = mpmath.power(b, self.shape-1)
            c = mpmath.exp(-mpmath.power(((x-self.location)/self.scale), self.shape))

            return a * b *c
    ##
    # Return the probability P(X <= x) or the probability that a 
    # random variable is less than or equal to the parameter x.
    # 
    # The returned value represents the probability that there is 
    # a 'failure' at or before x, which is most likely a time.
    #
    # @param x: variable, most likely a time
    # @return: probability of failure before or at x
    #
    def cdf_eval(self, x):
        x = mpf(x)

        if x < self.location:
            return 0

        return mpf(1) - mpmath.exp(-mpmath.power(((x-self.location)/self.scale), self.shape))

    ##
    # Return the hazard rate at x.  The hazard rate is interpreted 
    # as the instantaneous failure rate at x.
    #
    # Note: If shape == 1, then this value will be constant for all x
    #
    # @param x: variable, most likely a time
    # @return instantaneous failure rate at x 
    #
    def hazard_rate(self, x):
        if x < self.location:
            return 0
        elif self.shape == 1:
            return mpf(1) / self.scale
        else:
            return abs(self.pdf_eval(x) / (mpf(1) - self.cdf_eval(x)))
            
    ##
    # When the shape parameter is not 1, then the hazard rate will
    # change with time.  When simulating a semi-Markov process using
    # uniformization (which is a method we use), the maximum failure
    # rate for every distribution over all possible times is required.
    # This function evaluates the hazard rate at discrete points and
    # returns the maximum hazard rate for a given mission time (i.e. 
    # max simulation time).
    #
    # @param mission_time: expected maximum simulation time
    # @return maximum possible hazard rate for [0, mission_time]
    #
    #
    def get_max_hazard_rate(self, mission_time):
        max = mpf(0)

        if self.shape == 1:
            return mpf(1) / self.scale

        for i in range(1, mission_time, int(0.1 * mission_time)):
            curr_h_rate = self.hazard_rate(i)
            if curr_h_rate > max:
                max = curr_h_rate
            elif curr_h_rate == mpf('nan'):
                break
    
        return max

    ##
    # Get min hazard rate
    #    
    def get_min_hazard_rate(self, mission_time):
        min = mpf(1)

        if self.shape == 1:
            return mpf(1) / self.scale

        for i in range(0, mission_time, int(0.1 * mission_time)):
            curr_h_rate = self.hazard_rate(i)
            if curr_h_rate < min:
                max = curr_h_rate
            elif curr_h_rate == mpf('nan'):
                break
    
        return min

    ##
    # Draw a random value from this distribution
    #
    def draw(self):
        return random.weibullvariate(self.scale, self.shape) + self.location
    
    ##
    # Draw from a lower truncated Weibull distribution
    # (Reject all samples less than 'lower')
    #
    def draw_truncated(self, lower):
        val = self.draw()
        while val <= lower:
            print "%.2f %.2f" % (lower, val)
            val = self.draw()
            
        return val
    
    ##
    # Draw using the inverse transform method. 
    # This method draws from the waiting time 
    # based on the CDF built from the distribution's
    # hazard rates
    #
    def draw_inverse_transform(self, curr_time):
        U = random.uniform(0,1)
        while U == 0:
            U = random.uniform(0,1)
        draw = ((-(self.scale**self.shape)*ln(U)+((curr_time)**self.shape))**(1/self.shape) - (curr_time))
        
        return abs(draw)
            
##
# This class encapsulates the state of a component under simulation.
# Each component is given failure and repair distributions for 
# the entire component.
#
# A component may be in one of two states: OK (operational, no failures) or
# FAILED (entire component is failed). 
#
#
class Component:

    ##
    # The three possible states
    #
    STATE_OK = "state ok"
    STATE_FAILED = "state failed"

    ##
    # Possible failure events
    #
    EVENT_COMP_FAIL = "component failure"
    EVENT_COMP_REPAIR = "component repair"

    ##
    # A component is constructed by specifying the appropriate failure/repair distributions.
    #
    # The component fail/repair distributions must be specified.
    #
    # This function will set the component state to OK and set all clocks to 0.
    # 
    # init_clock *must* first be called in order to use this object in simulation.
    #
    def __init__(self, component_fail_distr, component_repair_distr):
        # Current state 
        self.state = self.STATE_OK

        # Last "global" clock update
        self.last_time_update = mpf(0)

        # Global begin time of this component
        self.begin_time = mpf(0)

        # Local repair time of this component
        self.repair_clock = mpf(0)

        # Local (relative) clock of this component
        self.clock = mpf(0)
        
        self.repair_start = mpf(0)

        # Failure and repair distributions
        self.component_fail_distr = component_fail_distr
        self.component_repair_distr = component_repair_distr
        

    ##
    # Set the last clock update to the current simulation time and initialize 
    # t_0 for this component (begin_time).
    #
    # @param curr_time: t_0 of this component
    #
    def init_clock(self, curr_time):
        self.last_time_update = curr_time
        self.begin_time = curr_time
        self.clock = mpf(0)
        self.repair_clock = mpf(0)
        self.repair_start = mpf(0)

    ##
    # Set the state of this component to OK
    #
    def init_state(self):
        self.state = self.STATE_OK

    ##
    # Update component clocks.  There are three main clocks to update:
    # the component clock, the repair clock (if there is an ongoing repair)
    # and the time of the last clock update (used to update the component clock)
    #
    # The clock member variable is used to get instantaneous failure 
    # rate, while the repair clock is used to get the instantaneous "repair" rate.
    #
    # @param curr_time: current simulation time
    #
    def update_clock(self, curr_time):
        self.clock += (curr_time - self.last_time_update)
        if self.state == self.STATE_FAILED:
            self.repair_clock = (curr_time - self.repair_start)
        else:
            self.repair_clock = mpf(0)

        
        self.last_time_update = curr_time
    ##
    # Get this component's current time    
    #
    # @return current component clock reading
    # 
    def read_clock(self):
        return self.clock
    
    ##
    # Get this component's current time    
    #
    # @return current component clock reading
    # 
    def read_repair_clock(self):
        return self.repair_clock

    ##
    # Get component state
    #
    # @return component state
    #
    def get_curr_state(self):
        return self.state

    ##
    # Fail this component.  Reset list of failed sub-components
    #
    def fail_component(self, curr_time):
        self.state = self.STATE_FAILED
        self.repair_clock = mpf(0)
        self.repair_start = curr_time
    
    ##
    # Repair this component.  
    #
    def repair_component(self):
        self.begin_time = self.last_time_update
        self.clock = mpf(0)
        self.repair_clock = mpf(0)
        self.state = self.STATE_OK
    
    ##
    # Get instantaneous failure rate of this component
    #
    # @return instantaneous whole-component failure rate
    #
    def curr_component_fail_rate(self):
        if self.state == self.STATE_FAILED:
            return mpf(0)

        return self.component_fail_distr.hazard_rate(self.clock)

    ##
    # Get instantaneous repair rate of this component
    #
    # @return instantaneous whole-component repair rate
    #
    def curr_component_repair_rate(self):
        if self.state == self.STATE_OK:
            return mpf(0)
        
        return self.component_repair_distr.hazard_rate(self.repair_clock)

    ##
    # Return sum of instantaneous fail/repair rates
    #
    def inst_rate_sum(self):
        return self.curr_component_fail_rate() + self.curr_component_repair_rate() 
    
class InverseTransformHomogeneousFailRepairRates:
    
    def __init__(self, components):
        self.components = components
        self.fail_shape = self.components[0].component_fail_distr.shape
        self.fail_scale = self.components[0].component_fail_distr.scale
        self.repair_shape = self.components[0].component_repair_distr.shape
        self.repair_scale = self.components[0].component_repair_distr.scale
        
        self.fail_scale_to_shape = self.fail_scale**self.fail_shape
        self.repair_scale_to_shape = self.repair_scale**self.repair_shape
        
        self.curr_avail_devices = None
        self.curr_failed_devices = None
        self.curr_uniform_variate = None
        
        for comp in self.components:
            if comp.component_fail_distr.shape != self.fail_shape:
                print "Failure distribution shapes are different!"
                exit()
            if comp.component_fail_distr.scale != self.fail_scale:
                print "Failure distribution scales are different!"
                exit()
            if comp.component_repair_distr.shape != self.repair_shape:
                print "Repair distribution shapes are different!"
                exit()
            if comp.component_repair_distr.scale != self.repair_scale:
                print "Repair distribution scales are different!"
                exit()
    
    def draw_waiting_time(self, avail_devices, failed_devices):
        self.curr_avail_devices = avail_devices
        self.curr_failed_devices = failed_devices
        self.curr_uniform_variate = random.uniform(0,1)
        
        while self.curr_uniform_variate == 0:
                self.curr_uniform_variate = random.uniform(0,1)
        
        try:
            wait_time = findroot(self.func, [0,100], solver='secant')
        except:
            wait_time = findroot(self.func, [0,100], solver='secant', maxsteps=1000)
            #wait_time = findroot(self.func, [0, 87600] , maxsteps=1000, solver='secant', verify=False)
            
        return abs(wait_time)
    
    def get_clock_value(self, clock, location):
        if clock - location < 0:
            return 0
        else:
            return clock - location
    
    def func(self, x):
    	
        lhs_const = -ln(self.curr_uniform_variate)
        
        for a in self.curr_avail_devices:
            lhs_const += ((self.get_clock_value(self.components[a].read_clock(),self.components[a].component_fail_distr.location) ** self.fail_shape) / -self.fail_scale_to_shape)
        for f in self.curr_failed_devices:
            lhs_const += ((self.get_clock_value(self.components[f].read_repair_clock(),self.components[f].component_repair_distr.location) ** self.repair_shape) / -self.repair_scale_to_shape)
        
        rhs_var=0

        for a in self.curr_avail_devices:
                rhs_var += (((x+self.get_clock_value(self.components[a].read_clock(),self.components[a].component_fail_distr.location))**self.fail_shape) / -self.fail_scale_to_shape)
        for f in self.curr_failed_devices:
                rhs_var += (((x+self.get_clock_value(self.components[f].read_repair_clock(),self.components[f].component_repair_distr.location))**self.repair_shape) / -self.repair_scale_to_shape)
        
        
        return rhs_var-lhs_const

        

def test():
    # Basic test of the Weibull functions
    w = Weibull(shape=mpf(2.0), scale=mpf(12), location=6)

    print "Weibull(%s,%s,%s): " % (w.shape, w.scale, w.location)

    sum = 0

    for i in range(10000):
        sum += w.draw_truncated(6)

    print "MEAN: ", mpf(sum) / 10000.
    
    print w.draw_inverse_transform(0)
    print w.draw_inverse_transform(0)
    print w.draw_inverse_transform(0)
    print w.draw_inverse_transform(0)

    print "Max hazard rate is %e\n" % w.get_max_hazard_rate(100)
    
    for i in range(0,200,5):
        print "CDF at time %d is %f\n" % (i, w.cdf_eval(i))
    
    w = Weibull(shape=mpf(1.0), scale=mpf(120000))

    print "Bunch of draws:"
    for i in range(10):
        print w.draw_inverse_transform(1000000)

    print "Weibull(%s,%s,%s): " % (w.shape, w.scale, w.location)

    print "Max hazard rate is %e\n" % w.get_max_hazard_rate(1000)
    
    for i in range(0,1000,100):
        print "Hazard rate at time %d is %e\n" % (i, w.hazard_rate(i))
    

if __name__ == "__main__":
    test()
