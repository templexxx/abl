from poisson_process import *
from simulation import *
from numpy import random

class UniformizationBFBOpt(Simulation):
     
    ##
    # __init__() from Simulation 
    #
     
    ##
    # Initialize the simulation
    #
    def init(self):
         ##
         # Get the maximum outgoing rate of any state over the mission time (beta)
         #
         self.poisson_rate = mpf(0)

         ##
         # The minimum hazard rate for rebuilds
         #
         self.min_rebuild_rate = mpf(1)

         logging.debug("Finding maximum rate for the Poisson process")
         
         for component in self.components:
              dist = component.component_fail_distr     
              self.poisson_rate += dist.get_max_hazard_rate(self.mission_time)

              logging.debug("Current maximum rate is %e" % self.poisson_rate)

              dist = component.component_repair_distr     
              self.poisson_rate += dist.get_max_hazard_rate(self.mission_time)

              logging.debug("Current maximum rate is %e" % self.poisson_rate)
      
         self.poisson_rate = self.components[0].component_repair_distr.get_max_hazard_rate(self.mission_time) * self.eras_code.m
         self.poisson_rate += self.components[0].component_fail_distr.get_max_hazard_rate(self.mission_time) * self.eras_code.k
         self.beta = mpf(1) / self.poisson_rate
          
         logging.debug("Inverse of maximum rate is %e" % self.beta)

         # Initialize the Poisson process via uniformization
         self.pp = PoissonProcess(self.beta, self.mission_time)
         self.fb_prob = self.is_parms.fb_prob
          
         # Initialize the state of the system
         self.state = State([i for i in range(self.num_components)])

         # Likelihood ratio
         self.lr = mpf(1)
     
    ##
    # Reset the simulator
    #
    def reset(self):
         # Generate next failure event

          # Reset clocks and state
          for component in self.components:
                component.init_clock(0)
                component.init_state()
          
          # Reset system state
          self.state = State(self.components)

          # Reset LR
          self.lr = mpf(1) 

    ##
    # Get the current event rate
    #
    def get_event_rate(self, poisson_rate):
          
        event_rate = mpf(0)
       
        ##
        # If repairs are ongoing, and we are doing BFB, then 
        # event rate is the same as w.o. BFB
        #
        if self.state.get_sys_state() == self.state.CURR_STATE_DEGRADED:
            for component in self.components:
                event_rate += component.inst_rate_sum()

            return ((event_rate / poisson_rate), (event_rate / poisson_rate))
    
    ##
    # Balanced failure biasing failure rate
    #
    def bfb_fail_rate(self, total_event_rate):
        return self.fb_prob * total_event_rate     
     
    ##
    # Balanced failure biasing repair rate
    #
    def bfb_repair_rate(self, total_event_rate):
        return (1-self.fb_prob) * total_event_rate

    ##
    # Get the next event
    #  
    def get_next_event(self, curr_time):
         
         # If not in a failed state, then draw for next device failure
         if self.state.get_sys_state() == self.state.CURR_STATE_OK:
              avail_comps = self.state.get_avail_components()
              #clock_val = self.components[0].read_clock()
              #next_event_time = self.components[0].component_fail_distr.draw_truncated(clock_val) - clock_val
              next_event_comp = 0
              next_event_time = self.components[0].component_fail_distr.draw() + curr_time
              
              for comp_idx in range(1, len(self.components)):
                    #clock_val = self.components[comp_idx].read_clock()
                    #event_time = self.components[comp_idx].component_fail_distr.draw_truncated(clock_val) - clock_val
                    event_time = self.components[comp_idx].component_fail_distr.draw() + curr_time
                    
                    if event_time < next_event_time:
                         next_event_time = event_time
                         next_event_comp = comp_idx
                         
              # Update internal component state
              self.components[next_event_comp].fail_component()
              
              return (next_event_time, Component.EVENT_COMP_FAIL, next_event_comp)
                         
         elif self.state.get_sys_state() == self.state.CURR_STATE_DEGRADED:
             draw = random.uniform()
        
             (act_event_rate, is_event_rate) = self.get_event_rate(self.poisson_rate)
             
             next_event_time = self.pp.draw() + curr_time

             # Determine if it is a "real" event or "pseudo" event
             if draw > is_event_rate:
                 # It is a pseudo event 

                 # Update the LR to deal with this
                 self.lr *= 1 # The event rate is the same for IS vs. "real world"

                 # Return nothing, since we are staying in the current state
                 return (next_event_time, None, None)

             else:
                fail_rate = self.bfb_fail_rate(act_event_rate)
                draw = random.uniform()

                # Failure
                if draw <= (fail_rate / act_event_rate):
                    avail_comps = self.state.get_avail_components()
                    comp_id = avail_comps[random.randint(0, len(avail_comps)-1)]
                    event_type = Component.EVENT_COMP_FAIL
                    
                    total_fail_rate = 0
                    for component_id in self.state.get_avail_components():
                        total_fail_rate += self.components[component_id].curr_component_fail_rate()

                    # Update the LR to deal with this
                    self.lr *= (self.components[comp_id].curr_component_fail_rate()) / ((fail_rate)/len(avail_comps))
                    
                    # Update internal component state
                    self.components[comp_id].fail_component()
                # Repair a component
                else:
                    repair_rate = self.bfb_repair_rate(act_event_rate)
                    draw = random.uniform()
                    total_repair_rate = mpf(0)
                    repair_rates = []
                    
                    for component_id in self.state.get_failed_components():
                        total_repair_rate += self.components[component_id].curr_component_repair_rate()

                    repair_prob_sum = mpf(0)
                    for component_id in self.state.get_failed_components():
                        repair_prob_sum += (self.components[component_id].curr_component_repair_rate() / total_repair_rate)

                        if draw < repair_prob_sum:
                            comp_id = component_id
                            break
                    
                    # Update the LR to deal with this
                    self.lr *= ((self.components[comp_id].curr_component_repair_rate())/((repair_rate) / len(self.state.get_failed_components())))
                  

                    # Update internal component state
                    self.components[comp_id].repair_component()
                    event_type = Component.EVENT_COMP_REPAIR
                    
                return (next_event_time, event_type, comp_id)
        
    ##
    # Run an iteration of the simulator
    #
    def run_iteration(self):
        self.reset()
        curr_time = 0
         
        while 1:

            (event_time, event_type, component_id)  = self.get_next_event(curr_time)
            curr_time = event_time
            
            if event_time > self.mission_time:
                break
            
            # Update all component clocks
            for component in self.components:
                component.update_clock(event_time)

            if event_type is not None:
                logging.debug("TIME %s, EVT TYPE: %s, COMP ID: %s, CURR COMP FAILED: %s\n" % (event_time, event_type, component_id, self.state.get_num_component_fail()))

            # Update the current state of the simulation
            state_summary = self.state.update_state(event_type, (component_id,))

            # Put check in to see if #erasures == HD-1 and see if there are any sector failures

             # Check to see if #erasures  >= minimum disk fault tolerance of erasure code
            if self.eras_code.min_disk_failures <= self.state.get_num_component_fail():

                failed_comps = self.state.get_failed_components()
                
                logging.debug("Components failed : %s" % failed_comps)

                # Check to see if we are in the "failed" state
                if self.eras_code.is_failure(failed_comps) is True:
                    logging.debug("LR : %e" % self.lr)
                    return self.lr

        return 0