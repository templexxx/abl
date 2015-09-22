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
         
         #for component in self.components:
         #     dist = component.component_repair_distr     
         #self.poisson_rate = self.components[0].component_repair_distr.get_max_hazard_rate(self.components[0].component_repair_distr.scale*2) * (self.eras_code.min_disk_failures-1) * 2
         self.poisson_rate = self.components[0].component_repair_distr.get_max_hazard_rate(self.components[0].component_repair_distr.scale*3) * 2
         
         self.fb_prob = self.is_parms.fb_prob
         
         # Set up structure for component repairs
         self.component_repairs = [0 for i in range(len(self.components))]
         
         # Set up structure for component repairs
         self.component_repair_start = [0 for i in range(len(self.components))]
          
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
    # Get failure rate
    #
    def get_fail_rate(self):
        fail_rate = mpf(0)
        for component in self.components:
                fail_rate += component.curr_component_fail_rate()

        return fail_rate
    ##
    # Set new component repair time for component comp_idx
    #
    def set_comp_repair(self, comp_idx, curr_time):
        self.component_repairs[comp_idx] = self.components[comp_idx].component_repair_distr.draw() + curr_time
        self.component_repair_start[comp_idx] = curr_time 
         
   ##
    # Get the next repair event
    #
    def get_next_repair(self, failed_comps):
        
        if len(failed_comps) == 0:
            return (None, None)
        
        comp_idx = failed_comps[0]
        repair_time = self.component_repairs[failed_comps[0]]
        
        for i in failed_comps[1:]:
            if self.component_repairs[i] < repair_time and self.component_repairs[i] > 0:
                comp_idx = i
                repair_time = self.component_repairs[i]
        
        return (repair_time, comp_idx)

    ##
    # Get the next event
    #  
    def get_next_event(self, curr_time):
         
         # If not in a failed state, then draw for next device failure
         if self.state.get_sys_state() == self.state.CURR_STATE_OK:
              avail_comps = self.state.get_avail_components()
              next_event_comp = 0
              clock_val = self.components[0].read_clock()
              next_event_time = self.components[0].component_fail_distr.draw_inverse_transform(clock_val) + curr_time
              
              for comp_idx in range(1, len(self.components)):
                    clock_val = self.components[comp_idx].read_clock()
                    event_time = self.components[comp_idx].component_fail_distr.draw_inverse_transform(clock_val) + curr_time
                    
                    if event_time < next_event_time:
                         next_event_time = event_time
                         next_event_comp = comp_idx
                         
              # Update internal component state
              self.components[next_event_comp].fail_component(next_event_time)
              
              # Schedule repair time
              self.set_comp_repair(next_event_comp, next_event_time)
              
              return (next_event_time, Component.EVENT_COMP_FAIL, next_event_comp)
                         
         elif self.state.get_sys_state() == self.state.CURR_STATE_DEGRADED:

             next_event_time = random.exponential(1/self.poisson_rate) + curr_time
             
             (repair_time, comp_idx) = self.get_next_repair(self.state.get_failed_components())
             
             if repair_time < next_event_time:
                 # Update internal component state
                self.components[comp_idx].repair_component()
                event_type = Component.EVENT_COMP_REPAIR
                    
                return (repair_time, Component.EVENT_COMP_REPAIR, comp_idx)
             
             # Update all component clocks
             for component in self.components:
                component.update_clock(next_event_time)
             
             draw = random.uniform()
             # Determine if it is a "real" event or "pseudo" event
             if draw > self.fb_prob:
                 # It is a pseudo event 
                 # Update the LR to deal with this
                 self.lr *=  ((1 - (self.get_fail_rate()/self.poisson_rate)) / (1 - self.fb_prob))

                 # Return nothing, since we are staying in the current state
                 return (next_event_time, None, None)

             else:
                 
                avail_comps = self.state.get_avail_components()
                comp_idx = avail_comps[random.randint(0, len(avail_comps))]
                
                
                # Update the LR to deal with this
                self.lr *= ((self.components[comp_idx].curr_component_fail_rate()/self.poisson_rate) / (self.fb_prob /len(avail_comps)))
                     
                # Update internal component state
                self.components[comp_idx].fail_component(next_event_time)
                
                # Schedule repair time
                self.set_comp_repair(comp_idx, next_event_time)
                
                return (next_event_time, Component.EVENT_COMP_FAIL, comp_idx)
