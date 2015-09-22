from simulation import *
from numpy import random

class BFBOpt(Simulation):
     
    ##
    # __init__() from Simulation 
    #
     
    ##
    # Initialize the simulation
    #
    def init(self):
        
         self.fb_prob = self.is_parms.fb_prob
          
         # Initialize the state of the system
         self.state = State([i for i in range(self.num_components)])

         # Likelihood ratio
         self.lr = mpf(1)
         
         self.inv_transform_variates = InverseTransformHomogeneousFailRepairRates(self.components)
     
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
    def get_event_rate(self):
          
        event_rate = mpf(0)
       
        for component in self.components:
            event_rate += component.inst_rate_sum()

        return event_rate
    
    ##
    # Get repair rate
    #
    def get_repair_rate(self):
        failed_components = self.state.get_failed_components()
        repair_rate = 0
        
        for comp in failed_components:
            repair_rate += self.components[comp].curr_component_repair_rate()
        
        return repair_rate
    ##
    # Get failure rate
    #
    def get_fail_rate(self):
        avail_components = self.state.get_avail_components()
        fail_rate = 0
        
        for comp in avail_components:
            fail_rate += self.components[comp].curr_component_fail_rate()
        
        return fail_rate
    
    def get_time_to_repair(self):
        failed_components = self.state.get_failed_components()
        pending_repairs = []
        
        for comp in failed_components:
            pending_repairs.append(self.components[comp].component_repair_distr.draw_inverse_transform(self.components[comp].repair_clock))
                
        return min(pending_repairs)
    
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
    # Get the next event - Make sure times are correctly being updated!
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
              
              # Update all component clocks
              for component in self.components:
                  component.update_clock(next_event_time)
              
              return (next_event_time, Component.EVENT_COMP_FAIL, next_event_comp)
                
         elif self.state.get_sys_state() == self.state.CURR_STATE_DEGRADED:
             
             next_event_time = self.inv_transform_variates.draw_waiting_time(self.state.get_avail_components(), self.state.get_failed_components()) + curr_time
             
             draw = random.uniform()

             # Update all component clocks
             for component in self.components:
                 component.update_clock(next_event_time)
                     
             event_rate = self.get_event_rate()
             
             fail_rate = self.bfb_fail_rate(event_rate)
             
             # Failure
             if draw <= (fail_rate / event_rate):
                 avail_comps = self.state.get_avail_components()
                 comp_id = avail_comps[random.randint(0, len(avail_comps)-1)]
                 event_type = Component.EVENT_COMP_FAIL
                 
                 total_fail_rate = 0
                 for component_id in self.state.get_avail_components():
                     total_fail_rate += self.components[component_id].curr_component_fail_rate()

                 # Update the LR to deal with this
                 self.lr *= (self.components[comp_id].curr_component_fail_rate()/event_rate) / (self.fb_prob/len(avail_comps))
                 
                 # Update internal component state
                 self.components[comp_id].fail_component(next_event_time)
                
             # Repair a component
             else:
                 repair_rate = self.bfb_repair_rate(event_rate)
                 draw = random.uniform()
                 total_repair_rate = self.get_repair_rate()
                 repair_rates = []
                 
                 repair_prob_sum = mpf(0)
                 for component_id in self.state.get_failed_components():
                     repair_prob_sum += (self.components[component_id].curr_component_repair_rate() / total_repair_rate)

                     if draw < repair_prob_sum:
                         comp_id = component_id
                         break
                    
                 # Update the LR to deal with this
                 self.lr *= (self.components[comp_id].curr_component_repair_rate()/event_rate)/((1-self.fb_prob) * (self.components[comp_id].curr_component_repair_rate()/total_repair_rate))
                  
                 
                 # Update internal component state
                 self.components[comp_id].repair_component()
                 event_type = Component.EVENT_COMP_REPAIR
             return (next_event_time, event_type, comp_id)
        