from simulation import *

class RegularSimulation(Simulation):
    
    ##
    # __init__() from Simulation 
    #
    
    ##
    # Initialize the simulation
    #
    def init(self):
        # Initialize the state of the system
        self.state = State([i for i in range(self.num_components)])
        
        # Set up structure for component repairs
        self.component_repair_start = [0 for i in range(len(self.components))]
        self.component_repairs = [0 for i in range(len(self.components))]
        
         # Reset LR
        self.lr = mpf(1) 
          
    ##
    # Reset the simulation
    #
    def reset(self):
         # Reset system state
        self.state = State(self.components)
        
        # Draw first component fail times
        self.component_failures = [0 for i in range(len(self.components))]
        
        i=0
        for component in self.components:
            self.component_failures[i] = component.component_fail_distr.draw()
            i+=1
        
        # Set up structure for component repairs
        self.component_repairs = [0 for i in range(len(self.components))]
        
        self.sim_time = 0
        
        for component in self.components:
            component.init_clock(0)
            component.init_state()
        
         # Reset LR
        self.lr = mpf(1) 
    
    ##
    # Set new component failure time for component comp_idx
    #
    def set_comp_fail(self, comp_idx, curr_time):
        self.component_failures[comp_idx] = self.components[comp_idx].component_fail_distr.draw() + curr_time
       
    ##
    # Set new component repair time for component comp_idx
    #
    def set_comp_repair(self, comp_idx, curr_time):
        self.component_repairs[comp_idx] = self.components[comp_idx].component_repair_distr.draw() + curr_time 
        self.component_repair_start[comp_idx] = curr_time 
        
    ##
    # Get the next failure event
    #
    def get_next_failure(self, avail_comps):
        
        if len(avail_comps) == 0:
            return (None, None)
        
        comp_idx = avail_comps[0]
        fail_time = self.component_failures[avail_comps[0]]
        
        for i in avail_comps[1:]:
            if self.component_failures[i] < fail_time and self.component_failures[i] > 0:
                comp_idx = i
                fail_time = self.component_failures[i]
                
        return (fail_time, comp_idx)
    
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
        failed_comps = self.state.get_failed_components()
        avail_comps = self.state.get_avail_components()
        
        if len(failed_comps) == 0:
            (fail_time, comp_idx) = self.get_next_failure(avail_comps)
            
            ##
            # Update internal component state
            #
            self.components[comp_idx].fail_component(fail_time)
            
            ##
            # Draw component repair time
            #
            self.set_comp_repair(comp_idx, fail_time)
            
            return (fail_time, Component.EVENT_COMP_FAIL, comp_idx)
        
        else:
            (fail_time, fail_comp_idx) = self.get_next_failure(avail_comps)
            (repair_time, repair_comp_idx) = self.get_next_repair(failed_comps)
            
            if fail_time < repair_time:
                ##
                # Update internal component state
                #
                self.components[fail_comp_idx].fail_component(fail_time)
                
                ##
                # Draw component repair time
                #
                self.set_comp_repair(fail_comp_idx, fail_time)
                
                return (fail_time, Component.EVENT_COMP_FAIL, fail_comp_idx)
            else:
                ##
                # Update internal component state
                #
                self.components[repair_comp_idx].repair_component()
                
                ##
                # Draw component fail time
                #
                self.set_comp_fail(repair_comp_idx, repair_time)
                
                return (repair_time, Component.EVENT_COMP_REPAIR, repair_comp_idx)
        
    
    