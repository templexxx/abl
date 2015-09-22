##
# Base class for any simulation
#

from erasure_code import *
from smp_data_structures import *
from state import *
import logging
from mpmath import *
from numpy import random

##
# Just a container for importance sampling parameters 
#
class ISParms:
    def __init__(self, forcing_prob=0.8, fb_prob=0.3):
        self.forcing_prob = forcing_prob
        self.fb_prob = fb_prob

class Simulation:
    REGULAR="regular"
    IS_UNIF_BFB_NO_FORCING_OPT = "uniformization_balanced_failure_biasing_without_forcing_with_optimization"
    IS_BFB_NO_FORCING_OPT = "balanced_failure_biasing_without_forcing_with_optimization"
    FAILURE="failure"
    REPAIR="repair"
    
    def __init__(self,code_desc, num_components, mission_time, is_parms, sector_fail_model, component_fail_dists, component_repair_dists, fail_check_type=None, critical_region_check=True):
        
        # Load the erasure code to use within this simulation
        self.eras_code = ErasureCode(code_desc, fail_check_type)

        # Number of components in the system
        self.num_components = num_components

        # Mission time of the simulation
        self.mission_time = mission_time

        self.is_parms = is_parms

        # Component failure and repair distributions
        self.component_fail_dists = component_fail_dists
        self.component_repair_dists = component_repair_dists
        
        self.components = [Component(self.component_fail_dists[i], self.component_repair_dists[i]) for i in range(num_components)]

        self.component_repairs = None
        
        self.sector_failure_model = sector_fail_model
        
        self.component_repairs = None
        
        self.critical_region_flg = critical_region_check
        
        
    ##
    # Init the simulation 
    #
    def init(self):
        return None
    
    ##
    # Reset the simulator
    #
    def reset(self):
        return None 
    
    ##
    # Get the next event
    #  
    def get_next_event(self, curr_time):
        return None
    
    ##
    # Run an iteration of the simulator
    #
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


            if event_type is None or event_type == Component.EVENT_COMP_REPAIR:
                continue
            
            # Check to see if #erasures  >= minimum disk fault tolerance of erasure code
            if self.eras_code.min_disk_failures <= self.state.get_num_component_fail():

                failed_comps = self.state.get_failed_components()
                critical_region = 0
                
                logging.debug("Components failed : %s" % failed_comps)

                # Check to see if we are in the "failed" state
                if self.eras_code.is_failure(failed_comps) is True:
                    logging.debug("LR : %e" % self.lr)
                    if self.component_repairs is None and self.critical_region_flg is True:
                        critical_region = (mpf(1) / (1 << (self.state.get_num_component_fail()-1))) * self.sector_failure_model.total_num_sectors
                        #None
                    elif self.component_repairs is not None and self.critical_region_flg is True:
                        next_repair = self.component_repairs[failed_comps[0]]
                        next_repair_idx = failed_comps[0]
                        for i in failed_comps[1:]:
                            if next_repair > self.component_repairs[i]:
                                next_repair = self.component_repairs[i]
                                next_repair_idx = i
                        
                        if self.critical_region_flg is True and self.sector_failure_model is not None:
                            critical_region = ((next_repair - curr_time) / (next_repair - self.component_repair_start[next_repair_idx])) * self.sector_failure_model.total_num_sectors
                        else:
                            critical_region = 0
                        
                    return (self.lr, "(%d, %d)" % (self.state.get_num_component_fail(), 0), critical_region)
            
            # Put check in to see if #erasures == HD-1 and see if there are any sector failures
            if (self.eras_code.min_disk_failures-1) <= self.state.get_num_component_fail():
                if self.sector_failure_model is not None:
                    failed_comps = self.state.get_failed_components()
                    
                    # Compute the critical region 
                    if self.critical_region_flg is True and self.component_repairs is None:
                        max_time = self.components[failed_comps[0]].repair_clock
                        max_idx = failed_comps[0]
                        for i in failed_comps[1:]:
                            if self.components[i].repair_clock > max_time:
                                max_time = self.components[i].repair_clock
                                max_idx = i
                        critical_region = (mpf(1) / (1 << (self.state.get_num_component_fail()-1))) * self.sector_failure_model.total_num_sectors 
                      
                    elif self.critical_region_flg is True and self.component_repairs is not None:
                        next_repair = self.component_repairs[failed_comps[0]]
                        next_repair_idx = failed_comps[0]
                        for i in failed_comps[1:]:
                            if next_repair > self.component_repairs[i]:
                                next_repair = self.component_repairs[i]
                                next_repair_idx = i
                                
                        critical_region = ((next_repair - curr_time) / (next_repair - self.component_repair_start[next_repair_idx])) * self.sector_failure_model.total_num_sectors
                    else:
                        critical_region = self.sector_failure_model.total_num_sectors-1
                    
                    # Check to see if there is a sector failure on any disk
                    # For each sector failure, check for data loss
                    avail_components = self.state.get_avail_components()
                    sector_failures = [[] for i in range(self.num_components)]
                
                    for comp in avail_components:
                        draw = random.uniform()
                        if draw < self.sector_failure_model.prob_of_bad_sector():
                            sector_index = random.randint(0, self.sector_failure_model.total_num_sectors-1)
                            if sector_index < critical_region:
                                sector_failures[comp].append(sector_index)
                
                    if self.eras_code.is_failure(failed_comps, sector_failures) is True:
                        logging.debug("LR : %e" % self.lr)
                        return (self.lr, "(%d, %d)" % (self.state.get_num_component_fail(), 1), critical_region)

        return (0, "(0, 0)", 0)
