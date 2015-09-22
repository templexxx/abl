##
# This module contains a container class for simulating 
# a reliable system.
#
# Kevin Greenan (kmgreen@cs.ucsc.edu)
#

from smp_data_structures import *
from state import *
from erasure_code import *
import logging
from mpmath import *
from sim_analysis_functions import *
from numpy import random
from poisson_process import *
from simulation import Simulation, ISParms
from unif_bfb_gen_repair import UniformizationBFBOpt
from bfb_optimization import BFBOpt
from regular_simulation import RegularSimulation
from sector_fail_model import *

#
# Imports for profiling
#
import cProfile
import pstats

##
# Print debug messages
#
#logging.basicConfig(level=logging.DEBUG)
#logging.basicConfig(level=logging.INFO)

global profile
profile = False

class Simulate:
    def __init__(self, code_desc, num_components, mission_time, is_type, is_parms, sector_fail_model, component_fail_dists, component_repair_dists, fail_check_type=None, critical_region_check=False):
        
        # Importance sampling type
        self.is_type = is_type
        
        if self.is_type == Simulation.IS_UNIF_BFB_NO_FORCING_OPT:
            self.sim = UniformizationBFBOpt(code_desc, num_components, mission_time, is_parms, sector_fail_model, component_fail_dists, component_repair_dists, fail_check_type,critical_region_check)
            self.sim.init()
        elif self.is_type == Simulation.IS_BFB_NO_FORCING_OPT:
            self.sim = BFBOpt(code_desc, num_components, mission_time, is_parms, sector_fail_model, component_fail_dists, component_repair_dists, fail_check_type,critical_region_check)
            self.sim.init()
        elif self.is_type == Simulation.REGULAR:
            self.sim = RegularSimulation(code_desc, num_components, mission_time, None, sector_fail_model, component_fail_dists, component_repair_dists, fail_check_type,critical_region_check)
            self.sim.init()
            
    def run_simulation(self, num_iterations=10):
        run_samples = []
        run_patterns = []
        if self.sim.sector_failure_model is not None:
            num_sectors_per_disk = self.sim.sector_failure_model.total_num_sectors
        else:
            num_sectors_per_disk = 585937500

            
        bytes_per_sector = 512
        avg_bytes_lost = 0
        distinct_patterns = {}
        pattern_probs = {}
        
        for i in range(num_iterations):
            #print "Iteration: ", i
            (sample, pattern, critical_region) = self.sim.run_iteration()
            if not distinct_patterns.has_key(pattern):
                distinct_patterns[pattern] = 0
                pattern_probs[pattern] = 0
            distinct_patterns[pattern] += 1
            pattern_probs[pattern] += sample
            
            if sample != 0:
                (num_disks, num_sectors) = eval(pattern)
                if num_sectors == 0:
                    avg_bytes_lost += (critical_region*sample)
                else:
                    avg_bytes_lost += 1
                    
            run_samples.append(sample)
            run_patterns.append(pattern)
    
        for pattern in pattern_probs.keys():
            pattern_probs[pattern] /= num_iterations
            
        avg_bytes_lost = ((avg_bytes_lost * bytes_per_sector)/ num_iterations) 
        return (run_samples, avg_bytes_lost, distinct_patterns, pattern_probs)

def test():
    num_components = 8
    mission_time = 87600
    iterations = 1000
    component_fail_dists = [Weibull(shape=1.12, scale=461386.) for i in range(num_components)]
    component_repair_dists = [Weibull(shape=2.0, scale=12., location=6.0) for i in range(num_components)]
    num_sectors_per_disk = 2147483648

    sector_fail_model = DeterministicScrubSectorFailModel(num_sectors_per_disk, num_sectors_per_disk, 168, 4.096e-11, 0.0045, 1)

    simulation = Simulate("7_1_mds", num_components, mission_time, Simulation.REGULAR, ISParms(forcing_prob=0.8, fb_prob=0.3), sector_fail_model, component_fail_dists, component_repair_dists, ErasureCode.CHECK_FTV,critical_region_check=False)

    (run_samples, avg_bytes_lost, distinct_patterns, pattern_probs) = simulation.run_simulation(iterations)

    samples = Samples(run_samples)
        
    mean = samples.calcMean()

    (low_ci, high_ci) = samples.calcConfInterval("0.90")

    relative_error = 100*samples.calcRE()

    print "Estimated reliability: %e +/- %f Percent, CI (%e,%e), num_zeroes: %d" % (mean, relative_error, low_ci, high_ci, samples.get_num_zeroes())
    print "Average bytes lost: ", avg_bytes_lost
    print "DL patterns: ", distinct_patterns
    print "Pattern probabilities: ", pattern_probs
    
if __name__ == "__main__":
    if profile is True:
        cProfile.run('test()', 'profile.out')
        p = pstats.Stats('profile.out')
        p.sort_stats('time')
        p.print_stats()
    else:
        test()
        


