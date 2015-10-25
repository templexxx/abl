import os
import sys
sys.path.append(os.getenv('PWD') + '/lib')

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
import getopt

class Simulate:
    def __init__(self, code_desc, num_components, mission_time, is_type, 
                 is_parms, sector_fail_model, component_fail_dists, 
                 component_repair_dists, fail_check_type=None,
                 critical_region_check=True):

        # Importance sampling type
        self.is_type = is_type

        if self.is_type == Simulation.IS_UNIF_BFB_NO_FORCING_OPT:
            self.sim = UniformizationBFBOpt(code_desc, num_components, mission_time, 
                                            is_parms, sector_fail_model, component_fail_dists, 
                                            component_repair_dists, fail_check_type,
                                            critical_region_check)
            self.sim.init()

        elif self.is_type == Simulation.IS_BFB_NO_FORCING_OPT:
            self.sim = BFBOpt(code_desc, num_components, mission_time, is_parms, sector_fail_model, 
                              component_fail_dists, component_repair_dists, fail_check_type,
                              critical_region_check)
            self.sim.init()

        elif self.is_type == Simulation.REGULAR:
            self.sim = RegularSimulation(code_desc, num_components, mission_time, None, sector_fail_model, 
                                         component_fail_dists, component_repair_dists, fail_check_type,
                                         critical_region_check)
            self.sim.init()

    def run_simulation(self, num_iterations=10):
        run_samples = []
        run_patterns = []
        if self.sim.sector_failure_model is not None:
            num_sectors_per_disk = self.sim.sector_failure_model.total_num_sectors
        else:
            num_sectors_per_disk = 1000000000.0

        bytes_per_sector = 4096.0
        avg_bytes_lost = 0
        distinct_patterns = {}
        pattern_probs = {}

        for i in range(num_iterations):
            (sample, pattern, critical_region) = self.sim.run_iteration()
            if not distinct_patterns.has_key(pattern):
                distinct_patterns[pattern] = 0                
                pattern_probs[pattern] = 0
            distinct_patterns[pattern] += 1
            pattern_probs[pattern] += sample

            if sample != 0:
                (num_disks, num_sectors) = eval(pattern)
                if self.sim.sector_failure_model is None:
                    avg_bytes_lost += num_sectors_per_disk
                elif num_sectors == 0:
                    avg_bytes_lost += (critical_region*sample)
                else:
                    avg_bytes_lost += 1

            run_samples.append(sample)
            run_patterns.append(pattern)

        for pattern in pattern_probs.keys():
            pattern_probs[pattern] /= num_iterations

        avg_bytes_lost = ((avg_bytes_lost * bytes_per_sector)/ num_iterations)
        return (run_samples, avg_bytes_lost, distinct_patterns, pattern_probs)

def usage(arg):
    print arg, ": -h [--help] -m <mission_time> [--mission_time <mission_time>]"
    print "-n <num_components> [--num_components <num_components>] -i <num_iterations> [--iterations <num_iterations>]"
    print "-C <code_file> [--code_file <code_file>]"
    print "-k <num_data_symbols>"
    print ""

    sys.exit(2)

def get_parms():

    sim_mode = None
    mission_time = 35040
    iterations = 10000
    fault_check = ErasureCode.CHECK_FTV
    critical_region_check = True
    code_file = None
    num_components = None
    use_is = False
    is_type = None
    is_forcing_prob = None
    is_fb_prob = None
    sector_failure_model = None
    component_fail_dist = None
    component_repair_dist = None
    bad_opt = None
    kt = None

    try:
        (opts, args) = getopt.getopt(sys.argv[1:], "hs:m:n:i:f:cC:S:F:R:k:", ["help",  "sim_mode", "mission_time",
                                                                             "num_components", "iterations", "fault_check", 
                                                                             "critical_check", "code_file", "sector_failure_model",
                                                                             "component_fail_dist", "component_repair_dist", "kt"])
    except:
        usage(sys.argv[0])
        print "getopts excepted"
        sys.exit(1)

    for o, a in opts:
        if o in ("-h", "--help"):
            print usage(sys.argv[0])
            sys.exit(0)
        if o in ("-F", "--component_fail_dist"):
            if len(eval(a)) == 1:
                component_fail_dist = Weibull(1, eval(a), 0)

            elif len(eval(a)) == 2:
                (shape, scale) = eval(a)
                component_fail_dist = Weibull(shape, scale, 0)

            elif len(eval(a)) == 3:
                (shape, scale, location) = eval(a)
                component_fail_dist = Weibull(shape, scale, location)

            else:
                bad_opt = o + " : " + a
                break
        
        elif o in ("-R", "--component_repair_dist"):
            if len(eval(a)) == 1:
                component_repair_dist = Weibull(1, eval(a), 0)

            elif len(eval(a)) == 2:
                (shape, scale) = eval(a)
                component_repair_dist = Weibull(shape, scale, 0)

            elif len(eval(a)) == 3:
                (shape, scale, location) = eval(a)
                component_repair_dist = Weibull(shape, scale, location)

            else:
                bad_opt = o + " : " + a
                break

        elif o in ("-S", "--sector_failure_model"):
            sec_fail_tup = eval(a)
            if len(sec_fail_tup) == 2: # total_num_sectors, sector_fail_prob
                (total_num_sectors, sector_fail_prob) = eval(a)
                total_num_sectors = mpf(total_num_sectors)
                sector_fail_prob = mpf(sector_fail_prob)
                sector_failure_model = BERSectorFailModel(total_num_sectors, None, None, sector_fail_prob, None, None)

            elif len(sec_fail_tup) == 6: # type, total_num_sectors, sector_fail_prob, sectors_per_region, scrub_interval, request_rate, write_ratio
                (type, total_num_sectors, sector_fail_prob, sectors_per_region, scrub_interval, request_rate, write_ratio) = eval(a)

                total_num_sectors = mpf(total_num_sectors)
                sector_fail_prob = mpf(sector_fail_prob)
                sectors_per_region = mpf(sectors_per_region)
                scrub_interval = mpf(scrub_interval)
                request_rate = mpf(request_rate)
                write_ratio = mpf(write_ratio)

                if type == "random":
                    sector_failure_model = NoScrubSectorFailModel(total_num_sectors, sectors_per_region, scrub_interval, sector_fail_prob, request_rate, write_ratio)
                elif type == "deterministic":
                    sector_failure_model = DeterministicScrubSectorFailModel(total_num_sectors, sectors_per_region, scrub_interval, sector_fail_prob, request_rate, write_ratio)
                else:
                    bad_opt = o + " : " + type
                    break
            else:
                bad_opt = o + " : num args must be 2 or 6"
                break

        elif o in ("-s", "--sim_mode"):
            if a == "reg":
                sim_mode = Simulation.REGULAR
            else: 
                (is_type, is_forcing_prob, is_fb_prob) = eval(a) 
                if is_type == "bfb":
                    sim_mode = Simulation.IS_BFB_NO_FORCING_OPT
                elif is_type == "unif":
                    sim_mode = Simulation.IS_UNIF_BFB_NO_FORCING_OPT
                else:
                    bad_opt = o + " : " + is_type
                    break

                is_forcing_prob = mpf(is_forcing_prob)
                is_fb_prob = mpf(is_fb_prob)

        elif o in ("-m", "--mission_time"):
            mission_time = float(a)

        elif o in ("-k", "--k"):
            kt = 3.7253 * int(a)

        elif o in ("-i", "--iterations"):
            iterations = int(a) 

        elif o in ("-f", "--fault_check"):
            if a == "ftv":
                fault_check = ErasureCode.CHECK_FTV
            elif a == "mel":
                fault_check = ErasureCode.CHECK_MEL
            elif a == "rank":
                fault_check = ErasureCode.CHECK_RANK
            elif a == "dscft":
                fault_check = ErasureCode.CHECK_DSCFT
            else:
                bad_opt = o
                break

        elif o in ("-c", "--critical_check"):
            critical_region_check = True
        
        elif o in ("-C", "--code_file"):
            code_file = a 

        elif o in ("-n", "--num_components"):
            num_components = int(a) 
                 

    if sim_mode is None:
        sim_mode = Simulation.IS_UNIF_BFB_NO_FORCING_OPT
        is_forcing_prob = mpf(0.5)
        is_fb_prob = mpf(0.3)

    if code_file is None:
        code_file = "rs_10_4"
        num_components = 14

    if kt is None:
        kt = 37.253

    if num_components is None:
        print "Must give number of components when a code_file is specified."
        sys.exit(2)

    if component_fail_dist is None:
        component_fail_dist = [Weibull(shape=1.12, scale= 281257.0) for i in range(num_components)]
    else:
        tmp_fail = []
        for i in range(num_components):
            tmp_fail.append(component_fail_dist)

        component_fail_dist = tmp_fail

    if component_repair_dist is None:
        component_repair_dist = [Weibull(shape=2.0, scale=120.0, location=72.0) for i in range(num_components)]
    else:
        tmp_repair = []
        for i in range(num_components):
            tmp_repair.append(component_repair_dist)

        component_repair_dist = tmp_repair

    if sector_failure_model is None:
        total_num_sectors = mpf(1000000000)
        sfp = 3.2768e-10
        sector_fail_prob = mpf(sfp)
        sector_failure_model = BERSectorFailModel(total_num_sectors, None, None, sector_fail_prob, None, None)



    return (sim_mode, mission_time, iterations, fault_check, critical_region_check, 
            code_file, num_components, ISParms(forcing_prob=is_forcing_prob, fb_prob=is_fb_prob), 
            sector_failure_model, component_fail_dist, component_repair_dist, kt)
            
def do_it():

    (sim_mode, mission_time, iterations, fault_check, critical_region_check, 
     code_file, num_components, is_parms, sector_failure_model, 
     component_fail_dists, component_repair_dists, kt) = get_parms()

    simulation = Simulate(code_file, num_components, mission_time, sim_mode, is_parms, 
                          sector_failure_model, component_fail_dists, component_repair_dists, 
                          fault_check, critical_region_check)

    (run_samples, avg_bytes_lost, distinct_patterns, pattern_probs) = simulation.run_simulation(iterations)

    samples = Samples(run_samples)

    mean = samples.calcMean()

    (low_ci, high_ci) = samples.calcConfInterval("0.90")

    relative_error = 100*samples.calcRE()

    abl = avg_bytes_lost / kt


    print "\n*******************\n"
    print "Average bytes lost per usable TB: %.5f" % abl
    print "\n*******************\n"
    print "DL patterns: (disks, sectors): num instances"
    for key in distinct_patterns.keys():
        print key, ": ", distinct_patterns[key]
    print "\n*******************\n"


if __name__ == "__main__":
    do_it()


