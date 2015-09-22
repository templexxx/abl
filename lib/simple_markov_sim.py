##
# Importance sampling hack...  Everything is exponentially distributed, 
# though this is not always the case...  
#

from numpy import *
from numpy.random import *
from mpmath import *
mp.dps += 12
mp.prec += 75
from markov_chain import *
from sim_analysis_functions import Samples
from build_transition_matrix import build_markov_chain
import random

class SimpleMCSim:
    def __init__(self, model_def, mission_time, num_iterations, is_type, fail_bias_prob=None):
        if is_type == MarkovChain.IS_BFB or is_type == MarkovChain.IS_SFB:
            self.mc = build_markov_chain(model_def, is_type)
            self.orig_mc = build_markov_chain(model_def, MarkovChain.REGULAR)
            
        else:
            self.mc = build_markov_chain(model_def, MarkovChain.REGULAR)
            self.orig_mc = None

        self.mission_time = mission_time
        self.num_iterations = num_iterations
        self.model_def = model_def

    def run_iteration(self):
        curr_time = 0
        curr_state = 0
        lr = 1

        while curr_time < self.mission_time:
            curr_time += random.expovariate(self.mc.holding_times[curr_state])
            old_state = curr_state
            curr_state = self.mc.get_next_state(curr_state)
            if not self.orig_mc is None:
                lr *= (self.orig_mc.jump_process[old_state][curr_state] / self.mc.jump_process[old_state][curr_state])
            if curr_state == (self.mc.num_states - 1):
                break

        if curr_time < self.mission_time:
            return lr
        else:
            return 0

    def run_simulation(self):
        samples = []

        for i in range(self.num_iterations):
            samples.append(self.run_iteration())

        return samples
        



def main():
    sims = []
    
    sims.append(SimpleMCSim('models/flat_5_3.disk.model', 87600, 10000, MarkovChain.IS_BFB, 0.99))
    sims.append(SimpleMCSim('models/flat_5_3.disk.model', 87600, 10000, MarkovChain.IS_SFB, 0.99))
    sims.append(SimpleMCSim('models/2DFT.disk.model', 87600, 10000, MarkovChain.IS_BFB, 0.99))
    sims.append(SimpleMCSim('models/2DFT.disk.model', 87600, 10000, MarkovChain.IS_SFB, 0.99))
    sims.append(SimpleMCSim('models/3DFT.disk.model', 87600, 10000, MarkovChain.IS_BFB, 0.99))
    sims.append(SimpleMCSim('models/3DFT.disk.model', 87600, 10000, MarkovChain.IS_SFB, 0.99))
    sims.append(SimpleMCSim('models/4DFT.disk.model', 87600, 10000, MarkovChain.IS_BFB, 0.99))
    sims.append(SimpleMCSim('models/4DFT.disk.model', 87600, 10000, MarkovChain.IS_SFB, 0.99))
    sims.append(SimpleMCSim('models/5DFT.disk.model', 87600, 10000, MarkovChain.IS_BFB, 0.99))
    sims.append(SimpleMCSim('models/5DFT.disk.model', 87600, 10000, MarkovChain.IS_SFB, 0.99))
    
    for sim in sims:
    
        print "Running %s under %s" % (sim.model_def, sim.mc.type)
    
        num_zeros = 0
        samples = []
        
        for i in range(sim.num_iterations):
            samples.append(sim.run_iteration())
            if samples[i] == 0:
                num_zeros += 1
    
        print "NUM ZEROS", num_zeros
        s = Samples(samples)
    
        mean = s.calcMean()
        confInt = s.calcConfInterval()
        re = s.calcRE()
        print "Num failures/mission time: %s +/- %f (%e,%e)" % (mean, re, confInt[0], confInt[1])
        
        if sim.orig_mc is None:
            exact_prob_failure = sim.mc.prob_of_failure(t=range(0, 87600, 1000))
            print "Prob. of failure determined by solving ODEs: ", exact_prob_failure[len(exact_prob_failure)-1][sim.mc.num_states-1]
            print "Est. prob. of failure determined by exp(-t/MTTDL):", sim.mc.estimate_prob_of_failure(t=sim.mission_time)
        else:
            exact_prob_failure = sim.orig_mc.prob_of_failure(t=range(0, 87600, 1000))
            print "Prob. of failure determined by solving ODEs: ", exact_prob_failure[len(exact_prob_failure)-1][sim.mc.num_states-1]
            print "Est. prob. of failure determined by exp(-t/MTTDL):", sim.orig_mc.estimate_prob_of_failure(t=sim.mission_time)
        print "\n***************\n"
    

if __name__ == "__main__":
    main()

