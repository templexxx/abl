from numpy import random

##
# A class that generates arrival times for a Poisson Process
# with rate parameter beta.
#
class PoissonProcess:

    ##
    # Constructor for the Poisson Process
    #
    # @param beta: rate parameter for process
    # @param mission_time: mission time of the process
    #
    def __init__(self, beta=1, mission_time=100):
        self.beta = beta
        self.mission_time = mission_time

    ##
    # Generate and return arrivals 
    #
    # @return list of arrival times 
    #
    def generate_arrivals(self):
        curr_time = 0
        arrivals = []

        while 1:
            draw = random.exponential(self.beta)
            curr_time += draw
            if curr_time > self.mission_time:
                break
            arrivals.append(curr_time)

        return arrivals
    
    ##
    # Draw the next arrival
    #
    def draw(self):
        return random.exponential(self.beta)