import mpmath
from mpmath import mpf
from mpmath import mp
from numpy import *
from numpy.random import *
from scipy import integrate

mp.dps += 12
mp.prec += 75

class MarkovChain:
	REGULAR = "homogeneous"
	IS_SFB = "simple_failure_biasing"
	IS_BFB = "balanced_failure_biasing"
	IS_FORCE = "forcing"
	FAILURE="failure"
	REPAIR="repair"
	BIAS_PROB=mpf(0.5)

	def __init__(self, num_states, transitions, tr_type_matrix, type, bias_prob=0):
		# Transition matrix 
		self.transition_matrix = []

		self.bias_prob = MarkovChain.BIAS_PROB

		# Embedded Discrete-time chain
		self.jump_process = []

		# Invertable matrix (matrix without the absorbing state)
		# This matrix is used to compute the MTTDL 
		self.mttdl_matrix = []

		self.num_states = num_states

		self.poss_transitions = {}

		self.holding_times = {}

		self.num_states = num_states

		self.tr_type_matrix = tr_type_matrix

		if type == self.REGULAR:
			self.reg_transition_matrix(transitions)
			self.reg_dtmc()
		elif type == self.IS_SFB:
			self.reg_transition_matrix(transitions)
			self.is_sfb_dtmc()
		elif type == self.IS_BFB:
			self.reg_transition_matrix(transitions)
			self.is_bfb_dtmc()
		elif type == self.IS_FORCE:
			self.reg_transition_matrix(transitions)
			self.is_force_dtmc()

		self.type = type
	
	def print_transition_matrix(self):
		print self.transition_matrix
	
	def print_jump_process(self):
		print self.jump_process
	
	def print_poss_transitions(self):
		print self.poss_transitions
	
	def print_mttdl_matrix(self):
		print self.mttdl_matrix

	def get_mttdl(self):
		if len(self.mttdl_matrix ) == 0:
			return None

		m_inv = mat(linalg.inv(transpose(mat(self.mttdl_matrix))))

		e = mat([-1 for i in range(len(self.mttdl_matrix))])

		init_vec = [[1]]

		for i in range(len(self.mttdl_matrix)-1):
			init_vec.append([0])

		init_vec = mat(init_vec)

		mttdl = (e * m_inv * init_vec) 

		return mttdl.tolist()[0][0]

	def estimate_prob_of_failure(self, t=0):
		mttdl = self.get_mttdl()
		return 1-mpmath.exp(mpf(-t)/mpf(mttdl))

	def prob_of_failure(self, t=0):

		# Vector containing initial probabilities 
		x0 = [0 for i in range(self.num_states)]
		x0[0] = 1

		# Build up ODEs for dp[i] / dt
		str="["

		for i in range(len(self.transition_matrix[0])):
			first_flg = False
			for j in range(len(self.transition_matrix)):
				if self.transition_matrix[j][i] != 0:
					if first_flg is False:
						str +=	"p[%d]*(%s)" % (j, self.transition_matrix[j][i])
						first_flg = True
					else:
						str +=	"+p[%d]*(%s)" % (j, self.transition_matrix[j][i])
			if i < len(self.transition_matrix[0]) - 1:
				str += ", "

		str += "]"

		x = integrate.odeint(lambda p, time, str: eval(str), x0, t, args=(str,))

		return abs(x)
		
	def is_sfb_dtmc(self):
		sum_fail_probs = []
		sum_repair_probs = []

		for i in range(self.num_states-1):
			self.jump_process.append([0 for k in range(self.num_states)])
			sum_fail_probs.append(mpf(0))
			sum_repair_probs.append(mpf(0))

			for j in range(self.num_states):
				if i != j and self.transition_matrix[i][j] > mpf(0):
					if self.tr_type_matrix[i][j] == self.FAILURE:
						sum_fail_probs[i] += (self.transition_matrix[i][j] / (-self.transition_matrix[i][i]))	
					elif self.tr_type_matrix[i][j] == self.REPAIR:
						sum_repair_probs[i] += (self.transition_matrix[i][j] / (-self.transition_matrix[i][i]))	

					self.jump_process[i][j] = self.transition_matrix[i][j] / (-self.transition_matrix[i][i])	
				else:
					self.jump_process[i][j] = mpf(0)

		for i in range(1, self.num_states-1):
			for j in range(self.num_states):
				if  i != j and self.jump_process[i][j] > 0:
					if sum_fail_probs[i] == 0 or sum_repair_probs[i] == 0:
						continue
					if self.tr_type_matrix[i][j] == self.FAILURE:
						self.jump_process[i][j] = self.bias_prob * (self.jump_process[i][j] / sum_fail_probs[i])
					elif self.tr_type_matrix[i][j] == self.REPAIR:
						self.jump_process[i][j] = (1-self.bias_prob)*(self.jump_process[i][j] / sum_repair_probs[i])
		
		# Fill in table of possible transitions and jump probs
		for i in range(self.num_states-1):
			self.poss_transitions[i] = []
			for j in range(self.num_states):
				if i != j and self.transition_matrix[i][j] > 0:
					# Append (to_state, probability)
					self.poss_transitions[i].append((j, self.jump_process[i][j]))
		
		return None

	def is_bfb_dtmc(self):
		sum_fail_events = []
		sum_repair_probs = []

		for i in range(self.num_states-1):
			self.jump_process.append([0 for k in range(self.num_states)])
			sum_fail_events.append(0)
			sum_repair_probs.append(mpf(0))

			for j in range(self.num_states):
				if i != j and self.transition_matrix[i][j] > mpf(0):
					if self.tr_type_matrix[i][j] == self.FAILURE:
						sum_fail_events[i] += 1
					elif self.tr_type_matrix[i][j] == self.REPAIR:
						sum_repair_probs[i] += (self.transition_matrix[i][j] / (-self.transition_matrix[i][i]))	

					self.jump_process[i][j] = self.transition_matrix[i][j] / (-self.transition_matrix[i][i])	
				else:
					self.jump_process[i][j] = mpf(0)

		for i in range(self.num_states-1):
			for j in range(self.num_states):
				if  i != j and self.jump_process[i][j] > 0:
					if sum_fail_events[i] == 0 or sum_repair_probs[i] == 0:
						continue
					if self.tr_type_matrix[i][j] == self.FAILURE:
						self.jump_process[i][j] = self.bias_prob * (mpf(1) / sum_fail_events[i])
					elif self.tr_type_matrix[i][j] == self.REPAIR:
						self.jump_process[i][j] = (1-self.bias_prob)*(self.jump_process[i][j] / sum_repair_probs[i])
		
		# Fill in table of possible transitions and jump probs
		for i in range(self.num_states-1):
			self.poss_transitions[i] = []
			for j in range(self.num_states):
				if i != j and self.transition_matrix[i][j] > 0:
					# Append (to_state, probability)
					self.poss_transitions[i].append((j, self.jump_process[i][j]))
		
	def is_force_dtmc(self):
		return None

	def reg_dtmc(self):
		# Fill in the embedded DTMC
		# The last state is the lone absorbing (failed state), so no row is necessary
		for i in range(self.num_states-1):
			self.jump_process.append([0 for k in range(self.num_states)])
			
			# The last state is the lone absorbing (failed state), so no row is necessary
			for j in range(self.num_states):
				if i != j and self.transition_matrix[i][j] > mpf(0):
					self.jump_process[i][j] = self.transition_matrix[i][j] / (-self.transition_matrix[i][i])	
				else:
					self.jump_process[i][j] = mpf(0)
		
		# Fill in table of possible transitions and jump probs
		for i in range(self.num_states-1):
			self.poss_transitions[i] = []
			for j in range(self.num_states):
				if i != j and self.transition_matrix[i][j] > 0:
					# Append (to_state, probability)
					self.poss_transitions[i].append((j, self.jump_process[i][j]))

	def reg_transition_matrix(self, transitions):
		# Fill in the transition matrix
		# The last state is the lone absorbing (failed state), so no row is necessary
		for i in range(self.num_states-1):
			self.transition_matrix.append([0 for k in range(self.num_states)])
			self.mttdl_matrix.append([mpf(0) for k in range(self.num_states-1)])
			sum = mpf(0)

			for j in range(self.num_states):
				if j < self.num_states-1 and i != j:
					self.mttdl_matrix[i][j] = mpf(transitions[i][j])
				if i != j:
					self.transition_matrix[i][j] = mpf(transitions[i][j])
					sum += mpf(transitions[i][j])
			self.transition_matrix[i][i] = (sum * -1)
			self.mttdl_matrix[i][i] = mpf(sum * -1)
		

		# Fill in holding times for each state
		for i in range(self.num_states-1):
			self.holding_times[i] = -(self.transition_matrix[i][i])

	def get_next_state(self, current_state):
		transitions = self.poss_transitions[current_state]
		draw = uniform(0, 1)
		curr_prob = 0
		
		for trans in transitions:
			if draw >= curr_prob and draw < (curr_prob+trans[1]):
				return trans[0]
			curr_prob += trans[1]
		
		return current_state
