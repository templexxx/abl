import sys
from markov_chain import *
from mpmath import *
import scipy

##
# The model spec contains transitions, variables and values.  
# We reduce all of the variables to values by performing
# substitution.  Substitution is driven by the assignments 
# (var=val) in the specification.
#
class Assignment:
	def __init__(self, var, val):
		self.var = var
		self.val = val

	def __str__(self):
		return "%s=%s" % (self.var, self.val) 


def build_markov_chain(model_def, type, mult_factor={}):

	transitions = [{}]
	matrix = [[]]
	tr_type_matrix = [[]]
	num_states = 0
	assign = []
	symbolic = False

	model_def_fd = open(model_def, "r")
	
	##
	# Generate the transition matrix based on the states,
	# their transitions and transition rates (determined 
	# through substitution)
	#
	while 1:
		line = model_def_fd.readline().strip()
		if line == '':
			break
		elif line == "[num states]":
			line = model_def_fd.readline().strip()
			num_states = int(line)
			transitions = [{} for i in range(num_states-1)]
			matrix = [[0 for i in range(num_states)] for i in range(num_states-1)]
			tr_type_matrix = [[None for i in range(num_states)] for i in range(num_states-1)]
		elif line == "[symbolic]":
			symbolic = True
		elif line == "[assign]":
			row = {}
			while 1:
				line = model_def_fd.readline().strip()
				
				if line == "[END]":
					break
				
				replace = line.split('=')
				
				if mult_factor.has_key(replace[0]):
					replace[1] = "%s" % (eval(replace[1]) * mult_factor[replace[0]])
						
				assign.append(Assignment(replace[0], replace[1]))
			assign.sort(lambda x,y : len(y.var)-len(x.var));
	
			remove_list = []
			for i in range(len(assign)):
				for j in range(len(assign)):
					if i != j:
						temp = assign[j].val.replace(assign[i].var, assign[i].val)
						if temp != assign[j].val:
							assign[j].val = temp 
							if i not in remove_list:
								remove_list.append(i) 
			temp = {}
			for i in range(len(assign)):
			   if i not in remove_list:
				   temp[assign[i].var] = assign[i].val
			assign = temp 
		else:
		   
			list = line.split(' ') 
			transitions[int(list[0])][int(list[1])] = (list[2], list[3])
	
	if symbolic is True:
		print "[symbolic]"
	
	for i in range(num_states-1):
		matrix[i][i] = []
		for key, (rate, tr_type) in transitions[i].items():
			if symbolic is False:
				#matrix[i][key] = "-(%s)" % assign[rate]
				matrix[i][key] = "(%s)" % assign[rate]
				tr_type_matrix[i][key] = tr_type
				matrix[i][i].append("(" + assign[rate] + ")")
			else:
				#matrix[i][key] = "-(%s)" % rate
				matrix[i][key] = "(%s)" % rate
				matrix[i][i].append(rate)
		if len(matrix[i][i]) > 0:
			sum = "(%s" % matrix[i][i][0]
		for j in range(1, len(matrix[i][i])):
			sum = sum + "+%s" %  matrix[i][i][j]
		if len(matrix[i][i]) > 0:
			sum = sum + ")" 
		matrix[i][i] = "-(%s)" % sum
	
	matrix_out = []
	for i in range(len(matrix)):
		matrix_out.append([])
		for entry in matrix[i]:
			if entry != 0:
				matrix_out[i].append(mpf(eval(entry)))
			else:
				matrix_out[i].append(0)

	mc = MarkovChain(num_states, matrix_out, tr_type_matrix, type)

	return mc


def main():
	mc = build_markov_chain('../models/2DFT.disk.ls.critical.model', MarkovChain.REGULAR)

	print mc.prob_of_failure(t=range(0, 1*(87600+1), 1*(87600)))
	
	for i in range(1, 11):
		print "%e" % mc.estimate_prob_of_failure(t=i*8760)
	
	print mc.get_mttdl()

if __name__ == "__main__":
	main()
