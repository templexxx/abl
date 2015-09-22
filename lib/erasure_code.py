from bit_matrix import *
from numpy import random
from util import kset_set

##
# This module encapsulates information and functions for 
# an arbitrary erasure code.  The code description information
# is read from a file (given as an argument to the constructor).
#
# Kevin M. Greenan (kmgreen@cs.ucsc.edu)
#

##
# Static location of all code description files
#
code_dir="./codes/"

##
# This class encapsulates the number of data/parity symbols,
# the code type, the tanner graph (for XOR-based codes) and 
# the generator matrix (for XOR-based codes).
#
class ErasureCode:
	TYPE_MDS = "mds"
	TYPE_FLAT_XOR = "flat xor"
	TYPE_ARRAY_XOR = "array xor"
	CHECK_RANK = "rank_check"
	CHECK_MEL = "mel_check"
	CHECK_FTV = "ftv_check"
	CHECK_DSCFT = "disk sector conditional fault tolerance check"

	def __init__(self, code_file, fail_check_type=None):
		file = open(code_dir+code_file, "r")
		
		if fail_check_type is None:
			self.fail_check_type = self.CHECK_RANK
		else:
			self.fail_check_type = fail_check_type
		
		self.hd = None

		line = file.readline().strip()

		while line != '':

			if line == '[type]':
				line = file.readline().strip()
				self.type = line
			elif line == '[k]':
				line = file.readline().strip()
				self.k = int(line)
			elif line == '[m]':
				line = file.readline().strip()
				self.m = int(line)
			elif line == '[hd]':
				line = file.readline().strip()
				self.hd = int(line)
			elif line == '[min disk failures]':
				line = file.readline().strip()
				self.min_disk_failures = int(line)
			elif line == '[tanner graph]':
				line = file.readline().strip()
				self.tanner_graph = eval(line)
			elif line == '[raw layout]':
				line = file.readline().strip()
				self.layout = eval(line)
			elif line == '[minimal fault sets]':
				self.mel_bm = []
				while 1:
					line = file.readline().strip()
					if line ==  '[END]':
						break
					self.mel_bm.append(list_to_bm(eval(line)))
			elif line == '[Disk sector conditional fault tolerance]':
				line = file.readline().strip()
				self.dsft = eval(line)
			elif line == '[fault tolerance vector]':
				line = file.readline().strip()
				self.ftv = eval(line)
				
			line = file.readline().strip()

		if self.type == self.TYPE_FLAT_XOR or self.type == self.TYPE_ARRAY_XOR:
			self.generator_matrix = build_generator(self.k, self.m, self.tanner_graph)

		# We assume that the HD must be at least 2.  If it is not specified, then set to 2
		if self.hd is None:
			self.hd = 2

	##
	# Determine if this set of failed symbols results in a failure.
	#
	# For MDS codes, this is a simple check to see if the number of 
	# failed symbols is greater than the number of parity symbols.
	#
	# For XOR-based codes, we must resort to either a lookup to the 
	# minimal erasures list or perform a rank test on a modified 
	# generator matrix.
	#
	# @param failed_symbols: list of failed symbol ids
	# @return True if no failure or False if there is a failure 
	#
	def is_failure(self, failed_disks, failed_sectors=[]):
		symbol_errors = []
		unique_sectors = {}
		
		if len(failed_sectors) == 0:
			unique_sectors[0] = []
		
		# Determine the symbol errors
		if self.type == self.TYPE_MDS or self.type == self.TYPE_FLAT_XOR:
			
			symbol_errors = failed_disks[:]
			
			for comp in range(len(failed_sectors)):
				if comp in failed_disks: continue
				for sector in failed_sectors[comp]:
					if not unique_sectors.has_key(sector):
						unique_sectors[sector] = []
					
					unique_sectors[sector].append(comp)
					
		elif self.type == self.TYPE_ARRAY_XOR:
			for i in range(len(failed_disks)):
				for sym_idx in self.layout[failed_disks[i]]:
					symbol_errors.append(sym_idx)
			
			for comp in range(len(failed_sectors)):
				if comp in failed_disks: continue
				for sector in failed_sectors[comp]:
					stripe_num = sector/len(self.layout[comp])
					if not unique_sectors.has_key(stripe_num):
						unique_sectors[stripe_num] = []

					unique_sectors[stripe_num].append(self.layout[i][sector % len(self.layout[i])])
		
		# Determine if there is a failure	
		if self.type == self.TYPE_MDS:
			for sectors in unique_sectors.values():
				if self.m < len(symbol_errors+sectors):
					return True
			return False
		
		if self.fail_check_type == self.CHECK_RANK:
			
			for sectors in unique_sectors.values():
			
				temp_matrix = self.generator_matrix.copy()
			
				temp_matrix.zero_cols(symbol_errors + sectors)
			
				if self.k > get_rank(temp_matrix):
					return True
				else:
					return False
			
		elif self.fail_check_type == self.CHECK_MEL:
			for sectors in unique_sectors.values():
				
				symbol_errors_bm = list_to_bm(symbol_errors+sectors)
			
				for me_pattern in self.mel_bm:
					if bm_intersection(symbol_errors_bm, me_pattern) == me_pattern:
						return True
					
			return False
		
		elif self.fail_check_type == self.CHECK_FTV:
			
			for sectors in unique_sectors.values():
			
				draw = random.uniform()

				if draw < self.ftv[len(symbol_errors+sectors)-1]:
					return True
			
			return False
		 
		elif self.fail_check_type == self.CHECK_DSCFT:
			if len(failed_disks) >= len(self.dsft):
				return True
			
			draw = random.uniform()
			
			if draw < self.dsft[len(failed_disks)][0]:
				return True
			
			for sectors in unique_sectors.values():
				if len(sectors) == 0: 
					continue
				draw = random.uniform()
				if draw < self.dsft[len(failed_disks)][len(sectors)]:
					return True
			return False
		
					
def test():
	ec = ErasureCode("5_3_flat", ErasureCode.CHECK_FTV)

	print "Tanner graph: ", ec.tanner_graph
	print "Generator Matrix:"
	print ec.generator_matrix

	failure_list = kset_set([i for i in range(ec.k+ec.m)], 3, 0, [], [])

	for failure in failure_list:
		print ec.is_failure(failure)



if __name__ == "__main__":
	test()
