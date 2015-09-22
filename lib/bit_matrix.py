import math
import random
import sys
from bm_ops import *
from big_bm import *
from util import *

class bit_matrix:
    def __init__(self, num_rows, num_cols):
        self.num_rows = num_rows
        self.num_cols = num_cols
        self.matrix = []
        for i in range(num_rows):
            self.matrix.append(big_bm(num_cols)) 

    def copy(self):
        obj = bit_matrix(self.num_rows, self.num_cols)
        for i in range(obj.num_rows):
            obj.matrix[i] = self.matrix[i].copy()
        return obj

    def __str__(self):
        str=""
        for row in self.matrix:
            str+="%s, (%s)\n" % (row, row.get_weight())
        return str

    def swap_rows(self, row1, row2):
        tmp = self.matrix[row1]
        self.matrix[row1] = self.matrix[row2]
        self.matrix[row2] = tmp

    def set_row(self, row_idx, bit_idxs):
        self.matrix[row_idx].set_bits(bit_idxs)
    
    def get_row(self, row_idx):
        return self.matrix[row_idx]

    def is_set(self, row_idx, bit_idx):
        return self.matrix[row_idx].is_set(bit_idx)

    def zero_cols(self, col_idxs):
        mask = big_bm(self.num_cols)
        mask.set_all_ones()
        mask.unset_bits(col_idxs)

        for i in range(self.num_rows):
            self.matrix[i].band(self.matrix[i], mask)

    def zero_rows(self, row_idxs):
        for row_idx in row_idxs:
            self.matrix[row_idx].set_all_zeroes()

    def randomize(self):
        for i in range(self.num_rows):
            self.matrix[i].randomize()

    def band_store(self, row1, row2):
        self.matrix[row1].band(self.matrix[row1], self.matrix[row2])

    def xor_store(self, row1, row2):
        self.matrix[row1].xor(self.matrix[row1], self.matrix[row2])

def get_rank(matrix):
    rank = 0
    curr_col =  matrix.num_cols-1
    i=0

    while i < matrix.num_rows and curr_col >= 0:
        swap = False
        for j in range(i, matrix.num_rows):
            if matrix.matrix[j].is_set(curr_col) == 1:
                swap = True
                # Swap the rows                 
                matrix.swap_rows(i, j)

                rank+=1

                break 

        if swap is True:
            for j in range(i, matrix.num_rows):
                if j != i and matrix.matrix[j].is_set(curr_col) == 1:
                    matrix.xor_store(j, i)
                    
        if swap is False:
            curr_col-=1
        else:
            curr_col-=1
            i += 1
    return rank

# Make a systematic generator matrix  G=[I_num_data | P]
def build_generator(num_data, num_parity, parity_eqns):
    generator_matrix = bit_matrix(num_data, num_parity+num_data)

    # Fill in row by row
    for i in range(num_data):
        row_els = [i]

        for j in range(num_parity):
            if i in parity_eqns[j]:
                row_els.append(j+num_data)

        generator_matrix.set_row(i, row_els)
    return generator_matrix

def unit_tests():

    # For an m x n matrix, zeroing out (n-m)+1 cols 
    # should result in rank < m
    ms = [24, 32, 48, 64, 72, 96, 128]
    n = 128

    for m in ms:
        matrix = bit_matrix(m, n)
        matrix.randomize()
        zero_cols = []
        while len(zero_cols) < (n-m)+1:
            col = random.randint(0, n-1)
            if col not in zero_cols:
                zero_cols.append(col)

        matrix.zero_cols(zero_cols)
    
        rank = get_rank(matrix)

        if rank <= m:
            print "Passed test"
        else:
            print "Failed test"

    # Construct matrix that should survive any k zeroed columns
    # with the exception of a few failure patterns.
    # Test to ensure rank == m for any k zeroed columns and 
    # the specific set of failure patterns.
    k=2
    n=8
    m=5
    matrix = bit_matrix(m, n)
    matrix.set_row(0, [0,5,6,7])
    matrix.set_row(1, [1,5])
    matrix.set_row(2, [2,5,6,7])
    matrix.set_row(3, [3,6,7])
    matrix.set_row(4, [4,7])

    will_fail_if = [list_to_bm([1,5]), list_to_bm([4,7]), list_to_bm([0,2])]
    
    all_els=[i for i in range(n)]
    fail_combs = []
    generate_combinations(all_els, k, fail_combs)

    for fail_patt in fail_combs:
        tmp_matrix = matrix.copy()
        tmp_matrix.zero_cols(fail_patt)
        rank = get_rank(tmp_matrix)
        print "Failure pattern : ", fail_patt
        print "Rank : ", rank
        print tmp_matrix
        if rank >= m:
            print "Passed test"
        elif list_to_bm(fail_patt) in will_fail_if:
            print "Passed test"
        else:
            print "Failed test"
            print "Failure pattern %s" % fail_patt

