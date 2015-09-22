import numpy
import math
import random
import sys
from bm_ops import *
from util import *

class big_bm:
    global word_size
    word_size=32
    def __init__(self, size):
        self.size = size
        self.size_in_words = int(math.ceil(float(size)/(word_size)))
        self.num_full_words = self.size_in_words
        self.shortened_leader = False
        self.array = numpy.ndarray(self.size_in_words)
        self.array.dtype = long
        self.array.fill(0)
        if size%word_size != 0:
            self.num_full_words -= 1
            self.shortened_leader = True
        self.leading_word_size = self.size-(self.num_full_words*word_size)

    def copy(self):
        obj = big_bm(self.size)
        obj.size = self.size 
        obj.size_in_words = self.size_in_words 
        obj.num_full_words = self.num_full_words 
        obj.shortened_leader = self.shortened_leader 
        obj.array = numpy.ndarray(obj.size_in_words)
        obj.array.dtype = long
        for i in range(self.array.size):
            obj.array[i] = self.array[i] 
        obj.leading_word_size = self.leading_word_size 

        return obj

    def __str__(self):
        str=""
        mask = 1

        if self.shortened_leader is True:
            for j in range(self.leading_word_size):
                str+="%d" % ((int(self.array[self.size_in_words-1]) & (mask << (self.leading_word_size-1-j))) >> (self.leading_word_size-1-j))

        for i in range(self.num_full_words):
            for j in range(word_size):
                str+="%d" % ((int(self.array[self.num_full_words-1-i]) & (mask << (word_size-1-j))) >> (word_size-1-j))

        return str

    def set_all_ones(self):
        self.array.fill((1 << word_size) - 1)
            
    def set_all_zeroes(self):
        self.array.fill(0)

    def band(self, result, x):
        if x.size_in_words != self.size_in_words:
            return None
 
        for i in range(self.size_in_words):
            result.array[i] = self.array[i] & x.array[i]

        return 0

    def xor(self, result, x):
        if x.size_in_words != self.size_in_words:
            return None
 
        for i in range(self.size_in_words):
            result.array[i] = self.array[i] ^ x.array[i]

        return 0

    def set_bits(self, bit_index):
        for idx in bit_index:
            self.array[idx/word_size] |= (1 << (idx%word_size)) 
    
    def unset_bits(self, bit_index):
        for idx in bit_index:
            self.array[idx/word_size] ^= (1 << (idx%word_size)) 
    
    def is_set(self, bit_index):
        word_index = bit_index / word_size
        if word_index < self.num_full_words-1:
            index = bit_index % word_size
        elif self.shortened_leader is True:
            index = bit_index % word_size
        elif self.shortened_leader is False:
            index = bit_index % word_size
        return ((self.array[word_index] & (1 << index)) >> index)

    def randomize(self):
        for i in range(self.num_full_words):
            self.array[i] = random.getrandbits(word_size)
        for i in range(self.size_in_words-self.num_full_words):
            self.array[i] = random.getrandbits(self.leading_word_size)

    def get_weight(self):
        weight = 0
        mask = 1 
        for i in range(self.size_in_words):
            for j in range(word_size):
                weight+=((self.array[i] & (mask << j)) >> j)
        return weight


    def get_leading_bit(self):
        word_num=self.size_in_words-1
        bit_num = self.size_in_words*word_size-1
        while 1:
            if self.array[word_num] != 0:
                break
            bit_num-=word_size
            word_num-=1
            if word_num < 0:
                return -1

        mask = 1 << (word_size-1)

        while 1:  
            if mask & self.array[word_num] != 0:
                break
            bit_num -= 1
            mask >>= 1

        return bit_num
    

def test():
    bm = big_bm(50)
    
    bm.set_bits([1,47])
    
    print bm.is_set(43)
    print bm.is_set(47)
    
    print bm.get_leading_bit()
    
    print bm

if __name__ == "__main__":
    test()
