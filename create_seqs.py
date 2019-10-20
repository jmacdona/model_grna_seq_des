from itertools import product
import sys


bases = ['A','T','G','C']

bases_list = [bases, bases]

print(bases_list)


nn_pairs = product(*bases_list)


for nn in nn_pairs:
    print(nn)

print(nn_pairs)