from itertools import product
import sys
import copy
import numpy as np
import numpy.random
import numpy.linalg

# get indicies of nn_pairs_list, where the first base matches first_base
def get_first_base_indices(nn_pairs_list, first_base):
    # use list comprehension
    return [i for (i, j) in enumerate(nn_pairs_list) if j[0] == first_base]


# filter indicies based on picked_dict
def filter_list(picked_dict, indicies, max_num = 0):
    filt_indicies = []
    for ii in indicies:
        if picked_dict[ii] <= max_num:
            filt_indicies.append(ii)
    return filt_indicies


def sample_seq(nn_pairs_list, picked_list):
    seq = []
    spacer_len = len(picked_list)
    indicies = range(0, 16)

    for pos in range(0, spacer_len):
        picked_dict = picked_list[pos]
        if pos > 0:
            first_base = nn_pairs_list[seq[pos-1]][1]   #  = last base of last nn
            indicies = get_first_base_indices(nn_pairs_list, first_base)
        ok = False
        filt_indicies = []
        max_num = 0
        while not ok:
            filt_indicies = filter_list(picked_dict, indicies, max_num)
            if len(filt_indicies) > 0:
                ok = True
            else:
                max_num += 1
                print("WARNING: increasing max_num, at pos=" + str(pos))

        nnpair = numpy.random.choice(filt_indicies)
        seq.append(nnpair)

    return seq


def update_picked_list(picked_list, seq):
    for pos in range(0, len(seq)):
        picked_list[pos][seq[pos]] += 1
    return


def idx_to_str(nn_pairs_list, seq):
    seqstr = ""
    for pos in range(0, len(seq)):
        this_pair = nn_pairs_list[seq[pos]]
        if pos == 0:
            seqstr += this_pair[0] +  this_pair[1]
        else:
            seqstr += this_pair[1]
    return seqstr


def add_matrix_row(row, seqs_mat, seq):
    for pos in range(0, len(seq)):
        mat_col = (pos * 16) + seq[pos]
        seqs_mat[row, mat_col] += 1
    return


spacer_len = 20

picked_dict = {
    0: 0,
    1: 0,
    2: 0,
    3: 0,
    4: 0,
    5: 0,
    6: 0,
    7: 0,
    8: 0,
    9: 0,
    10: 0,
    11: 0,
    12: 0,
    13: 0,
    14: 0,
    15: 0
}

picked_list = []
for ii in range(0, spacer_len):
    picked_list.append(copy.deepcopy(picked_dict))

blank_picked_list = copy.deepcopy(picked_list)

# picked_list[0][0] += 1 # test

print(picked_list)

bases = ['A','T','G','C']

bases_list = [bases, bases]

print(bases_list)

nn_pairs = product(*bases_list)
nn_pairs_list = list(nn_pairs)

for nn in nn_pairs_list:
    print(nn)

print(nn_pairs_list)

print(get_first_base_indices(nn_pairs_list, 'C'))

gen_rounds = 20

seqs_mat = numpy.zeros((gen_rounds * 16, 320))

all_seqs = []

row = 0

for rdn in range(0, gen_rounds):
    for ii in range(0, 16):
        seq = sample_seq(nn_pairs_list, picked_list)
        update_picked_list(picked_list, seq)
        print(str(seq))
        print(idx_to_str(nn_pairs_list, seq))
        all_seqs.append(seq)
        add_matrix_row(row, seqs_mat, seq)
        row += 1
    picked_list = copy.deepcopy(blank_picked_list)

print(seqs_mat.shape)

mat_rank = numpy.linalg.matrix_rank(seqs_mat)

print(str(mat_rank))
