from itertools import product
import sys
import copy
import numpy as np
import numpy.random
import numpy.linalg

from scipy.linalg import null_space

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
    nn_pairs_list_len = len(picked_list[0])
    indicies = range(0, nn_pairs_list_len)             #16

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


def add_matrix_row(row, seqs_mat, seq, nn_pairs_list_len):
    for pos in range(0, len(seq)):
        mat_col = (pos * nn_pairs_list_len) + seq[pos]
        seqs_mat[row, mat_col] += 1
    return


spacer_len = 21     # 22 bp long spacer

bases = ['A','T','G','C']    # ['A','B'] # ['A','T','G','C']    # ['A','B'] #

num_bases = len(bases)

bases_list = [bases, bases]

max_rank = (num_bases*num_bases*spacer_len) - ( (spacer_len - 1) * num_bases )

#print(bases_list)

nn_pairs = product(*bases_list)
nn_pairs_list = list(nn_pairs)

#for nn in nn_pairs_list:
#    print(nn)

#print(nn_pairs_list)

nn_pairs_list_len = len(nn_pairs_list)

picked_dict = {}

for ii in range(0, nn_pairs_list_len):
    picked_dict[ii] = 0

#print(picked_dict)

picked_list = []
for ii in range(0, spacer_len):
    picked_list.append(copy.deepcopy(picked_dict))

blank_picked_list = copy.deepcopy(picked_list)

#print(picked_list)

#print(get_first_base_indices(nn_pairs_list, 'C'))

gen_rounds = spacer_len + 4     # plus 4 additionals

seqs_mat = numpy.zeros((gen_rounds * nn_pairs_list_len, spacer_len*nn_pairs_list_len))

best_score = spacer_len*nn_pairs_list_len + 1

all_seqs = []
all_seqs_str = []

row = 0
rdn = 0

score = best_score
max_rank_reached = False

last_rank = 0

mat_info = []

#for rdn in range(0, gen_rounds):
while rdn < gen_rounds:
    all_uniq = True
    rdn_seqs = []
    rdn_seqs_str = []
    for ii in range(0, nn_pairs_list_len):
        seq = sample_seq(nn_pairs_list, picked_list)
        update_picked_list(picked_list, seq)
        seq_str = idx_to_str(nn_pairs_list, seq)
        #print(seq_str)
        if seq_str in all_seqs_str:
            all_uniq = False
        rdn_seqs.append(seq)
        rdn_seqs_str.append(seq_str)
    picked_list = copy.deepcopy(blank_picked_list)
    if all_uniq:
        # rdn accepted
        for rdn_ii in range(0, len(rdn_seqs)):
            seq = rdn_seqs[rdn_ii]
            seq_str = rdn_seqs_str[rdn_ii]
            all_seqs.append(seq)
            all_seqs_str.append(seq_str)
            #update_picked_list(picked_list, seq)
            add_matrix_row(row, seqs_mat, seq, nn_pairs_list_len)
            mat_rank = numpy.linalg.matrix_rank(seqs_mat)
            this_mat_info = str(rdn) + "\t" + str(rdn_ii)
            if mat_rank > last_rank:
                this_mat_info += "\t(+)"
            else:
                this_mat_info += "\t(-)"
            #print("Matrix_rank:\t" + str(row) + "\t" + str(mat_rank) + "\t" + str(max_rank))
            if mat_rank == max_rank and max_rank_reached == False:
                score = row
                max_rank_reached = True
            mat_info.append(this_mat_info)
            last_rank = mat_rank
            row += 1
        rdn += 1

np.set_printoptions(threshold=np.inf)

#print(all_seqs_str)

print(str(len(all_seqs_str)) + " were generated")
all_seqs_str_set = set(all_seqs_str)
print(str(len(all_seqs_str_set)) + " were unique")

#print(all_seqs_str_set)

# print(seqs_mat)

# print(seqs_mat.shape)

print("Max_Matrix_rank:\t" + str(max_rank))

mat_rank = numpy.linalg.matrix_rank(seqs_mat)
print("Final_Matrix_rank:\t" + str(mat_rank))

print("score:\t" + str(score))

print("outputting sequences")

f = open('designed_spacers.fa', 'w')

for ii in range(0, len(all_seqs_str)):
    f.write("> spacer" + str(ii) + "\t" + mat_info[ii] + "\n")
    f.write(all_seqs_str[ii] + "\n")

f.close()

ns = null_space(seqs_mat)

#print(ns)
print(ns.shape)
