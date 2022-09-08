# model_grna_seq_des

The purpose of this script is to produce spacer sequences for a training set for nearest neighbour models of gRNA binding to a cognate binding sequence.

The script produces a set of 16 sequences by sampling nearest neighbour pairs in a stepwise manner for left to right 16 times until all nearest neighbour pairs at each position have been sampled. The script avoids sampling the same NN paur twice by recording previously sampled pairs.

The script repeats the production of sets of 16 sequences gen_round (=21+4=25) times producing 400 unique sequences. 

The generated sequences are also stored in a matrix format, enabling the calculation of matrix rank. 

One can express the testing of set of sequences in matrix form, A: 

Ax = b

where A is the set of sequences in matrix form where each row is a sequence and the columns are the nearest-neighbour pairs at each position, x is the vector of unknown parameters, b is the vector of binding free-energies.

The aim is to design matrix A such that it has full rank so that each sequences gives new information for training the model.



