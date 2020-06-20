# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 14:59:19 2019

@author: Amaury Leroy
"""

import numpy as np
from scipy.cluster.hierarchy import average

'''seqs = [SeqRecord('TAACCCCAAAAGAACCAGA'), SeqRecord('TTCTGTAGTAGAGATGGAATTAAGAAAAAACCATCAACTATAACCCCAAAAGAACCAGA'),
        SeqRecord('TTCTGTAGTAGAGATGGAATTAAGAAAAAACCATCAACTATAACCCCAAAAGAACCAGA'),
        SeqRecord('GTAGTAGAGATGGAATTAAGAAAAAACCATCAACTATAACCCCAAGAGAACCAGA'),
        SeqRecord('GAGCCGGATGAGAAGAAACTCTCATGTCCGGTTCTGTAGTAGAGATGGAATTAAGAAAAAACCATCAACTATAACCCCAAGAGAACCAGA'),
        SeqRecord(
            'TTTTCATTCGCGAGGAGCCGGATGAGAAGAAACTCTCATGTCCGGTTCTGTAGTAGAGATGGAATTAAGAAAAAACCATCAACTATAACCCCAAGAGAACCAGA')]
'''
# Arguments : une liste de séquences SeqRecord, une fonction d'alignement
score_gap = -8


def tree_build(seqs, align):
    len_seqs = [len(seq) for seq in seqs]
    M = np.zeros((len(seqs), len(seqs)))
    for i in range(len(seqs)):
        for j in range(i):
            M[i][j], a, b = align(seqs[i], seqs[j])
            M[i][j] = 1 - abs(M[i][j] / (max(len_seqs) * (11 - score_gap)))
    print(type(M), M)
    #    tree1 = treecluster(None, method = 'a', dist = 'e', distancematrix=M)
    tree = average(M)
    return tree


# evolutionaryInference = tree_build(seqs, alignementseq.align2steps)
# print(evolutionaryInference)
# dn = dendrogram(evolutionaryInference)
