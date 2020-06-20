#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 15:02:16 2019

@author: ariane.delrocq
"""

import numpy as np

import Bio.SubsMat.MatrixInfo
import Bio.SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import Bio.pairwise2

# seqs=[]
# handle = open("balibase/RV11.unaligned/BBS11001.fasta")
# for seq in Bio.SeqIO.parse(handle, "fasta"):
#    seqs.append(seq)
# handle.close()

mat = Bio.SubsMat.MatrixInfo.blosum62


def score(a, b, gap=-8):
    """Default scoring function for align
    """
    if not a or not b:
        return gap
    else:
        a = a.upper()
        b = b.upper()
        try:
            return mat[a, b]
        except KeyError:
            return mat[b, a]


def align(seq1, seq2, score=score):
    """Main implementation of align, with blosum scoring and no specific gap
    extension penalties.

    Args:
        seq1: A SeqRecord object
        seq2: Same.
    """

    x = seq1.seq
    y = seq2.seq
    n = len(x)
    m = len(y)
    mm = [[0 for _ in range(m + 1)] for _ in range(n + 1)]  # scores matrix
    path = [[0 for _ in range(m + 1)] for _ in range(n + 1)]  # last step for optimal score

    # first line and col
    for i in range(1, n + 1):
        mm[i][0] = mm[i - 1][0] + score('', x[i - 1])
        path[i][0] = (i - 1, 0)
    for j in range(1, m + 1):
        mm[0][j] = mm[0][j - 1] + score('', y[j - 1])
        path[0][j] = (0, j - 1)

    # fill table
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            a = x[i - 1]
            b = y[j - 1]
            s1 = mm[i - 1][j - 1] + score(a, b)
            s2 = mm[i - 1][j] + score('', a)
            s3 = mm[i][j - 1] + score('', b)
            if s1 >= s2:
                if s1 >= s3:
                    mm[i][j] = s1
                    path[i][j] = (i - 1, j - 1)
                else:
                    mm[i][j] = s3
                    path[i][j] = (i, j - 1)
            else:  # s2 > s1
                if s2 >= s3:
                    mm[i][j] = s2
                    path[i][j] = (i - 1, j)
                else:
                    mm[i][j] = s3
                    path[i][j] = (i, j - 1)

    bestScore = mm[n][m]
    alx = ""
    aly = ""
    i, j = n, m
    while i > 0 and j > 0:
        i2, j2 = path[i][j]
        if i == i2:
            alx += '-'
        else:
            alx += x[i - 1]
        if j == j2:
            aly += '-'
        else:
            aly += y[j - 1]
        i, j = i2, j2
    while i > 0:
        aly += '-'
        alx += x[i - 1]
        i -= 1
    while j > 0:
        alx += '-'
        aly += y[j - 1]
        j -= 1
    aseq1 = SeqRecord(seq1)
    aseq1.seq = Seq(alx[::-1])
    aseq2 = SeqRecord(seq2)
    aseq2.seq = Seq(aly[::-1])
    return bestScore, aseq1, aseq2


def dict_to_weight(subs_dict):
    """Converts blosum dict to indexable numpy array.

    The score function is the critical op here, we do not need a dict lookup in
    the innermost loop.

    Returns:
        26 x 26 numpy array, [0,0] being s('A', 'A')
    """

    weight = np.zeros((26, 26))
    for (r1, r2), w in subs_dict.items():
        r1 = ord(r1) - ord('A')
        r2 = ord(r2) - ord('A')
        weight[r1, r2] = weight[r2, r1] = w
    return weight


def vec_align(seq1, seq2, mat=dict_to_weight(Bio.SubsMat.MatrixInfo.blosum62), gap=-8):
    """Proof of concept vector-style implementation of align, written mainly in numpy
    arrays operations.

    M(i,j) is the score of the best alignement containing exactly i letters of
    x and j letters of y.
    To allow vector operations, we rewrite M as S(s, k) with s=i+j and
    k=i-max(0,s-m)
    Hence:
        M(i,j) = S(i+j, i-max(0,i+j-m))
        S(s, k) = M(k+max(0,s-m), min(s,m)-k)
    The relation become:
        M(i,j) = max( M(i,j-1)+gap, M(i-1,j)+gap, M(i-1,j-1)+s(i,j) )
        S(s,k) = max(
           S(s-1,k+1(s>m))+gap,
           S(s-1,k-1+1(s>m))+gap,
           S(s-2,k-1+1(s>m)+1(s>m+1)) + s(k+max(0,s-m),min(s,m)-k)
        )
    Note that S(s,k) does not depend anymore on S(k,.), allowing vector operation
    Last index of the diag s can be written (useful for simplifications):
        min(n,s)-max(0,s-m)
        = min(n, m, s, n+m-s)
        = min(m,s)-max(0,s-n)
    """

    n = len(seq1.seq)
    m = len(seq2.seq)
    x = np.array(list(ord(c) - ord('A') for c in seq1.seq))
    y = np.array(list(ord(c) - ord('A') for c in seq2.seq))
    mat = mat.flatten()
    S = [np.zeros(min(s, n, m, n + m - s) + 1, dtype=np.int32) for s in range(n + m + 1)]

    S[0][0] = 0
    S[1][0] = gap
    S[1][1] = gap

    x *= 26
    for s in range(2, n + m + 1):
        # First/last cell must be treated differently when first line/col of M
        if (s <= m):
            S[s][0] = S[s - 1][0] + gap
        if (s <= n):
            S[s][s - max(0, s - m)] = S[s - 1][s - 1 - max(0, s - 1 - m)] + gap

        # The offsets in the formula become begin/end of slices.
        gapped = S[s - 1] + gap
        S[s][int(s <= m):len(S[s]) - int(s <= n)] = np.maximum(np.maximum(
            gapped[1:],
            gapped[:-1]),
            S[s - 2][int(s > m + 1):len(S[s - 2]) - int(s > n + 1)]
            + mat[
                x[max(0, s - m - 1):min(n, s - 1)] + y[max(0, s - n - 1):min(m, s - 1)][::-1]
                ]
        )

    # Trace back
    score = S[n + m][0]
    i, j = n, m
    alx, aly = "", ""
    while i > 0 or j > 0:
        if j > 0:
            prev_score = S[i + j - 1][i - max(0, i + j - m - 1)]
            if prev_score + gap == score:
                score = prev_score
                j -= 1
                aly += seq2[j]
                alx += '-'
                continue
        if i > 0:
            prev_score = S[i + j - 1][i - max(0, i + j - m - 1) - 1]
            if prev_score + gap == score:
                score = prev_score
                i -= 1
                alx += seq1[i]
                aly += '-'
                continue
        score = S[i + j - 2][i - max(0, i + j - m - 2) - 1]
        i -= 1
        j -= 1
        alx += seq1[i]
        aly += seq2[j]

    aseq1 = SeqRecord(seq1)
    aseq1.seq = Seq(alx[::-1])
    aseq2 = SeqRecord(seq2)
    aseq2.seq = Seq(aly[::-1])
    return S[n + m][0], aseq1, aseq2


def align2steps(seq1, seq2, d=8, e=4):
    """Main implementation of align, handling also gap extensions.
    """
    x = seq1.seq
    y = seq2.seq
    n = len(x)
    m = len(y)
    mm = [[0 for _ in range(m + 1)] for _ in range(n + 1)]  # scores matrix
    mx = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    my = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    pathm = [[0 for _ in range(m + 1)] for _ in range(n + 1)]  # last step for optimal score
    pathx = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    pathy = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

    lower_bound = - (n + m) * max(e, d)

    # first line and col
    for i in range(1, n + 1):
        mx[i][0] = -d - (i - 1) * e
        pathx[i][0] = pathx
        # Best alignment ending with no gap but with 0 letters of one seq does
        # not exist... but is still used in the formulas below
        mm[i][0] = lower_bound
        # Same for best alignment with 0 letter of Y ending on a gap on the X
        # side with a letter of Y
        my[i][0] = lower_bound

    for j in range(1, m + 1):
        my[0][j] = -d - (j - 1) * e
        pathy[0][j] = pathy
        mm[0][j] = lower_bound
        mx[0][j] = lower_bound

    # fill table
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            a = x[i - 1]
            b = y[j - 1]

            # find max for M
            s1 = mm[i - 1][j - 1] + score(a, b)
            s2 = mx[i - 1][j - 1] + score(a, b)
            s3 = my[i - 1][j - 1] + score(a, b)

            if s1 >= s2:
                if s1 >= s3:
                    mm[i][j] = s1
                    pathm[i][j] = pathm
                else:
                    mm[i][j] = s3
                    pathm[i][j] = pathy
            else:  # s2 > s1
                if s2 >= s3:
                    mm[i][j] = s2
                    pathm[i][j] = pathx
                else:
                    mm[i][j] = s3
                    pathm[i][j] = pathy

            # find max for I_x
            s4 = mm[i - 1][j] - d
            s5 = mx[i - 1][j] - e
            if s4 >= s5:
                mx[i][j] = s4
                pathx[i][j] = pathm
            else:
                mx[i][j] = s5
                pathx[i][j] = pathx

            # find max for I_y
            s6 = mm[i][j - 1] - d
            s7 = my[i][j - 1] - e
            if s6 >= s7:
                my[i][j] = s6
                pathy[i][j] = pathm
            else:
                my[i][j] = s7
                pathy[i][j] = pathy

    bestScore, bestPath = max(
        [(mm, pathm), (mx, pathx), (my, pathy)],
        key=(lambda t: t[0][n][m]))
    bestScore = bestScore[n][m]
    alx = ""
    aly = ""
    i, j = n, m
    path = bestPath
    while i > 0 and j > 0:
        prev_path = path[i][j]
        if path == pathm:
            i -= 1
            j -= 1
            alx += x[i]
            aly += y[j]
        elif path == pathy:
            j -= 1
            alx += '-'
            aly += y[j]
        else:
            i -= 1
            alx += x[i]
            aly += '-'
        path = prev_path

    while i > 0:
        aly += '-'
        alx += x[i - 1]
        i -= 1
    while j > 0:
        alx += '-'
        aly += y[j - 1]
        j -= 1
    alx = alx[::-1]
    aly = aly[::-1]
    aseq1 = SeqRecord(seq1)
    aseq1.seq = Seq(alx)
    aseq2 = SeqRecord(seq2)
    aseq2.seq = Seq(aly)
    return bestScore, aseq1, aseq2


def test():
    """Obsolete test function, to be moved to unit tests in a separate file.
    """
    #    print(align("chat", "cat"))
    #    print(align("chat", "cgat"))
    #    print(align("chat", "at"))
    s, x, y = align2steps(seqs[0], seqs[1])
    print(s)
    print(x.seq)
    print(y.seq)
    with open("balibase/us/1.fasta", "w") as fd:
        Bio.SeqIO.write((x, y), fd, "fasta")
    print("Ref:")
    ref = Bio.pairwise2.align.globalds(
        seqs[0].seq, seqs[1].seq,
        Bio.SubsMat.MatrixInfo.blosum62, -8, -4)

    for x, y, s, *_ in ref:
        print(s)
        print(x)
        print(y)
