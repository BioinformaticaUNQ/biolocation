from src.evolutionaryInference.alignementseq import score
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def multiple_align_from_tree(tree):
    """Multiple alignmnent from a SeqRecord evolutionaryInference

    Arguments:
        tree: nested structure of python tuples, with SeqRecord objects as leafs
    """

    def _join(node):
        """Recursive helper to browse evolutionaryInference

        Arguments:
            node: (node, node) | SeqRecord
        """
        if not isinstance(node, tuple):
            return [node]
        left, right = node
        aleft = _join(left)
        aright = _join(right)
        return join_alignments(aleft, aright)

    return _join(tree)


def multiple_align_from_linkage_matrix(sequences, tree):
    """Multiple alignmnent from a SeqRecord evolutionaryInference

    Arguments:
        sequences: indexable container of SeqRecord objects to align
        tree: scipy.hierarchy linkage matrix
    """

    n = len(sequences)

    def _join(inode):
        """Recursive helper to browse evolutionaryInference

        Arguments:
            inode: index of node row in linkage matrix
        """
        inode = int(inode)
        if inode < n:
            return [sequences[inode]]
        ileft, iright, dist, size = tree[inode - n]
        aleft = _join(ileft)
        aright = _join(iright)
        return join_alignments(aleft, aright)

    return _join(n + len(tree) - 1)


def _score(a, b, base_score=score):
    """Generalized one-letter score, with '-' as new letter.
    """

    if a == '-':
        if b == '-':
            return 0
        else:
            return score('', b)
    else:
        if b == '-':
            return score('', a)
        else:
            return score(a, b)


def _mscore(seqs1, seqs2, i, j, _score=_score):
    """Multiple score
    """
    return sum(_score(x[i], y[j]) for x in seqs1 for y in seqs2)


def _mscore_gap(seqs, i, gap):
    """Multiple score for a global single gap against aligned letters

    One gap penalty for each non-gap in seqs at pos i.
    """
    return gap * sum(x[i] != '-' for x in seqs)


def join_alignments(aleft, aright, _mscore=_mscore, gap=-8):
    """
    Takes two multiple alignments as lists of SeqRecords, and align them.

    Returns:
        Multiple align as list of aligned SeqRecords
    """

    n = len(aleft[0])
    m = len(aright[0])
    scores = [([0] * (m + 1)) for _ in range(n + 1)]
    paths = [([0] * (m + 1)) for _ in range(n + 1)]

    for i in range(1, n + 1):
        scores[i][0] = scores[i - 1][0] + _mscore_gap(aleft, i - 1, gap)
        paths[i][0] = (i - 1, 0)
    for j in range(1, m + 1):
        scores[0][j] = scores[0][j - 1] + _mscore_gap(aright, j - 1, gap)
        paths[0][j] = (0, j - 1)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            scores[i][j], paths[i][j] = max(
                (
                    scores[i - 1][j - 1] + _mscore(aleft, aright, i - 1, j - 1),
                    (i - 1, j - 1)
                ), (
                    scores[i - 1][j] + _mscore_gap(aleft, i - 1, gap),
                    (i - 1, j)
                ), (
                    scores[i][j - 1] + _mscore_gap(aright, j - 1, gap),
                    (i, j - 1)
                ),
                key=(lambda v: v[0])
            )

    best_score = scores[n][m]
    i, j = n, m
    aseqs = [""] * (len(aleft) + len(aright))

    while i > 0 and j > 0:
        nexti, nextj = paths[i][j]

        if nexti == i:
            for iseq in range(len(aleft)):
                aseqs[iseq] += '-'
        else:
            for iseq in range(len(aleft)):
                aseqs[iseq] += aleft[iseq][nexti]
        if nextj == j:
            for jseq in range(len(aright)):
                aseqs[len(aleft) + jseq] += '-'
        else:
            for jseq in range(len(aright)):
                aseqs[len(aleft) + jseq] += aright[jseq][nextj]
        i, j = nexti, nextj

    ret = [SeqRecord(s) for s in aleft] + [SeqRecord(s) for s in aright]
    for rec, seq in zip(ret, aseqs):
        rec.seq = Seq(seq[::-1])
    return ret

