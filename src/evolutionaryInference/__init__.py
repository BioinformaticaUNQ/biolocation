import Bio
from scipy.cluster.hierarchy import dendrogram

from src.evolutionaryInference import alignementseq, alignementseq_multiple, tree_generator


def fasta_to_tree(filename):
    file = open(filename, "r")
    records = [seqrec for seqrec in Bio.SeqIO.parse(file, "fasta")]
    tree = tree_generator.tree_build(records, alignementseq.align2steps)
    aseqs = alignementseq_multiple.multiple_align_from_linkage_matrix(records, tree)
    alignedSeqs = [s.seq for s in aseqs]
    file.close()
    # ESTA ALINEACION DIFIERE DE LA DE CLUSTAL. REVISAR
    for s in aseqs:
        print(s.seq + '\n')
    dendrogram(tree)
    return alignedSeqs