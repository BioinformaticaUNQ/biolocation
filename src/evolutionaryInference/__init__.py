import os
import subprocess
from pathlib import Path

from Bio import SeqIO, Phylo
from Bio.Align.Applications import MuscleCommandline
from Bio.Cluster import treecluster
from Bio.Phylo import PhyloXML
from Bio.Phylo.Applications import PhymlCommandline
from Bio.Phylo.PhyloXML import Phylogeny
from io import StringIO
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
from ete3 import Tree, phyloxml, PhyloTree, PhyloNode, TreeStyle, NodeStyle
from matplotlib import get_label
import numpy as np
from scipy.cluster.hierarchy import average, dendrogram

from src.evolutionaryInference import clustalo


def fasta_to_tree(filename):
    '''
    file = open(filename, "r")
    records = [seqrec for seqrec in SeqIO.parse(file, "fasta")]
    project_dir = Path(__file__).parent.parent.parent
    muscle_exe = os.path.join(project_dir, "resources/muscle3.8.31_i86/muscle3.8.31_i86win32.exe")
    # muscle_exe = os.path.join(project_dir, "resources/muscle3.8.31_i86/muscle3.8.31_i86linux32")

    muscle_cline = MuscleCommandline(muscle_exe, input=filename, clwstrict=True)
    stdout, stderr = muscle_cline()
    # cmdline = MuscleCommandline(muscle_exe, input=filename, out="egfr-family.aln", clw=True)
    # cmdline()

    # align = AlignIO.read(StringIO(stdout), "clustal")
    # align = AlignIO.read(StringIO(stdout), "fasta")
    AlignIO.convert(StringIO(stdout), "clustal",  "egfr-family.phy", "phylip-relaxed")
    # alignments = AlignIO.parse("egfr-family.aln", format="clustal", alphabet="phylip-relaxed")
    # calculator = DistanceCalculator('blastn')
    # dm = calculator.get_distance(align)
    # constructor = DistanceTreeConstructor()
    # tree = constructor.nj(dm)
    # fig = plt.figure(figsize=(10, 20), dpi=100)
    # axes = fig.add_subplot(1, 1, 1)
    # Phylo.draw(tree, axes=axes, do_show=False)

    iqtree_exe = os.path.join(project_dir, 'resources/iqtree-1.6.12-Windows/bin/iqtree.exe')
    # p = subprocess.run([iqtree_exe, "-s", "egfr-family.phy", "-m", "MF", "-nt", "AUTO", "-alrt", "1000", "-bb", "1000", "-bo", "100", "-pre", "iq_", "-redo"], stdout=subprocess.PIPE)
    # subprocess.run([iqtree_exe, "-s", "egfr-family.phy", "-m", "TIM3+R5", "-alrt", "1000", "-bb", "1000",
    #                "-bo", "100", "-pre", "iq_", "-wsr", "-wbt", "-nt", "AUTO"])
    subprocess.run([iqtree_exe, "-s", "egfr-family.phy", "-m", "TIM3+R5", "-alrt", "1000", "-bo", "100","-wbtl", "-nt", "AUTO", "-redo"])
    os.remove('egfr-family.phy')
    seqTree = open("egfr-family.phy.bionj", "r")
    handle = StringIO(str(seqTree.readlines().__getitem__(0)))
    tree = Phylo.read(handle, "newick")
    # phy_exe = os.path.join(project_dir, "resources/PhyML-3.1/PhyML-3.1_win32.exe")
    # cmdline = PhymlCommandline(phy_exe, input='egfr-family.phy', datatype='nt', alpha='e', bootstrap=10, search='BEST')
    # out_log, err_log = cmdline()
    # tree = Phylo.read("egfr-family.phy_phyml_tree.txt", "newick")

    tree.rooted = True
    egfr_phy = tree.as_phyloxml()
    egfr_phy = Phylogeny.from_tree(egfr_phy)
    egfr_phy.root.color = (128, 128, 128)
    egfr_phy.root.color = "#808080"  # This is one alternative
    egfr_phy.root.color = "gray"
    egfr_phy.clade[0].color = "blue"
    egfr_phy.clade[1].color = "red"
    # egfr_phy.clade[2, 0, 0].color = "green"
    mrca = egfr_phy.common_ancestor({"name": "JX912364.1"})
    mrca.color = "salmon"
    # Phylo.draw(egfr_phy)

    file = open('egfr-family.phy.mldist', 'r')
    size = int(file.readline())
    M = np.zeros((size, size))
    labels = []
    for i in range(size):
        row = file.readline().split(' ')
        labels.append(row[0])
        for j in range(i):
            M[i][j] = row[j+1]
    tree1 = treecluster(None, method='a', dist='e', distancematrix=M)
    tree = average(M)
    dendrogram(tree, orientation='left', labels=np.array(labels))
    plt.show()
    return tree
    '''

    # out_file = os.path.basename(filename) + '-aligned.fasta'
    # clustalo.runClustalO("grupo6@bioinformatica.com", filename, outfilename=out_file, fmt='clustal')
    # AlignIO.convert(out_file, "clustal", "egfr-family.phy", "phylip-relaxed")
    seqTree = open("egfr-family.phy.bionj", "r")
    t = Tree(str(seqTree.readlines().__getitem__(0)))
    Tree.convert_to_ultrametric(t)
    # nstyle = NodeStyle()
    # nstyle["fgcolor"] = "red"
    # t.get_leaves_by_name("JX912364.1")[0].set_style(nstyle)
    t.show()

    # Make a lookup table for sequences
    # lookup = dict((rec.id, str(rec.seq)) for rec in records)
    # for clade in egfr_phy.get_terminals():
    #    key = clade.name
        # accession = PhyloXML.Accession(key, 'NCBI')
        # mol_seq = PhyloXML.MolSeq(lookup[key], is_aligned=True)
        # sequence = PhyloXML.Sequence(type='aa', accession=accession, mol_seq=mol_seq)
        # clade.sequences.append(sequence)
    # Save the annotated phyloXML file
    # Phylo.write(egfr_phy, 'egfr-family-annotated.xml', 'phyloxml')
    # egfr = Phylo.read('egfr-family-annotated.xml', 'phyloxml')
    # print(egfr)
    # Phylo.draw_ascii(egfr)
