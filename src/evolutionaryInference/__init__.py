import os
import subprocess
from pathlib import Path
from Bio import AlignIO
from ete3 import Tree

from src.evolutionaryInference import clustalo

project_dir = Path(__file__).parent.parent.parent

def fasta_to_tree(filename, aligned, bootstrap):
    out_file = filename
    if not aligned:
        out_file = os.path.basename(filename) + '-aligned.fasta'
        clustalo.runClustalO("grupo6@bioinformatica.com", filename, outfilename=out_file, fmt='clustal')
    AlignIO.convert(out_file, "clustal", "egfr-family.phy", "phylip-relaxed")
    iqtree_exe = os.path.join(project_dir, 'resources/iqtree-1.6.12-Windows/bin/iqtree.exe')
    subprocess.run([iqtree_exe, '-s', "egfr-family.phy", '-bb', str(bootstrap)])
    # os.remove('egfr-family.phy')
    seqTree = open("egfr-family.phy.treefile", "r")
    tree_phy = Tree(str(seqTree.readlines().__getitem__(0)))
    seqTree.close()
    os.remove("egfr-family.phy.bionj")
    os.remove("egfr-family.phy.ckp.gz")
    os.remove("egfr-family.phy.mldist")
    return tree_phy


