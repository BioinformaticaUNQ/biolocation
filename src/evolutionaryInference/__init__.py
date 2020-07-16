import os
import shutil
import subprocess
from pathlib import Path
from Bio import AlignIO
from ete3 import Tree

from src.evolutionaryInference import clustalo

project_dir = Path(__file__).parent.parent.parent

def fasta_to_tree(filename, aligned, bootstrap, quantitySequences):
    nameFile = os.path.basename(filename).split('.')[0]
    out_file = filename
    if not aligned:
        out_file = nameFile + '-aligned.fasta'
        clustalo.runClustalO("grupo6@bioinformatica.com", filename, outfilename=out_file, fmt='clustal')
    iqtree_exe = os.path.join(project_dir, 'resources/iqtree-1.6.12-Windows/bin/iqtree.exe')
    dirName = 'IQtree_results-' + nameFile
    try:
        os.mkdir(dirName)
        print("Directory ", dirName, " Created ")
    except FileExistsError:
        print("Directory ", dirName, " already exists")
    cwd = os.getcwd()
    shutil.move(cwd + '/' + out_file, cwd + '/' + dirName)
    os.chdir(cwd + '/' + dirName)
    AlignIO.convert(out_file, "clustal", "egfr-family.phy", "phylip-relaxed")
    os.remove(out_file)
    # IQTree need 3 sequences
    if quantitySequences > 3:
        subprocess.run([iqtree_exe, '-s', "egfr-family.phy", '-bb', str(bootstrap), '-redo'])
    else:
        subprocess.run([iqtree_exe, '-s', "egfr-family.phy", '-redo'])
    os.remove('egfr-family.phy')
    seqTree = open("egfr-family.phy.treefile", "r")
    tree_phy = Tree(str(seqTree.readlines().__getitem__(0)))
    seqTree.close()
    # os.remove("egfr-family.phy.bionj")
    # os.remove("egfr-family.phy.ckp.gz")
    # os.remove("egfr-family.phy.mldist")
    return tree_phy


