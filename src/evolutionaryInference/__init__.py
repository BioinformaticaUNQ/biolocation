import os
import shutil
import subprocess
from pathlib import Path
from Bio import AlignIO
from ete3 import Tree
from src.evolutionaryInference import clustalo

project_dir = Path(__file__).parent.parent.parent


def fasta_to_tree(filename, aligned, bootstrap, quantitySequences):
    os.chdir(project_dir)
    nameFile = os.path.basename(filename).split('.')[0]
    out_file = filename
    if not aligned:
        out_file = nameFile + '-aligned.fasta'
        clustalo.runClustalO("grupo6@bioinformatica.com", filename, outfilename=out_file, fmt='clustal')
    iqtree_exe = os.path.join(project_dir, 'resources/iqtree-1.6.12-Windows/bin/iqtree.exe')
    dirName = 'IQtree_results-' + nameFile
    cwd = os.getcwd()
    try:
        os.mkdir(dirName)
        print("Directory ", dirName, " Created ")
    except FileExistsError:
        print("Directory ", dirName, " already exists")
        os.remove(cwd + '/' + dirName + '/' + out_file)
    shutil.move(cwd + '/' + out_file, cwd + '/' + dirName)
    os.chdir(cwd + '/' + dirName)
    AlignIO.convert(out_file, "clustal", "egfr-family.phy", "phylip-relaxed")
    # IQTree need 5 sequences
    subprocess.run([iqtree_exe, '-s', "egfr-family.phy", '-bb', str(bootstrap), '-redo'])
    os.remove('egfr-family.phy')
    seqTree = open("egfr-family.phy.treefile", "r")
    tree_phy = Tree(str(seqTree.readlines().__getitem__(0)))
    seqTree.close()
    return tree_phy


