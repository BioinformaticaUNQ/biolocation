import os
import sys
from pathlib import Path

from Bio import SeqIO, Phylo
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo import PhyloXML
from Bio.Phylo.Applications import PhymlCommandline
from Bio.Phylo.PhyloXML import Phylogeny
from io import StringIO
from Bio import AlignIO


def fasta_to_tree(filename):
    file = open(filename, "r")
    records = [seqrec for seqrec in SeqIO.parse(file, "fasta")]
    project_dir = Path(__file__).parent.parent.parent
    muscle_exe = os.path.join(project_dir, "resources/muscle3.8.31_i86/muscle3.8.31_i86win32.exe")
    # muscle_exe = os.path.join(project_dir, "resources/muscle3.8.31_i86/muscle3.8.31_i86linux32")
    muscle_cline = MuscleCommandline(muscle_exe, input=filename, clwstrict=True)
    stdout, stderr = muscle_cline()
    # align = AlignIO.read(StringIO(stdout), "clustal")
    # align = AlignIO.read(StringIO(stdout), "fasta")
    AlignIO.convert(StringIO(stdout), "clustal",  "egfr-family.phy", "phylip-relaxed")
    # alignments = AlignIO.parse(StringIO(stdout), format="clustal", alphabet="phylip-relaxed")
    phyml_exe = os.path.join(project_dir, "resources/PhyML-3.1/PhyML-3.1_win32.exe")
    # phyml_exe = os.path.join(project_dir, "resources/PhyML-3.1/PhyML-3.1_linux32")
    phyml_cline = PhymlCommandline(phyml_exe, input='egfr-family.phy', datatype='aa', model='WAG', alpha='e',
                                   bootstrap=1)
    out_log, err_log = phyml_cline()
    os.remove('egfr-family.phy')
    egfr_tree = Phylo.read("egfr-family.phy_phyml_tree.txt", "newick")
    os.remove('egfr-family.phy_phyml_tree.txt')
    # Phylo.draw_ascii(egfr_tree)
    egfr_tree.rooted = True
    # Phylo.draw(egfr_tree)
    # Promote the basic tree to PhyloXML
    egfr_phy = egfr_tree.as_phyloxml()
    egfr_phy = Phylogeny.from_tree(egfr_phy)
    egfr_phy.root.color = (128, 128, 128)
    egfr_phy.root.color = "#808080"  # This is one alternative
    egfr_phy.root.color = "gray"  # This is another
    # mrca = egfr_phy.common_ancestor({"name": "E"}, {"name": "F"})
    # mrca.color = "salmon"
    egfr_phy.clade[0].color = "blue"
    Phylo.draw(egfr_phy)
    os.remove('egfr-family.phy_phyml_boot_stats.txt')
    os.remove('egfr-family.phy_phyml_boot_trees.txt')
    os.remove('egfr-family.phy_phyml_stats.txt')
    # Make a lookup table for sequences
    '''lookup = dict((rec.id, str(rec.seq)) for rec in records)
    for clade in egfr_phy.get_terminals():
        key = clade.name
        accession = PhyloXML.Accession(key, 'NCBI')
        mol_seq = PhyloXML.MolSeq(lookup[key], is_aligned=True)
        sequence = PhyloXML.Sequence(type='aa', accession=accession, mol_seq=mol_seq)
        clade.sequences.append(sequence)
    # Save the annotated phyloXML file
    Phylo.write(egfr_phy, 'egfr-family-annotated.xml', 'phyloxml')'''
