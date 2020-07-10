import subprocess
from pathlib import Path
from tkinter import ttk, filedialog, Tk
import os

from Bio import SeqIO
from Bio.Alphabet import generic_protein, generic_nucleotide, generic_dna

from src import geolocation, evolutionaryInference
from src.evolutionaryInference import clustalo


###########################################################################
# Cambiar arcchivo py por pyw para poder ejecutarlo directamente del archivo
# raiz=Tk()
# raiz.title('Ingrese archivo .fasta')
# raiz.resizable(True, False)
# .ico en source para cambiar la imagen de la pluma
# raiz.iconbitmap('...ico')
# raiz.geometry("650x350")
# raiz.config(bg='red')
# miFrame=Frame()
# miFrame.pack(fill='both', expand='False')
# miFrame.config(bg='blue')
# miFrame.config(width='650', height='350')
# miFrame.config(bd='35')
# miFrame.config(relief='groove')
# miFrame.config(cursor='hand2')
# raiz.mainloop()
#########################################################################

def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

class Root(Tk):
    def __init__(self):
        super(Root, self).__init__()
        self.title('Biolocalizacion')
        self.minsize(500, 200)
        self.wm_iconbitmap('../resources/favicon.ico')
        self.labelFrame = ttk.LabelFrame(self, text='Abrir archivo a geolocalizar')
        self.labelFrame.grid(column=0, row=1, padx=20, pady=20)
        self.button()

    def button(self):
        self.button = ttk.Button(self.labelFrame, text='Ingrese archivo', command=self.fileDialog)
        self.button.grid(column=1, row=1, padx=30, pady=10)

    def fileDialog(self):
        self.fileName = filedialog.askopenfilename(initialdir='/', title='Seleccionar archivo',
                                                   filetype=(('fasta', '*.fasta'), ('All Files', '*.*')))
        if is_fasta(self.fileName):
        # if self.fileName.endswith('.fasta'):
            self.button.configure(text=os.path.basename(self.fileName))
            geolocation.dataset(self.fileName, generic_dna)
            # evolutionaryInference.fasta_to_tree(self.fileName) # TEST
        else:
            self.label = ttk.Label(self.labelFrame, text='')
            self.label.grid(column=1, row=2)
            self.label.configure(text='El formato del archivo es invalido')


if __name__ == '__main__':
    root = Root()
    root.mainloop()

##########################################################################
# Desde la terminal --> python tpFinal .py SecuenciasCytocromoC.fasta
# inFile = sys.argv[1]
# inFile = '../bioinformatica/SecuenciasCytocromoC.fasta'
# filename, file_extension = os.path.splitext(inFile)
# if file_extension.__eq__('.fasta'):
#    a = open('../bioinformatica/SecuenciasCytocromoC.fasta', 'r').read()
# else:
#    print('El formato del archivo es invalido')
