import subprocess
import threading
from pathlib import Path
from tkinter import ttk, filedialog, Tk, IntVar, Entry, StringVar, OptionMenu, Label, Frame
import os
from Bio import SeqIO
from Bio.Alphabet import generic_protein, generic_nucleotide
from src import geolocation, evolutionaryInference
from src.evolutionaryInference import clustalo
import matplotlib.pyplot as plt

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

SEQUENCE_TYPE = {"Nucleotido": generic_nucleotide, "Proteina": generic_protein}


class Root(Tk):
    def __init__(self):
        super(Root, self).__init__()
        self.title('Biolocalizacion')
        self.minsize(500, 200)
        self.wm_iconbitmap('../resources/favicon.ico')
        self.labelFrameWaiting = ttk.LabelFrame(self, text='')
        self.labelFrameWaiting.grid(column=0, row=0, padx=30, pady=10)
        self.labelFrameOptionMenu = ttk.LabelFrame(self, text='Tipo de secuencia')
        self.labelFrameOptionMenu.grid(column=0, row=1, padx=20, pady=20)
        self.optionMenu()
        self.labelFrameOpenFile = ttk.LabelFrame(self, text='Abrir archivo a geolocalizar')
        self.labelFrameOpenFile.grid(column=1, row=1, padx=20, pady=20)
        self.button()
        self.labelFrameBootstrap = ttk.LabelFrame(self, text='Bootstrap')
        self.labelFrameBootstrap.grid(column=2, row=1, padx=20, pady=20)
        self.input()

    def optionMenu(self):
        self.typeSequence = StringVar(self.labelFrameOptionMenu)
        self.typeSequence.set(next(iter(SEQUENCE_TYPE.keys())))
        w = OptionMenu(self.labelFrameOptionMenu, self.typeSequence, *SEQUENCE_TYPE.keys())
        w.pack()

    def button(self):
        self.button = ttk.Button(self.labelFrameOpenFile, text='Ingrese archivo', command=self.fileDialog)
        self.button.grid(column=1, row=1, padx=30, pady=10)

    def fileDialog(self):
        self.fileName = filedialog.askopenfilename(initialdir='/', title='Seleccionar archivo',
                                                   filetype=(('fasta', '*.fasta'), ('All Files', '*.*')))
        try:
            self.check_fasta()
            self.button.configure(text=os.path.basename(self.fileName))
            self.waitingLabel()
            geolocation.dataset(self.fileName, alphabet=SEQUENCE_TYPE.get(self.typeSequence.get()),
                                bootstrap=self.bootstrap.get(), aligned=self.aligned)
            # evolutionaryInference.fasta_to_tree(self.fileName) # TEST
            self.waitingLabel.config(text="Procesado con exito")
            self.update()
            threading.Thread(target=lambda: os.system('egfr-family.phy.log')).start()
            threading.Thread(target=lambda: plt.show()).run()
        except FileNotFoundError:
            None
        except Exception:
            None

    def check_fasta(self):
        self.aligned = False
        self.label = ttk.Label(self.labelFrameOpenFile, text='')
        self.label.grid(column=1, row=2)
        self.is_fasta()
        self.is_content_valid_fasta()

    def is_fasta(self):
        with open(self.fileName, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            if not any(fasta):  # False when `fasta` is empty, i.e. wasn't a FASTA file
                msg = 'El formato del archivo debe ser fasta'
                self.label.configure(text=msg)
                raise Exception(msg)

    def is_content_valid_fasta(self):
        with open(self.fileName, "r") as handle:
            lines = iter(handle.read().splitlines(True))
            nextLine = next(lines, None)
            seqs = []
            seq = []
            if nextLine.startswith('>'):
                nextLine = next(lines, None)
                while nextLine is not None:
                    if not nextLine.startswith('>'):
                        seq.append(nextLine)
                        nextLine = next(lines, None)
                    else:
                        seq = "".join(seq)
                        seq = seq.replace("\n", "")
                        seq = seq.replace("\r", "")
                        seqs.append(seq)
                        seq = []
                        nextLine = next(lines, None)
                seq = "".join(seq)
                seq = seq.replace("\n", "")
                seq = seq.replace("\r", "")
                seqs.append(seq)
                if len(set(list(map(lambda sequence: len(sequence), seqs))))==1:
                    self.aligned = True
                    raise Exception("Secuencias alineadas, se salteara ese paso")
            else:
                msg = 'El formato del contenido debe comenzar con > para cada header por secuencia'
                self.label.configure(text=msg)
                raise Exception(msg)

    def waitingLabel(self):
        self.waitingLabel = ttk.Label(self.labelFrameWaiting, text='')
        self.labelFrameOpenFile.destroy()
        self.labelFrameBootstrap.destroy()
        self.labelFrameOptionMenu.destroy()
        self.waitingLabel.grid(column=0, row=0)
        self.waitingLabel.config(text="Procesando...")
        self.update()

    def input(self):
        vcmd = (self.register(self.onValidate), '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        self.bootstrap = IntVar(self, value="1000")
        self.entry = Entry(self.labelFrameBootstrap, validate="key", validatecommand=vcmd, textvariable=self.bootstrap)
        self.entry.grid(column=1, row=1, padx=30, pady=10)
        self.entry.pack()

    def onValidate(self, action, index, value_if_allowed,
                   prior_value, text, validation_type, trigger_type, widget_name):
        if value_if_allowed:
            try:
                int(value_if_allowed)
                return True
            except ValueError:
                return False
        else:
            return False


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
