import platform
import shutil
import subprocess
import threading
from importlib import reload
from pathlib import Path
from tkinter import ttk, filedialog, Tk, IntVar, Entry, StringVar, OptionMenu, Label, Frame, messagebox
import os
from Bio import SeqIO
from Bio.Alphabet import generic_protein, generic_nucleotide
from src import geolocation
import matplotlib.pyplot as plt
import python_script

PROJECT_DIR = Path(__file__).parent.parent

SEQUENCE_TYPE = {"Nucleotido": generic_nucleotide, "Proteina": generic_protein}

BOOTSTRAP_QUANTITY = [1000, 2000, 3000]

SEQUENCES_LIMIT = 200

SEQUENCES_MIN = 5


class TypeSequenceException(Exception):
    pass


class Root(Tk):
    def __init__(self):
        super(Root, self).__init__()
        self.title('Biolocalizacion')
        self.minsize(500, 200)
        ico = os.path.join(PROJECT_DIR, 'resources/favicon.ico')
        # self.wm_iconbitmap(ico)
        self.labelFrameWaiting = ttk.LabelFrame(self, text='')
        self.labelFrameWaiting.grid(column=0, row=1, padx=30, pady=10)
        self.createLabel()

    def createLabel(self):
        self.labelFrameOptionMenu = ttk.LabelFrame(self, text='Tipo de secuencia')
        self.labelFrameOptionMenu.grid(column=0, row=0, padx=20, pady=20)
        self.optionMenu()
        self.labelFrameOpenFile = ttk.LabelFrame(self, text='Abrir archivo a geolocalizar')
        self.labelFrameOpenFile.grid(column=1, row=0, padx=20, pady=20)
        self.label = ttk.Label(self.labelFrameOpenFile, text='')
        self.label.grid(column=1, row=2)
        self.button()
        self.labelFrameBootstrap = ttk.LabelFrame(self, text='Bootstrap')
        self.labelFrameBootstrap.grid(column=2, row=0, padx=20, pady=20)
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
        if platform.system() == 'Windows':  # Windows
            self.fileName = filedialog.askopenfilename(initialdir='/', title='Seleccionar archivo',
                                                       filetype=(('fasta', '*.fasta'), ('All Files', '*.*')))
        else:  # linux variants
            self.fileName = filedialog.askopenfilename(initialdir='/', title='Seleccionar archivo',
                                                       filetypes=(('fasta', '*.fasta'), ('All Files', '*.*')))
        try:
            self.check_fasta()
            self.button.configure(text=os.path.basename(self.fileName))
            self.waitingLabel()
            geolocation.run(self.fileName, alphabet=SEQUENCE_TYPE.get(self.typeSequence.get()),
                            bootstrap=self.bootstrap.get(), aligned=self.aligned,
                            quantitySequences=self.quantitySequences)
            self.waitingLabel.config(text='Procesado con exito')
            self.update()
            filepath = 'egfr-family.phy.log'
            # Windows
            # threading.Thread(target=lambda: os.system('egfr-family.phy.log')).start()
            # Linux
            # threading.Thread(target=lambda: subprocess.run(["xdg-open", 'egfr-family.phy.log'], check=True)).start()
            if platform.system() == 'Darwin':  # macOS
                threading.Thread(target=lambda: subprocess.run(['open', filepath])).start()
            elif platform.system() == 'Windows':  # Windows
                threading.Thread(target=lambda: os.startfile(filepath)).start()
            else:  # linux variants
                threading.Thread(target=lambda: subprocess.run(["xdg-open", filepath], check=True)).start()
            threading.Thread(target=lambda: plt.show()).run()
            self._update()
        except FileNotFoundError:
            pass
        except TypeSequenceException:
            pass
        except Exception:
            pass

    def check_fasta(self):
        self.aligned = False
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
                self.quantitySequences = len(seqs)
                if self.quantitySequences < SEQUENCES_MIN:
                    msg = 'Un mínimo de 5 secuencias es requerido'
                    self.label.configure(text=msg)
                    raise Exception(msg)
                msg = ''
                if self.is_aligned(seqs):
                    self.aligned = True
                    msg = msg + ' Secuencias alineadas, se salteara ese paso'
                self.label.configure(text=msg)
                if len(set(list(map(lambda sequence: len(list(set(sequence))), seqs)))) == 1:
                    for sequence in seqs:
                        seqDistinct = list(set(sequence))
                        seqDistinct.sort()
                        if seqDistinct == ['A', 'C', 'G', 'T'] or seqDistinct == ['A', 'C', 'G', 'U']:
                            if self.typeSequence.get().lower() != "nucleotido":
                                msg = 'Ha especificado tipo de secuencia: ' + self.typeSequence.get() + ' y son posibles secuencias de nucleotidos'
                                if not messagebox.askyesno(message=msg + ", ¿Desea continuar?", title="Advertencia"):
                                    raise TypeSequenceException(msg)
                                else:
                                    break
                        else:
                            if self.typeSequence.get().lower() != 'proteina':
                                msg = 'Ha especificado tipo de secuencia: ' + self.typeSequence.get() + ' y son posibles secuencias de proteinas'
                                if not messagebox.askyesno(message=msg + ", ¿Desea continuar?", title="Advertencia"):
                                    raise TypeSequenceException(msg)
                                else:
                                    break
            else:
                msg = 'El formato del contenido debe comenzar con > para cada header por secuencia'
                self.label.configure(text=msg)
                raise Exception(msg)

    def is_aligned(self, seqs):
        return len(set(list(map(lambda sequence: len(sequence), seqs)))) == 1 and next(
            filter(lambda sequence: '-' in sequence, seqs), None) is not None

    def waitingLabel(self):
        self.waitingLabel = ttk.Label(self.labelFrameWaiting, text='')
        self.infoWaitingLabel = ttk.Label(self.labelFrameWaiting, text='')
        self.button.state(["disabled"])
        self.waitingLabel.grid(column=1, row=0)
        self.waitingLabel.config(text="Procesando...")
        self.update()

    def input(self):
        self.bootstrap = IntVar(self.labelFrameBootstrap)
        self.bootstrap.set(next(iter(BOOTSTRAP_QUANTITY)))
        w = OptionMenu(self.labelFrameBootstrap, self.bootstrap, *BOOTSTRAP_QUANTITY)
        w.pack()

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

    def _update(self):
        python_script.main_refresh(self, python_script)


def main_refresh(root, python_script):
    reload(python_script)
    root.destroy()
    python_script.main()


def main():
    root = Root()
    root.mainloop()


if __name__ == '__main__':
    main()
