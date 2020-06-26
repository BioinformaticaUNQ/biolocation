from tkinter import ttk, filedialog, Tk
import os
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
        if self.fileName.endswith('.fasta'):
            self.button.configure(text=os.path.basename(self.fileName))

            # geolocation.dataset(self.fileName)

            # out_file = os.path.basename(self.fileName) + '-aligned.fasta'
            # clustalo.runClustalO("martin@gmail.com", self.fileName, outfilename=out_file, fmt='fasta')
            evolutionaryInference.fasta_to_tree(self.fileName)
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