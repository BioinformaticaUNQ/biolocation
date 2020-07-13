my_file = "sarasa.fasta"
import os
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

limiteDeSecuencias = 3


def check_is_fasta(filename):
    ext = os.path.splitext(filename)[-1].lower()
    # Now we can simply use == to check for equality, no need for wildcards.
    if ext != ".fasta":
        raise Exception(
            'el formato del archivo debe ser fasta, pedazo de forre')


def valid(filename):
    check_is_fasta(filename)
    check_contenido(filename)


def check_contenido(filename):
    with open(filename) as handle:
        lines = handle.read().splitlines(True)

    lines = eliminarEspaciosEnBlanco(lines)
    largoDePrimerSecuencia = len(lines[0])

    lines = list(filter(lambda seq: comienzaConCaracter(seq), lines))
    try:
        checkearCantidadDeSecuencias(lines)
    except:
        print("se cortara la lista de secuencias --- debe abrir una ventana para mostrar el msj")
        lines = lines[:limiteDeSecuencias]

    comprobarSiSonTodosDelMismoTipo(lines)

    secuenciasTienenMismoLargo(lines,largoDePrimerSecuencia)



def checkearCantidadDeSecuencias(lines):
    if (len(lines) > limiteDeSecuencias):
        textoParseado = 'solo se van a utilizar las primeras {0} secuencias'.format(
            limiteDeSecuencias)
        raise Exception(textoParseado)


def comienzaConCaracter(lines):
    return any((line.startswith(">")) for line in lines)


def eliminarEspaciosEnBlanco(lines):
    return [line.replace("\n", "").replace("\r", "") for line in lines]


def secuenciasTienenMismoLargo(lines, largoDePrimerSecuencia):
    if any((len(line) != largoDePrimerSecuencia) for line in lines):
      raise Exception("Secuencias no alineadas, desea alinearlas?") ##deberia preguntar esto en un modal y al confirmar deberia llamar a algo que haga la alineacion

def comprobarSiSonTodosDelMismoTipo(lines):
  tipo= "proteina"
  if all (("gta" in seq or "ta" in seq) for seq in lines) :
    tipo= "nucleotido"
  else:
    raise Exception("Las secuencias ingresadas no son del mismo tipo")




if __name__ == "__main__":
    valid(my_file)
