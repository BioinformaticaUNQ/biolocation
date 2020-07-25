# Biolocalización

## Introducción
Este programa está realizado en Python. Su función es graficar árboles filogénicos mediante archivos fastas de nucleotidos o de proteínas con formato de cabecera por secuencia: '>nroGenbank' (signo mayor + número Genbank). En la carpeta 'resources' hay ejemplos de archivos. Si el archivo no se encuentra alineado, se alineará y se guardará el archivo alineado para su próxima utilización. Además,se cuenta con un mapa que nos ubica las secuencias, dependiendo si poseen un país de origen consultandolo con el número de Genbank a la base de datos del NCBI Entrez. 

## Requisitos previos
```
Anaconda por ejemplo wget https://repo.anaconda.com/archive/Anaconda3-5.3.1-Linux-x86_64.sh
bash Anaconda-5.3.1-Linux-x86_64.sh
Fix instalación anaconda en Linux https://www.youtube.com/watch?v=DY0DB_NwEu0:
 -gedit .bashrc
 -export PATH="/home/{user_name}/anaconda3/bin:$PATH"
 -source .bashrc
```

## Instalación y dependecias
`conda env create -f environment.yml`

Si el comando no trabaja, instalar las depencias manualmente:
```
conda create --name bioinformatica
conda update -n base -c defaults conda
conda activate bioinformatica
conda install git
sudo apt install python3-pip
pip3 install --upgrade pip
pip3 install --upgrade pip xmltramp2 requests
pip install biopython
pip install geopy
pip install matplotlib
pip install ete3
conda install basemap
```

 
 **Observación:** La libreria "ete3" puede traer conflicto a la hora de utilizarse, el cual podría resolverse de las siguientes maneras:
    - Instalar/actualizar la versión de matplotlib
    - Realizar la instalación mediante la consola de Anaconda con el comando ```conda install -c etetooltikit ete3```.
    

### Actualizar conda env
`conda env update -n bioinformatica --file environment.yml  --prune`


## Utilización:

Al correr el proyecto nos aparecerá una pequeña ventana en la cuál, podremos visualizar dos selectores al costado del botón para subir un archivo. El primero es para elegir el tipo de secuencia a utilizar (estas podrán ser nucleotidos o proteínas), el segundo es para elejir el Bootstrap que se querra utilizar para graficar el arbol.
 El tipo de archivo admitido es fasta, el cuál debe respetar el formato universalmente utilizado y deberá tener al mínimo 5 y como máximo 100 secuencias

### Bash command line
`python python_script.py`


![bioIf3 (1)](https://user-images.githubusercontent.com/31372437/88124537-be3d0800-cba3-11ea-9319-c5dda01de6cc.jpg)
 
 En el caso de elegir el tipo incorrecto de secuencia, nos aparecerá una advertencia la cuál nos indicará el error y con podrémos decidir si deseamos continuar o no (al seleccionar el dato erroneo no perderemos funcionalidad con respecto al arbol pero el mapa será visualizado sin datos)
 
 El proceso puede tardar unos minutos (todo depende de la cantidad de secuencias ingresadas y del hardware de la PC), mientras se va ejecutando se podrá visualizar en la consola la ejecución.
 
 Tambíen se generaran en el proyecto los archivos correspondientes a la alineacion, el log del arbol, la imagen en formato ".png" del arbol y la imagen en formato ".png" del mapa.
 
![bioIf2 (1)](https://user-images.githubusercontent.com/31372437/88124471-9d74b280-cba3-11ea-884d-18ff7a9dcee0.jpg)

 
 **Observación:** Si las secuencias ingresadas no poseen país de origen dentro  de los registros de Entrez, no se mostrarán en el mapa.
 
 Finalizado el proceso, contaremos con 3 ventanas abiertas:
    - El log del IQTree
    - El árbol filogénico
    - El mapa con la ubicación de las secuencias

 **Observación:** Los colores de los nodos no ancestrales coinciden con los puntos en el mapa para su identificación.
  
![bioIn (1)](https://user-images.githubusercontent.com/31372437/88124418-761de580-cba3-11ea-8fbe-d914f7817dc5.jpg)
 
 ## Interacción con el árbol:

![archivo100 (1)](https://user-images.githubusercontent.com/31372437/88124751-34416f00-cba4-11ea-84a0-46a09ef89ded.jpg)

 En el caso de contar con muchas secuencias, el no se podrá visualizar de manera óptima las partes que componen al árbol. Por lo cuál, será necesario hacerle zoom. La utilización de la herramienta zoom no es trivial, ya que primero se debe seleccionar con el cursor la parte que se quiere visualizar y luego tocar el icono correspondiente. 
 
![mapa2de100 (1)](https://user-images.githubusercontent.com/31372437/88124613-eb89b600-cba3-11ea-902d-09ef9f4296ec.jpg)
![mapa100de3 (1)](https://user-images.githubusercontent.com/31372437/88124680-0eb46580-cba4-11ea-928a-eab4c68ca73c.jpg)
