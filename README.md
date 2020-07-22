# bioinformatica

## Introducción
Este programa está realizado en Python. Su función es graficar árboles filogénicos mediante archivos fastas de nucleotidos o de proteínas. Si el archivo no se encuentra alineado, se alineará y se guardará el archivo alineado para su próxima utilización. Además,se cuenta con un mapa que nos ubica las secuencias, dependiendo si poseen un país de origen.

## Instalación y dependecias
### Linux:
```
wget https://repo.anaconda.com/archive/Anaconda3-5.3.1-Linux-x86_64.sh
bash Anaconda-5.3.1-Linux-x86_64.sh
https://www.youtube.com/watch?v=DY0DB_NwEu0
gedit .bashrc
export PATH="/home/{user_name}/anaconda3/bin:$PATH"
source .bashrc
anaconda-navigator
conda create --name   bioinformatica
conda update -n base -c defaults conda
conda activate bioinformatica
conda install git
sudo apt install python3-pip
pip3 install --upgrade pip
conda activate
git clone https://github.com/martinCastello/bioinformatica.git
pip3 install --upgrade pip xmltramp2 requests
pip install biopython
pip install geopy
pip install ete3
pip install matplotlib
conda install basemap
conda env export > environment.yml
```

Para intercalar entre windows y linux:
https://stackoverflow.com/questions/434597/open-document-with-default-os-application-in-python-both-in-windows-and-mac-os

 
 **Observación:** La libreria "et3" puede traer conflicto a la hora de utilizarse, el cual podría resolverse de las siguientes maneras:
    - Instalar/actualizar la versión de matplot
    - Realizar la instalación mediante la consola de Anaconda con el comando ```conda install -c etetooltikit et3```.
    
## Utilización:

Al correr el proyecto en el ide utilizado, nos aparecerá una pequeña ventana en la cuál, podremos visualizar dos selectores al costado del botón para subir un archivo. El primero es para elegir el tipo de secuencia a utilizar (estas podrán ser nucleotidos o proteínas), el segundo es para elejir el Bootstrap que se querra utilizar para graficar.
 El tipo de archivo admitido es fasta, el cuál debe respetar el formato universalmente utilizado y deberá tener al mínimo 8 y como máximo 300 secuencias

![bioIf3](https://user-images.githubusercontent.com/31372437/88123394-2fc78700-cba1-11ea-8c0c-5fcf26a9a80d.jpg)
 
 En el caso de elegir el tipo incorrecto de secuencia, nos aparecerá una advertencia la cuál nos indicará el error y con podrémos decidir si deseamos continuar o no (al seleccionar el dato erroneo no perderemos funcionalidad )
 
 El proceso puede tardar unos minutos (todo depende de la cantidad de secuencias ingresadas y del hardware de nuestra PC), mientras se va ejecutando se podrá visualizar en la consola la ejecución.
 
 Tambíen se generaran en el proyecto los archivos correspondientes a la alineacion, el log del arbol y la imagen en formato ".png" del mapa.
 
 ![bioIf2](https://user-images.githubusercontent.com/31372437/88123437-47067480-cba1-11ea-9696-5da6dac84f78.jpg)
 
 **Observación:** Si las secuencias ingresadas no poseen país de origen dentro  de los registros de Entrez, no se mostrarán en el mapa.
 
 Finalizado el proceso, contaremos con 3 ventanas abiertas:
    - El log del IQTree
    - El árbol filogénico
    - El mapa con las secuencias

 **Observación:** Los colores de los nodos no ancestrales coinciden con los puntos en el mapa para su identificación.
  
![bioIn](https://user-images.githubusercontent.com/31372437/88123572-87fe8900-cba1-11ea-97be-9b910d588cbc.jpg)
 
 ## Interacción con el árbol:

![archivo100](https://user-images.githubusercontent.com/31372437/88123672-c6944380-cba1-11ea-8d49-9c2fb45d2c8b.jpg)

 En el caso de contar con muchas secuencias, el no se podrá visualizar de manera óptima las partes que componen al árbol. Por lo cuál, será necesario hacerle zoom. La utilización de la herramienta zoom no es trivial, ya que primero se debe seleccionar con el cursor la parte que se quiere visualizar y luego tocar el icono correspondiente. 
 
![mapa2de100](https://user-images.githubusercontent.com/31372437/88123712-dca20400-cba1-11ea-98f4-6ce43569eb9f.jpg)
![mapa100de3](https://user-images.githubusercontent.com/31372437/88123714-ddd33100-cba1-11ea-84d8-5f5b23750566.jpg)
