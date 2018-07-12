######################################################################
#SCRIPT PARA ELEGIR LAS CONDICIONES PARA LA GENERACION DE LA RED MIXTA
#Autor: Andrés Rabinovich
#Creación: 01/06/2018
#Última modificación: XX/XX/XXXX (por XXX)
######################################################################

#Elige directorio de trabajo
setwd("/home/arabinov/doctorado/programacion/redes_mixtas/pipeline_archivos/")

#Elegimos las condiciones para la que vamos a generar la red
temperatura            = 12
temperatura_referencia = 22

#Guarda las condiciones para ser usadas en el resto de los scripts
save(temperatura, temperatura_referencia, file="1_seleccion_de_condiciones.Rdata")