#Elige directorio de trabajo
setwd("/home/arabinov/doctorado/programacion/redes_cosplicing/pipeline_archivos/")

#Elegimos las condiciones para la que vamos a generar la red
temperatura            = 12
temperatura_referencia = 22

#Guarda las condiciones para ser usadas en el resto de los scripts
save(temperatura, temperatura_referencia, file="1_seleccion_de_condiciones.Rdata")