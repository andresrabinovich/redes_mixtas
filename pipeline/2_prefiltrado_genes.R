######################################################################
#SCRIPT QUE FILTRA DATOS DE EXPRESION
#Autor: Andrés Rabinovich en base a un script de Ariel Chernomoretz
#Creación: 01/06/2018
#Última modificación: XX/XX/XXXX (por XXX)
######################################################################

#Librerías que necesita
require(edgeR)

#Elige directorio de trabajo
setwd("/home/arabinov/doctorado/programacion/redes_mixtas/")

#Levanta las cuentas 
(load("pipeline_archivos/cuentas.Rdata"))

#Cargamos la lista de reguladores generada a partir de varios repositorios de datos.
(load("pipeline_archivos/reguladores.Rdata"))

#Elegimos las condiciones para la que vamos a generar la red
(load("pipeline_archivos/1_seleccion_de_condiciones.Rdata"))

#Nos quedamos con los reguladores de expresión génica
reguladores_de_expresion <- union(union(poi[["TF"]],poi[["TCF"]]),reguladores[reguladores[,"tipo_de_regulador"]=="TF","gene_id"])

#Solo queremos las cuentas de la temperatura indicada
iTemp <- grep(paste0("at_",temperatura,"_"),colnames(cuentas_genes))
iTemp <- c(iTemp, grep(paste0("at_",temperatura_referencia,"_"),colnames(cuentas_genes)))

#Filtra las cuentas. Un gen se considera expresado si recibió en promedio más de un mínimo número de cuentas
#(parametro minGenReads, default=10) y si el average read density per condition supera un umbral en alguna de
#las condiciones (parametro minRds,default=0.05).
i1    <- apply(cuentas_genes[,iTemp],1,function(x){mean(x) > 10})
i22   <- apply(cuentas_genes[,iTemp]/cuentas_genes[,"effective_length"],1,function(x){mean(x) > 0.05})
ipass <- i1 & i22

genes_de_interes <- unique(c(unlist(poi), reguladores$gene_id))

#AGREGADO POR AR: Reduce la red únicamente a los genes de interes (reguladores))
ipass <- names(ipass) %in% genes_de_interes

#Filtramos las cuentas
cuentas_genes <- cuentas_genes[ipass, c(1:9, iTemp)]

#data.frame con las condiciones experimentales.
phenotype<-data.frame(condition    = c(paste0(temperatura_referencia, ".T", rep(seq(1:12), each=2)), 
                                       paste0(temperatura, ".T", rep(seq(1:12), each=2)) ), 
                      temperaturas = rep(c(temperatura_referencia, temperatura), each=24),
                      tiempos      = rep(1:12, each=2, times=2)
)


#Ordena los niveles de las condiciones experimentales
phenotype$condition    <- factor(phenotype$condition,levels=unique(phenotype$condition)) 
phenotype$tiempos      <- factor(phenotype$tiempos,levels=unique(phenotype$tiempos)) 
phenotype$temperaturas <- factor(phenotype$temperaturas,levels=unique(phenotype$temperaturas)) 
rownames(phenotype) <- paste0("at_", paste(rep(c(temperatura_referencia, temperatura), each=24), 
                                           paste(rep(1:12, each=2), c("A", "B"), sep="_"), 
                                           sep="_"))

#Arma el objeto DGEList con las cuentas de los bines y de los genes y las condiciones
y <- DGEList(counts=cuentas_genes[,10:ncol(cuentas_genes)], group = phenotype$condition, genes = cuentas_genes[,1:9])

#Calcula el factor de normalización para escalar los tamaños de las librerías
y <- calcNormFactors(y, method=c("TMM","RLE","upperquartile")[1])

#El diseño experimental es cada condición contra la referencia, que elegimos por defecto sea T22.
#Entonces, para decidir si un bin se expresó diferencialmente, se compara su expresión en cada tiempo contra la 
#del mismo tiempo en T22
design <- model.matrix(~condition + 0, data=phenotype)

#Para fitear una binomial negativa se necesita estimar el parámetro de dispersión de la binomial negativa (1/r, con r 
#la cantidad de veces que tiene que fallar la binomial negativa).
y      <- estimateDisp(y, design)

#Se ajusta con el modelo de binomial negativa las cuentas de los bines y se los compara con la referencia.
#En realidad el ajuste devuelve los beta que mejor ajustan y como es un glm los beta son logs. La resta del beta contra
#el beta de la referencia entonces es un ratio entre betas y eso es lo que termina siendo la comparación.
#Además, modificamos la función glmFit para que devuelva el error estandar. Necesitamos el error estandar para poder
#reducir los coeficientes que vienen de comparaciones con muy pocas cuentas y que por eso se disparan como si estuvieran
#cambiando mucho pero en realidad no cambian sino que están en cero y es el gen el que cambia.
fit    <- glmFit(y, design)

#Ajustamos por glmLRT para encontrar que genes se encuentran diferencialmente expresados entre las dos condiciones
lfchange_genes           <- matrix(NA,ncol=(0.5*ncol(fit$coef)),nrow=nrow(fit$coef))
rownames(lfchange_genes) <- rownames(fit$coef)
colnames(lfchange_genes) <- paste("condicion", paste(temperatura_referencia,temperatura, sep="."), 1:12, sep=".")
pvalues_genes            <- lfchange_genes

for(i in 1:(0.5*ncol(fit$coef))){

  #Generamos los contrastes
  contrastes         <- rep(0, 24)
  contrastes[i]      <- -1
  contrastes[i + 12] <- 1
  cat("===============================================\nEvaluando contrastes",
      paste(paste0(contrastes[1:12], collapse = ""), paste0(contrastes[13:24], collapse = ""), sep="|")
      ,"\n")
  
  ds                         <- glmLRT(fit, contrast = contrastes)
  genname                    <- rownames(ds$genes)
  lfchange_genes[genname, i] <- ds$table$logFC
  pvalues_genes[genname, i]  <- ds$table$PValue
}
lfchange_genes <- lfchange_genes[!apply(lfchange_genes, 1, function(x){any(is.na(x))}), ]
pvalues_genes  <- pvalues_genes[!apply(pvalues_genes, 1, function(x){any(is.na(x))}), ]

#ajusto TODOS los pv 
qvalues_genes           <- matrix(p.adjust(pvalues_genes,"fdr"), ncol=ncol(pvalues_genes), byrow=FALSE)        
rownames(qvalues_genes) <- rownames(pvalues_genes)


#data.frame con las condiciones experimentales.
phenotype<-data.frame(condition    = c(paste0("T", temperatura, ".t", rep(seq(1:12), each=2))), 
                      temperaturas = rep(c(temperatura), each=24),
                      tiempos      = rep(1:12, each=2)
)


#Ordena los niveles de las condiciones experimentales
phenotype$condition    <- factor(phenotype$condition,levels=unique(phenotype$condition)) 
phenotype$tiempos      <- factor(phenotype$tiempos,levels=unique(phenotype$tiempos)) 
phenotype$temperaturas <- factor(phenotype$temperaturas,levels=unique(phenotype$temperaturas)) 
rownames(phenotype) <- paste0("at_", paste(rep(c(temperatura), each=24), 
                                           paste(rep(1:12, each=2), c("A", "B"), sep="_"), 
                                           sep="_"))

genes_crudos   <- cuentas_genes[, grep(paste0("at_", temperatura), colnames(cuentas_genes))]

#Fiteamos las cuentas de los genes por condición para temperatura
design <- model.matrix(~condition+0, data=phenotype)

#Arma el objeto DGEList con las cuentas de los bines y de los genes y las condiciones
y_genes <- DGEList(counts=genes_crudos, group = phenotype$condition, genes = genes_crudos[,1:8])

#Calcula el factor de normalización para escalar los tamaños de las librerías
y_genes <- calcNormFactors(y_genes, method=c("TMM","RLE","upperquartile")[1])

#Para fitear una binomial negativa se necesita estimar el parámetro de dispersión de la binomial negativa (1/r, con r 
#la cantidad de veces que tiene que fallar la binomial negativa).
y_genes <- estimateDisp(y_genes, design)

#Ajustamos los genes y usamos los coeficientes del ajuste como perfiles
fit_genes <- glmFit(y_genes, design)

#Guarda las cuentas de los genes y los cambios (log fold change) entre condiciones
perfiles_genes <- fit_genes$coefficients
save(lfchange_genes, qvalues_genes, perfiles_genes, cuentas_genes, reguladores_de_expresion, file="pipeline_archivos/2_genes_prefiltrados.Rdata")
