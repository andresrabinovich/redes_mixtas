##############################################################################################
#SCRIPT QUE GENERA UNA RED DE COSPLICING PARA UNA DADA TEMPERATURA USANDO COMO REFERENCIA T22.
#Autor: Andrés Rabinovich en base al script red_de_cosplicing.R de Ariel Chernomoretz
#Creación: 16/03/2018
#Última modificación: XX/XX/XXXX (por XXX)
##############################################################################################

#Librerías que necesita
library(ASpli)
library(ashr)
setwd("/home/arabinov/doctorado/programacion/redes_cosplicing/pipeline_archivos/")

#Levanta las cuentas 
(load("cuentas.Rdata"))

#Elegimos las condiciones para la que vamos a generar la red
(load("1_seleccion_de_condiciones.Rdata"))

#Solo queremos las cuentas de la temperatura indicada y las de referencia
iTemp<-c(grep(paste0("at_",temperatura_referencia,"_"),colnames(cuentas_genes)), 
         grep(paste0("at_",temperatura,"_"),colnames(cuentas_genes)))
cuentas_genes <- cuentas_genes[, c(1:7, iTemp)]
iTemp<-c(grep(paste0("at_",temperatura_referencia,"_"),colnames(cuentas_bines)), 
         grep(paste0("at_",temperatura,"_"),colnames(cuentas_bines)))
cuentas_bines <- cuentas_bines[, c(1:9, iTemp)]

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

#Aplicamos los filtros para la referencia y para la temperatura indicada. Usamos solo los genes/bines que pasen los
#filtros en ambas condiciones.
df.gen <- cuentas_genes
df.bin <- cuentas_bines
for(temp in c(temperatura_referencia, temperatura)){
  
  #Usamos para filtrar el fenotipo para la temperatura actual
  pheno <- phenotype[phenotype$temperaturas == temp, ]
  
  #Filtra las cuentas de los genes. Un gen se considera expresado si recibió en promedio más de un mínimo número de cuentas
  #(parametro minGenReads, default=10) y si el average read density per condition supera un umbral en alguna de
  #las condiciones (parametro minRds,default=0.05).
  dfG0 <- ASpli:::.filterByReads(df0 = df.gen, targets = pheno, 
                                 min = 10, type = "any", contrast=NULL, onlyContrast = FALSE, 
                                 verbose=FALSE)
  
  dfGen <- ASpli:::.filterByRdGen(df0 = dfG0, targets = pheno, min = 0.05, 
                                  type = "any", contrast=NULL, onlyContrast=FALSE, verbose=FALSE)
  
  #Filtra las cuentas de los bines. Un bin se considera expresado si el gen al que pertenece lo está y si recibió en 
  #promedio más de un mínimo número de cuentas (parametro min, default=5) y con un mínimo de RdBinRATIO
  #(parametro min, default=0.05). Sacamos además los bines que son Intron original que no tienen nada
  #que ver y los external, que meten ruido porque pueden corresponder con otro gen.
  dfBin <- cuentas_bines[df.bin[, "locus"] %in%  row.names(dfGen), ]
  
  dfBin <- dfBin[which(dfBin[, "feature"] != "Io" & dfBin[, "event"] != "external"), ] 
  
  df1 <- ASpli:::.filterByReads(df0 = dfBin, targets = pheno, min = 5, 
                                type = "any", contrast=NULL, onlyContrast=FALSE, verbose=FALSE)
  
  df2 <- ASpli:::.filterByRdBinRATIO(dfBin = df1, dfGen = dfGen, targets = pheno, 
                                     min = 0.05, type = "any", contrast=NULL, onlyContrast = FALSE, 
                                     verbose = FALSE)
  
  #Generamos dos data.frame con la información de los bines y de los genes filtrados
  df.bin <- df2
  df.gen <- dfGen
  rm(dfBin, dfGen, df1, df2);gc()
}

#Cuantos bines y cuantos genes sobreviven el filtro
print(paste("Genes:", nrow(df.gen)))
print(paste("Bines:", nrow(df.bin)))

#Arma el objeto DGEList con las cuentas de los bines y de los genes y las condiciones
y <- DGEList(counts=df.bin[,10:ncol(df.bin)], group = phenotype$condition, genes = df.bin[,1:9])

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
lfchange_bines           <- matrix(NA,ncol=(0.5*ncol(fit$coef)),nrow=nrow(fit$coef))
rownames(lfchange_bines) <- rownames(fit$coef)
colnames(lfchange_bines) <- paste("condicion", paste(temperatura_referencia,temperatura, sep="."), 1:12, sep=".")
pvalues_bines            <- lfchange_bines

for(i in 1:(0.5*ncol(fit$coef))){
  #Generamos los contrastes
  contrastes <- rep(0, 24)
  contrastes[i]    <- -1
  contrastes[i + 12] <- 1
  cat("====================================\nEvaluando contrastes",
      paste(paste0(contrastes[1:12], collapse = ""), paste0(contrastes[13:24], collapse = ""), sep="|")
      ,"\n")
  
  ds       <- diffSpliceDGE(fit, contrast = contrastes, geneid="locus",exonid="exonid")
  
  binname                    <- rownames(ds$genes)
  # Ojo...notar que coeff -> logFC y exon.p.value -> P.Value
  lfchange_bines[binname, i] <-ds$coefficients
  pvalues_bines[binname, i]  <-ds$exon.p.value
  
}
lfchange_bines <- lfchange_bines[!apply(lfchange_bines, 1, function(x){any(is.na(x))}), ]
pvalues_bines  <- pvalues_bines[!apply(pvalues_bines, 1, function(x){any(is.na(x))}), ]

#ajusto TODOS los pv 
qvalues_bines           <- matrix(p.adjust(pvalues_bines,"fdr"), ncol=ncol(pvalues_bines), byrow=FALSE)        
rownames(qvalues_bines) <- rownames(pvalues_bines)

#Genero los perfiles de los bines para la temperatura elegida

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
bines_crudos   <- cuentas_bines[, grep(paste0("at_", temperatura), colnames(cuentas_bines))]

#El diseño experimental es cada condición contra la referencia, que elegimos por defecto sea el primer tiempo.
#Entonces, para decidir si un bin se expresó diferencialmente, se compara su expresión en cada tiempo contra la 
#del tiempo T1
design <- model.matrix(~condition+0, data=phenotype)  #referencia: T1

#Arma el objeto DGEList con las cuentas de los bines y de los genes y las condiciones
y_genes <- DGEList(counts=genes_crudos, group = phenotype$condition, genes = genes_crudos[,1:8])

#Calcula el factor de normalización para escalar los tamaños de las librerías
y_genes <- calcNormFactors(y_genes, method=c("TMM","RLE","upperquartile")[1])

#Para fitear una binomial negativa se necesita estimar el parámetro de dispersión de la binomial negativa (1/r, con r 
#la cantidad de veces que tiene que fallar la binomial negativa).
y_genes <- estimateDisp(y_genes, design)

#Se ajusta con el modelo de binomial negativa las cuentas de los bines y se los compara con la referencia.
#En realidad el ajuste devuelve los beta que mejor ajustan y como es un glm los beta son logs. La resta del beta contra
#el beta de la referencia entonces es un ratio entre betas y eso es lo que termina siendo la comparación.
#Además, modificamos la función glmFit para que devuelva el error estandar. Necesitamos el error estandar para poder
#reducir los coeficientes que vienen de comparaciones con muy pocas cuentas y que por eso se disparan como si estuvieran
#cambiando mucho pero en realidad no cambian sino que están en cero y es el gen el que cambia.
fit_genes <- glmFit(y_genes, design)

#Arma el objeto DGEList con las cuentas de los bines y de los genes y las condiciones
y_bines <- DGEList(counts=bines_crudos, group = phenotype$condition, genes = cuentas_bines[,1:9])

#Calcula el factor de normalización para escalar los tamaños de las librerías
y_bines <- calcNormFactors(y_bines, method=c("TMM","RLE","upperquartile")[1])

#Para fitear una binomial negativa se necesita estimar el parámetro de dispersión de la binomial negativa (1/r, con r 
#la cantidad de veces que tiene que fallar la binomial negativa).
y_bines <- estimateDisp(y_bines, design)

#Para fittear los bines y sacar la señal del gen necesitamos poner como offset los coeficientes de los genes.
#Además necesitamos ajustar por el tamaño de librería
offset <- apply(fit_genes$coefficients[y_bines$genes$symbol, ],1, mean)

y_bines$offset <- matrix(rep(offset, times=ncol(y_bines$counts)), ncol=ncol(y_bines$counts), byrow=FALSE) + 
  matrix(log(rep(y_bines$samples$lib.size*y_bines$samples$norm.factors, each=nrow(y_bines$counts))), ncol=ncol(y_bines$counts), byrow=F) 

#Se ajusta con el modelo de binomial negativa las cuentas de los bines y se los compara con la referencia.
#En realidad el ajuste devuelve los beta que mejor ajustan y como es un glm los beta son logs. La resta del beta contra
#el beta de la referencia entonces es un ratio entre betas y eso es lo que termina siendo la comparación.
#Además, modificamos la función glmFit para que devuelva el error estandar. Necesitamos el error estandar para poder
#reducir los coeficientes que vienen de comparaciones con muy pocas cuentas y que por eso se disparan como si estuvieran
#cambiando mucho pero en realidad no cambian sino que están en cero y es el gen el que cambia.
fit_bines              <- glmFit(y_bines, design)
rownames(fit_bines$se) <- rownames(fit_bines$coefficients)
colnames(fit_bines$se) <- colnames(fit_bines$coefficients)

#Shrinkeamos los coeficientes usando ashr. Esto ajusta los bines con pocas cuentas para que no se disparen por
#movimientos en el gen. Como no podemos ajustar los 90000 que tiene, reducimos el número suponiendo que no
#tiene sentido pedir un qvalue mayor a 0.9 en al menos un tiempo
#Pide que el bin tenga un qvalue < qvalue_limite en cantidad_de_tiempos_limite tiempos para validad el fold change
qvalue_limite              <- 0.9
cantidad_de_tiempos_limite <- 1
bines_con_uso_diferencial  <- rownames(qvalues_bines)[apply(qvalues_bines, 1, function(x){
  sum(x<qvalue_limite) >= cantidad_de_tiempos_limite
})]

perfiles_bines <- matrix(0, nrow=length(bines_con_uso_diferencial), ncol=ncol(fit_bines$coefficients))
rownames(perfiles_bines) <- bines_con_uso_diferencial
fit_bines$se[bines_con_uso_diferencial, ][which(is.infinite(fit_bines$se[bines_con_uso_diferencial, ]), arr.ind = TRUE)] <- 1e8
for(i in 1:12){
  a <- ash(fit_bines$coefficients[bines_con_uso_diferencial, i], fit_bines$se[bines_con_uso_diferencial, i])  
  perfiles_bines[, i] <- a$result$PosteriorMean
}
colnames(perfiles_bines) <- colnames(fit_bines$coefficients[bines_con_uso_diferencial, ])

#Se queda con la data de los que hayan pasado el filtro enorme de 0.9
lfchange_bines           <- lfchange_bines[bines_con_uso_diferencial, ]
qvalues_bines            <- qvalues_bines[bines_con_uso_diferencial, ]

#Guarda las cuentas de los bines y los cambios (log fold change) entre condiciones
save(lfchange_bines, qvalues_bines, perfiles_bines, file="4_bines_prefiltrados.Rdata")
