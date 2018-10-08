library(GenomicRanges)
datos_bines            <- data.frame(features@bins[names(ccor[ccor < -0.8]), ])
rownames(datos_bines)  <- names(ccor[ccor < -0.6])

datos_bines            <- data.frame(features@bins[c("AT3G54500:I005", "AT4G25500:I004"), ])
rownames(datos_bines)  <- c("AT4G25500:E005", "AT4G25500:I004")

secuencias <- apply(datos_bines,1, function(a){
  s <- subseq(fasta[a["seqnames"]], start=as.numeric(a["start"]), end=as.numeric(a["end"]))
  if(a["strand"] == "-") s <- reverseComplement(s)
  return(s <- as.character(s))
})

unique(c(grep("TAG", secuencias), grep("TGA", secuencias), grep("TAA", secuencias)))


#perfiles_bines       <- perfiles_bines

ccor <- sapply(rownames(perfiles_bines), function(x){
  gen <- strsplit2(x, ":")[, 1]
  if(gen %in% rownames(perfiles_genes)){
    return(cor(perfiles_bines[x, ], perfiles_genes[gen, ]))
  }else{
    return(0)
  }    
})
ccor[abs(ccor) > 0.8]

isoforma3<-"gaaaaaaatatctttgaaactaaaaagaaaagaaaataaaataaaaatccaaatgatctctcttgccttctttatctaaatcgtacctaactctgttttcagtttcacctatctctcgacgcgagcttcttcttcttctagggcttccacgcgactacgcctgccaaaatcattctacaggaagcatgaagccagtcttctgtgggaactttgagtatgatgcgcgcgaaggtgacctggaacgactattcaggaaatacggcaaggttgagagggttgatatgaaagctggatgtgtttgataatcttgggacgcctatactcacacctggccatgaatggaccatagggaaatccatctcttatccctcaaaatcagtcactttccatccgcattagcactgccatcttaattgcatttcattcctcatcactttgcacacttggacatgcctctgcaatggagaaggctctctattttcattcatatccccgtctgacttccacattttcagttgttccttgacatattaattccataatgcaagggtttgcttttgtatacatggaagatgaaagggatgcggaagatgccatccgagcacttgaccgctttgaatttgggcgtaagggacgcagacttcgtgttgaatggacaaagagtgaacgtggaggtgataaaagatctggtggtggttcaaggagatcctcatccagcatgagaccttccaagactctctttgtgattaactttgatgcggataatactaggacccgggatctagagaaacactttgagccgtatggaaagatcgtaaacgttaggatcaggaggaattttgcatttatccagtacgaggcacaagaggatgccaccagagcattggatgcttcaaataacagtaagctgatggataaggtgatctcggtggagtatgctgtgaaggatgatgatgctagaggaaatggacacagtcctgaaagacgccgtgataggtcacctgaaaggagaaggcgatcacctagtccttacaaaagagaaagaggaagccctgattatggccgaggagctagtcctgttgctgcctacagaaaggaaaggaccagtcctgactatggtcgaagacgtagcccaagtccttacaagaaatcaagacgtggcagtcccgagtatggtcgtgaccgcagaggcaatgatagccctcgcaggagggagagagtcgcaagccctactaagtacagccgcagtcccaacaacaagagagagaggatgagccctaatcacagcccgttcaagaaggagagtccgagaaatggggttggtgaagttgaaagtcccattgaaaggagagagagatcgaggtctagccccgagaatggccaagttgaaagccctgggtcaataggaagaagagacagtgatggtgggtatgatggtgcagagagcccaatgcagaagagccggtctcctcgttcgccaccagctgacgagtgataagagtggatccacaatctctatcaaagtaggatgttgtaactgtttgtagtcaacaacgctatgtcgtcgtgattagtttttgtcgtttggtttttgatagattcgaactcggatcacttttattgtcggattgaaaaacttatgtcaagttactacatttccgtttttttttatcacctctcagtttgaacccc"
isoforma3<-toupper(isoforma3)
grep(secuencias[1], isoforma3)
gregexpr(pattern = secuencias[1], isoforma3) 

gregexpr(pattern = "TGA", secuencias[1])
gregexpr(pattern = "TAG", secuencias[1])
gregexpr(pattern = "TAA", secuencias[1]) 
for(s in secuencias){
  for(i in 1:3){
    print(suppressWarnings(translate(DNAString(substr(s, i, nchar(s))))))
  }
}

s = subseq(fasta$`3`, 20177705, 20177934)
s <- reverseComplement(s)
for(i in 1:3){
  print(suppressWarnings(translate(DNAString(substr(s, i, nchar(s))))))
}

a            <- data.frame(features@genes[c("AT4G25500"), ])

secuencias <- apply(datos_bines,1, function(a){
  s <- subseq(fasta[a["seqnames"]], start=as.numeric(a["start"]), end=as.numeric(a["end"]))
  if(a["strand"] == "-") s <- reverseComplement(s)
  return(s <- as.character(s))
})

layout(matrix(1:4, ncol=2, nrow=2))
matplot(t(cuentas_bines_promedios[c("AT4G25500:E005", "AT4G25500:I004"), 13:24]), type="l", main="12 cruda")
#matplot(t(cuentas_bines_promedios[c("AT4G25500:E005", "AT4G25500:I004"), 1:12]), type="l", main="22")
matplot(t(rbind(fit_genes$fitted.values["AT4G25500", ], cuentas_genes_promedios[c("AT4G25500"), 13:24])), type="l", main="12 cruda gen")
#matplot(t(cuentas_genes_promedios[c("AT4G25500"), 1:12]), type="l", main="22")
matplot(t(perfiles_bines[c("AT4G25500:E005", "AT4G25500:I004"), ]), type="l", main="12 fiteada ")
matplot(perfiles_genes["AT4G25500", ], type="l", main="12 fiteada gen")

cor(t(rbind(perfiles_bines[c("AT4G25500:E005", "AT4G25500:I004"), ], perfiles_genes["AT4G25500", ])))
#--------------------------------------------
layout(matrix(1:4, ncol=2, nrow=2))
matplot(t(cuentas_bines_promedios[c("AT3G54500:I007"), 13:24]), type="l", main="12 cruda")
#matplot(t(cuentas_bines_promedios[c("AT4G25500:E005", "AT4G25500:I004"), 1:12]), type="l", main="22")
matplot(t(rbind(cuentas_genes_promedios[c("AT3G54500"), 13:24])), type="l", main="12 cruda gen")
#matplot(t(cuentas_genes_promedios[c("AT4G25500"), 1:12]), type="l", main="22")
matplot((perfiles_bines[c("AT3G54500:I007"), ]), type="l", main="12 fiteada ")
matplot(perfiles_genes["AT3G54500", ], type="l", main="12 fiteada gen")
#--------------------------------------------
layout(matrix(1:4, ncol=2, nrow=2))
matplot(t(cuentas_bines_promedios[c("AT2G21660:I002"), 13:24]), type="l", main="12 cruda")
#matplot(t(cuentas_bines_promedios[c("AT4G25500:E005", "AT4G25500:I004"), 1:12]), type="l", main="22")
matplot(t(rbind(cuentas_genes_promedios[c("AT2G21660"), 13:24])), type="l", main="12 cruda gen")
#matplot(t(cuentas_genes_promedios[c("AT4G25500"), 1:12]), type="l", main="22")
matplot((perfiles_bines[c("AT2G21660:I002"), ]), type="l", main="12 fiteada ")
matplot(perfiles_genes["AT2G21660", ], type="l", main="12 fiteada gen")
#--------------------------------------------
layout(matrix(1:4, ncol=2, nrow=2))
matplot(t(cuentas_bines_promedios[c("AT4G39260:I005"), 13:24]), type="l", main="12 cruda")
#matplot(t(cuentas_bines_promedios[c("AT4G25500:E005", "AT4G25500:I004"), 1:12]), type="l", main="22")
matplot(t(rbind(cuentas_genes_promedios[c("AT4G39260"), 13:24])), type="l", main="12 cruda gen")
#matplot(t(cuentas_genes_promedios[c("AT4G25500"), 1:12]), type="l", main="22")
matplot((perfiles_bines[c("AT4G39260:I005"), ]), type="l", main="12 fiteada ")
matplot(perfiles_genes["AT4G39260", ], type="l", main="12 fiteada gen")
#--------------------------------------------
layout(matrix(1:4, ncol=2, nrow=2))
matplot(t(cuentas_bines_promedios[c("AT5G11260:I005"), 13:24]), type="l", main="12 cruda")
#matplot(t(cuentas_bines_promedios[c("AT4G25500:E005", "AT4G25500:I004"), 1:12]), type="l", main="22")
matplot(t(rbind(cuentas_genes_promedios[c("AT5G11260"), 13:24])), type="l", main="12 cruda gen")
#matplot(t(cuentas_genes_promedios[c("AT4G25500"), 1:12]), type="l", main="22")
matplot((perfiles_bines[c("AT5G11260:I005"), ]), type="l", main="12 fiteada ")
matplot(perfiles_genes["AT5G11260", ], type="l", main="12 fiteada gen")
#--------------------------------------------
layout(matrix(1:4, ncol=2, nrow=2))
matplot(t(cuentas_bines_promedios[c("AT3G02830:I005"), 13:24]), type="l", main="12 cruda")
#matplot(t(cuentas_bines_promedios[c("AT4G25500:E005", "AT4G25500:I004"), 1:12]), type="l", main="22")
matplot(t(rbind(cuentas_genes_promedios[c("AT3G02830"), 13:24])), type="l", main="12 cruda gen")
#matplot(t(cuentas_genes_promedios[c("AT4G25500"), 1:12]), type="l", main="22")
matplot((perfiles_bines[c("AT3G02830:I005"), ]), type="l", main="12 fiteada ")
matplot(perfiles_genes["AT3G02830", ], type="l", main="12 fiteada gen")
#--------------------------------------------
layout(matrix(1:4, ncol=2, nrow=2))
matplot(t(cuentas_bines_promedios[c("AT4G31800:I005"), 13:24]), type="l", main="12 cruda")
#matplot(t(cuentas_bines_promedios[c("AT4G25500:E005", "AT4G25500:I004"), 1:12]), type="l", main="22")
matplot(t(rbind(cuentas_genes_promedios[c("AT4G31800"), 13:24])), type="l", main="12 cruda gen")
#matplot(t(cuentas_genes_promedios[c("AT4G25500"), 1:12]), type="l", main="22")
matplot((perfiles_bines[c("AT4G31800:I005"), ]), type="l", main="12 fiteada ")
matplot(perfiles_genes["AT4G31800", ], type="l", main="12 fiteada gen")
#--------------------------------------------
layout(matrix(1:4, ncol=2, nrow=2))
matplot(t(cuentas_bines_promedios[c("AT3G56290:I001", "AT3G56290:I002"), 13:24]), type="l", main="12 cruda")
#matplot(t(cuentas_bines_promedios[c("AT4G25500:E005", "AT4G25500:I004"), 1:12]), type="l", main="22")
matplot(t(rbind(cuentas_genes_promedios[c("AT3G56290"), 13:24])), type="l", main="12 cruda gen")
#matplot(t(cuentas_genes_promedios[c("AT4G25500"), 1:12]), type="l", main="22")
matplot(t(perfiles_bines[c("AT3G56290:I001", "AT3G56290:I002"), ]), type="l", main="12 fiteada ")
matplot(perfiles_genes["AT3G56290", ], type="l", main="12 fiteada gen")
cor(t(rbind(perfiles_bines[c("AT3G56290:I001", "AT3G56290:I002"), ], perfiles_genes["AT3G56290", ])))

, , , AT5G11260