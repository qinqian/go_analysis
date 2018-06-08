source("https://bioconductor.org/biocLite.R")
#biocLite("gage")
#biocLite("pathview")
biocLite("org.Mm.eg.db")
library(KEGG.db)
#install.packages("gagaData")

library("org.Mm.eg.db")
library(gage)
library(pathview)
library(gageData)

data("kegg.gs")
data("kegg.sets.mm")

data(go.sets.mm)
data(go.subs.mm)

for (i in list.files('.', pattern='.csv$')) {
  f = read.csv(i)
  f = f[f$p_val_adj<=0.05 & abs(f$avg_logFC) > 0.41,]
  print(dim(f))
  meanfc = f[,3]
  names(meanfc) = f[,1]
  names(meanfc) = id2eg(ids=names(meanfc), category = "SYMBOL", org="Mm")[,2] 
  result = gage(meanfc, gsets=kegg.sets.mm, ref=NULL, samp=NULL)
  write.table(rbind(result$greater, result$less),
              sep='\t', quote=F, file=paste0(i, '_kegg_mm.xls'))
  ##write.table(result$stats, sep='\t', quote=F, file=paste0(i, '_stats.xls'))
  
  result.bp = gage(meanfc, gsets=go.sets.mm[go.subs.mm$BP], ref=NULL, samp=NULL)
  write.table(rbind(result.bp$greater, result.bp$less),
              sep='\t', quote=F, file=paste0(i, '_BP.xls'))
  
  result.mf = gage(meanfc, gsets=go.sets.mm[go.subs.mm$MF], ref=NULL, samp=NULL)
  write.table(rbind(result.mf$greater, result.mf$less),
              sep='\t', quote=F, file=paste0(i, '_MF.xls'))  
  
  result.cc = gage(meanfc, gsets=go.sets.mm[go.subs.mm$CC], ref=NULL, samp=NULL)
  write.table(rbind(result.cc$greater, result.cc$less),
              sep='\t', quote=F, file=paste0(i, '_CC.xls'))    
}
