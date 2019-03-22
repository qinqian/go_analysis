library(clusterProfiler)
library(stringr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
#hs = scan('../../hg38.fixed.3000.symbols', what='')
#eg = bitr(hs, fromType='SYMBOL', toType='ENTREZID', OrgDb='org.Hs.eg.db')
hs = scan('../../hg38.fixed.1000.symbols', what='')
eg = bitr(hs, fromType='SYMBOL', toType='ENTREZID', OrgDb='org.Hs.eg.db')
universe.df = read.csv('~/lisa_explore/background_genes/analysis_test_gene_number/AR.symbol.H3K27ac.lisa_predicted_rp.csv')
universe = bitr(str_split(universe.df[, 1], ":", simplify=TRUE)[,5], fromType='SYMBOL', toType='ENTREZID', OrgDb='org.Hs.eg.db')

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
        toType = c("ENSEMBL", "SYMBOL"),
        OrgDb = org.Hs.eg.db)

egocc = enrichGO(gene = unique(eg[, 2]), 
               #gene = gene,
               universe=unique(universe[,2]),
               #universe=names(geneList),
               OrgDb = org.Hs.eg.db,                 
               ont = 'CC', 
               pAdjustMethod = "BH",    
               pvalueCutoff  = 0.9, qvalueCutoff=0.9, readable=T)
egobp = enrichGO(gene = unique(eg[, 2]),
               #gene = gene,
               universe=unique(universe[,2]),
               #universe=names(geneList),
               OrgDb = org.Hs.eg.db,
               ont = 'BP',
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.9, qvalueCutoff=0.9, readable=T)
egomf = enrichGO(gene = unique(eg[, 2]),
               #gene = gene,
               universe=unique(universe[,2]),
               #universe=names(geneList),
               OrgDb = org.Hs.eg.db,
               ont = 'MF',
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.9, qvalueCutoff=0.9, readable=T)

pdf('hs_background_1000_genes.pdf', width=10, height=4)
#barplot(egomf, drop=TRUE, showCategory=20)
#barplot(egocc, drop=TRUE, showCategory=20)
#barplot(egobp, drop=TRUE, showCategory=20)
dotplot(egomf, showCategory=20)
dotplot(egocc, showCategory=20)
dotplot(egobp, showCategory=20)
dev.off()

#mm = scan('../../mm10.fixed.3000.symbols', what='')
#meg = bitr(mm, fromType='SYMBOL', toType='ENTREZID', OrgDb='org.Mm.eg.db')
mm = scan('../../mm10.fixed.1000.symbols', what='')
meg = bitr(mm, fromType='SYMBOL', toType='ENTREZID', OrgDb='org.Mm.eg.db')

universe.df = read.csv('/data/home/qqin/lisa_explore/test_mouse/foxd3.symbol.DNase.lisa_predicted_rp.csv')
universe = bitr(str_split(universe.df[, 1], ":", simplify=TRUE)[,5], fromType='SYMBOL', toType='ENTREZID', OrgDb='org.Mm.eg.db')

megocc = enrichGO(gene = unique(meg[, 2]), 
               universe=unique(universe[,2]),
               OrgDb = org.Mm.eg.db,                 
               ont = 'CC', 
               pAdjustMethod = "BH",    
               pvalueCutoff  = 0.9, qvalueCutoff=0.9, readable=T)

megobp = enrichGO(gene = unique(meg[, 2]), 
               universe=unique(universe[,2]),
               OrgDb = org.Mm.eg.db,                 
               ont = 'BP', 
               pAdjustMethod = "BH",    
               pvalueCutoff  = 0.9, qvalueCutoff=0.9, readable=T)

megomf = enrichGO(gene = unique(meg[, 2]),
               universe=unique(universe[,2]),
               OrgDb = org.Mm.eg.db,
               ont = 'MF',
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.9, qvalueCutoff=0.9, readable=T)

pdf('mm_background_1000_genes.pdf', width=10, height=4)
dotplot(megomf, showCategory=20)
dotplot(megocc, showCategory=20)
dotplot(megobp, showCategory=20)
dev.off()

