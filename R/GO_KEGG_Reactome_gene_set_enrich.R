#Usage: Rscript GO_KEGG_Reactome_gene_set_enrich.R <File:gene.lst|gene.matrix> <Path:outdir> [Str:header]
args <- commandArgs(T)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(pathview)
library(topGO)
options(stringsAsFactors = FALSE)
options(bitmapType='cairo')

if (!is.na(args[3])){
  data <- read.table(args[1], header = T, row.names=NULL)
}else{
  data <- read.table(args[1], header = F)
}

outdir = args[2]
genes = unique(data[,1])
trans = bitr(genes, fromType="SYMBOL", toType=c("ENTREZID", "GENENAME"), OrgDb="org.Hs.eg.db")
# GO analysis
GO_analysis = function(genes = NULL, type = "BP", pval = 0.01, qval = 0.05, outdir = ".", Name = NULL){
  ego <- enrichGO(gene          = genes,
                  OrgDb         = org.Hs.eg.db,
                  keytype       = 'SYMBOL',
                  ont           = type,
                  pvalueCutoff  = pval,
                  qvalueCutoff  = qval)
  if (!is.null(ego) && nrow(ego@result) > 0){
    write.table(ego@result, file = paste(outdir, "/GO_enrich.xls",sep = ""),
                sep = "\t", quote = FALSE, row.names = T, col.names = NA)
    str_length = max(nchar(ego@result$Description))
    str_height = nrow(ego)
    if (str_height >15){
      pdf(paste(outdir, "/GO_", type, "_enrich_plot1.pdf",sep = ""), width = 3+(str_length)/10, height = str_height/3)
    }else{
      pdf(paste(outdir, "/GO_", type, "_enrich_plot1.pdf",sep = ""), width = 3+(str_length)/10, height = 5)
    }
    print(barplot(ego, showCategory = nrow(ego)))
    print(dotplot(ego, showCategory = nrow(ego)))
    dev.off()
    pdf(paste(outdir, "/GO_", type, "_enrich_plot2.pdf",sep = ""), width = 25, height = 25)
    enrichMap(ego)
    cnetplot(ego, categorySize="pvalue")#, foldChange=data$foldChange)
    plotGOgraph(ego)
    dev.off()
  }
}

GO_analysis(genes = genes, type = "BP", outdir = args[2], Name = Name)
GO_analysis(genes = genes, type = "MF", outdir = args[2], Name = Name)
GO_analysis(genes = genes, type = "CC", outdir = args[2], Name = Name)

# KEGG analysis
kk <- enrichKEGG(gene         = trans$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keytype = "ENTREZID")
if (nrow(kk@result) > 0){
  str_length = max(nchar(kk@result$Description))
  str_height = nrow(kk)
  write.table(kk@result, file = paste(outdir, "/KEGG_enrich.xls",sep = ""),
              sep = "\t", quote = FALSE, row.names = T, col.names = NA)
  if (str_height > 15){
    pdf(paste(outdir, "/KEGG_enrich_plot1.pdf",sep = ""), width = 3+(str_length)/10, height = str_height/3)
  }else{
    pdf(paste(outdir, "/KEGG_enrich_plot1.pdf",sep = ""), width = 3+(str_length)/10, height = 5)
  }
  print(barplot(kk, showCategory = nrow(kk)))
  print(dotplot(kk, showCategory = nrow(kk)))
  dev.off()
  pdf(paste(outdir, "/KEGG_enrich_plot2.pdf",sep = ""), width = 25, height = 25)
  enrichMap(kk)
  dev.off()
  for (i in 1:nrow(kk@result)){
    #browseKEGG(kk, 'hsa04392')
    pathway.id = kk@result$ID[i]
    pathview(gene.data  = genes,
             pathway.id = pathway.id,
             kegg.dir = outdir,
             gene.idtype ="SYMBOL",
             out.suffix = "enrich",
             species    = "hsa",
             kegg.native = TRUE)
    
    pathview(gene.data   = genes,
             pathway.id  = pathway.id,
             kegg.dir    = outdir,
             gene.idtype ="SYMBOL",
             out.suffix  = "enrich",
             species     = "hsa",
             kegg.native = F,
             same.layer  = F)
    system(paste("mv ", pathway.id, "* ", outdir, sep = ""))
    }
}

# Reactome analysis
re <- enrichPathway(gene         = trans$ENTREZID,
                    organism     = 'human',
                    pvalueCutoff = 0.05)
re <- setReadable(re, OrgDb = org.Hs.eg.db, keytype = "ENTREZID")
if (nrow(re@result) > 0){
  str_length = max(nchar(re@result$Description))
  str_height = nrow(re)
  write.table(re@result, file = paste(outdir, "/Reactome_enrich.xls",sep = ""),
              sep = "\t", quote = FALSE, row.names = T, col.names = NA)
  if (str_height > 15){
    pdf(paste(outdir, "/Reactome_enrich_plot1.pdf",sep = ""), width = 3+(str_length)/10, height = str_height/3)
  }else{
    pdf(paste(outdir, "/Reactome_enrich_plot1.pdf",sep = ""), width = 3+(str_length)/10, height = 5)
  }
  
  print(barplot(re, showCategory = nrow(re)))
  print(dotplot(re, showCategory = nrow(re)))
  dev.off()
  pdf(paste(outdir, "/Reactome_enrich_plot2.pdf",sep = ""), width = 25, height = 25)
  enrichMap(re)
  dev.off()
}

