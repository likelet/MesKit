# 1.The function and steps of this project: cluster analysis
# Main steps:
# - Reads data txt input.txt
# - Do the KEGG analysis
# - Save the outcome in output.pdf

# 2.See whether these pachages exist on comp. If not, install.

# Clean the working environment
# rm(list=ls())

# Load pachages
#library("BiocManager")
#library("BiocInstaller")
library("topGO")
#library("KEGG.db")
library("Rgraphviz")
library("clusterProfiler")
library("org.Hs.eg.db")
library("ggplot2")
options(warn = -1)

Args <- commandArgs()
cat("Args[6]=",Args[6],"\n")
cat("Args[7]=",Args[7],"\n")
maf<-read.table(Args[6],sep = "\t",header = TRUE)
KEGG1.pdf<-Args[7]

KEGGanalysis<-function(f){
  
  #Get raw data, put them in seperate lists
  genenames<-list()
  for (n in (1:length(levels(f$Tumor_Sample_Barcode)))){
    genenames[[n]]=as.vector(f$Hugo_Symbol[which (f$Tumor_Sample_Barcode == levels(f$Tumor_Sample_Barcode)[n])])
  }
  
  #Divide the genes if there is a "," and name the seperate vectors.
  new_genenames<-list()
  a=1
  while (a<length(genenames) || a==length(genenames)){
    for (b in 1:(length(genenames[[a]]))){
      new_genenames[[a]]=c(strsplit(genenames[[a]][b],split = ","),new_genenames[a])}
    new_genenames[[a]]<-unique(unlist(new_genenames[[a]]))
    a=a+1
  }
  names(new_genenames)<-c(levels(maf$Tumor_Sample_Barcode))
  
  
  #Do the KEGG analysis seperately 
  u=1
  while (u>0 && u<length(new_genenames)+1){
    currentdata_u<- bitr(new_genenames[[u]], fromType="SYMBOL", 
                         toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
    KEGG_u<- enrichKEGG(currentdata_u[,2], organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
                        qvalueCutoff = 0.1)
    
    #pick up the suitable results
    if (all(KEGG_u@result$p.adjust>0.05)){
      print("NO cluster.")
      u=u+1
    }else{
      needresult_u<-KEGG_u@result[c(2,3,6,9)]
      realresult_u<-needresult_u[which(needresult_u$p.adjust<0.05),]
      for (i in 1:length(realresult_u$GeneRatio)){
        realresult_u$GeneRatio[i]<-round(as.numeric(strsplit(realresult_u$GeneRatio[i],split = "/")[[1]][1])/
                                         as.numeric(strsplit(realresult_u$GeneRatio[i],split = "/")[[1]][2]),2)
      }
      
    #draw pictures
      plotdraw<-function(finalresult){
        ggplot(data=finalresult)+
        aes(GeneRatio,Description)+
        geom_point()+
        geom_point(aes(size=Count,color=p.adjust))+
        scale_color_gradient(low="pink",high = "blue")+
        labs(size="Count",x="GeneRatio",y="Pathway name",title="Pathway enrichment")
      }
      pdf(paste("KEGG_analysis.",names(new_genenames)[u],".pdf",sep = ""),width=7,height=8)
      print(plotdraw(realresult_u))
      dev.off()
      u=u+1
    }
  }
}
KEGGanalysis(maf)