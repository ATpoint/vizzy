library(vizzy)
library(patchwork)
library(ggrepel)

set.seed(1)
dds <- DESeq2::DESeq(DESeq2::makeExampleDESeqDataSet(5000,10))
res <- DESeq2::results(dds) %>% data.frame %>% na.omit
res$baseMean<-log2(res$baseMean+1)
theme_set(theme_classic(base_size = 12.5))

dir.create("plots")

png("plots/MAplots.png")
(
ggMAplot(xval = res$baseMean, yval = res$log2FoldChange, pval = res$pvalue,
         title = "MA-plot", subtitle = "most basic type of MA-plot") | 

ggMAplot(xval = res$baseMean, yval = res$log2FoldChange, pval = res$pvalue,
         y.ablines = c(log2(1.5),-log2(1.5)), title = "MA-plot", subtitle = "with ablines") ) /

(ggMAplot(xval = res$baseMean, yval = res$log2FoldChange, pval = res$pvalue,
         y.ablines = c(log2(1.5),-log2(1.5)), title = "MA-plot", subtitle = "with altered ylims",
         ylim=c(-2,2)) | 

ggMAplot(xval = res$baseMean, yval = res$log2FoldChange, pval = res$pvalue,
         title = "MA-plot", subtitle = "individual genes highlighted",
         labels=ifelse(rownames(res) %in% c("gene1153", "gene1828"), rownames(res), "")) +
         geom_text_repel(aes(label=labels), show.legend=FALSE, max.overlaps=Inf, min.segment.length=0) )
dev.off()

png("plots/Volcanos.png")
(
  ggMAplot(xval = res$log2FoldChange, yval = -log10(res$pvalue), pval = res$pvalue,
           title = "Volcano", subtitle = "most basic type of Volcano",
           preset="volcano", ylab="-log10(pvalue)") | 
    
    ggMAplot(xval = res$log2FoldChange, yval = -log10(res$pvalue), pval = res$pvalue,
             title = "Volcano", subtitle = "with ablines",
             x.ablines = c(-log2(1.5), log2(1.5)),
             preset="volcano", ylab="-log10(pvalue)") ) /
  
  (ggMAplot(xval = res$log2FoldChange, yval = -log10(res$pvalue), pval = res$pvalue,
            title = "Volcano", subtitle = "with altered xlims",
            x.ablines = c(-log2(1.5), log2(1.5)), xlim = c(-2,2),
            preset="volcano", ylab="-log10(pvalue)") | 
     
     ggMAplot(xval = res$log2FoldChange, yval = -log10(res$pvalue), pval = res$pvalue,
              title = "Volcano", subtitle = "individual genes highlighted",
              preset="volcano", ylab="-log10(pvalue)",
              labels=ifelse(rownames(res) %in% c("gene1153", "gene1828"), rownames(res), "")) +
     geom_text_repel(aes(label=labels), show.legend=FALSE, max.overlaps=Inf, min.segment.length=0) )
dev.off()



