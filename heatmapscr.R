betas <- coef(dds)
colnames(betas)


topGenes <- head(order(res$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 1 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)

library(pheatmap)


library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))
heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
dev.copy(png,'deseq2_Heatmap.png')
dev.off()

heatmap.2(mat, breaks = seq(-1.5, 1.5, 0.3),
          key.xtickfun = function(){breaks <- parent.frame()$breaks
          return(list(at = parent.frame()$scale01(seq(-1.5,1.5,0.5)),
                      labels = c("< -1.5", as.character(seq(-1,1,0.5)), "> 1.5")))},)



print(plotPCA(rld, intgroup=c('condition')))
dev.copy(pdf,'deseq2_pca.pdf')
dev.off()