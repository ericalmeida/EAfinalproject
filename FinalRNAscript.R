source('http://bioconductor.org/biocLite.R')
biocLite('DESeq2')
library('DESeq2')

library('DESeq2')
directory<-'~/Desktop/BIO594_FinalProject/RNAFinal'
#use grep to search for the 'treated' part of filename to collect files
HTS_Files<-grep('treated',x,value=TRUE)
HTS_Files<-grep('treated',list.files(directory),value=TRUE)
# sampleFiles
#[1] 'treated1.txt'   'treated2.txt' 'treated3.txt'  'untreated1.txt'
#[5] 'untreated2.txt' 'untreated3.txt'

SACondition<-c('treated','treated','untreated','untreated')
SATable<-data.frame(sampleName=HTS_Files, fileName=HTS_Files, condition=SACondition)
####
#sampleTable
#     sampleName       fileName condition
#1   treated1.txt   treated1.txt   treated
#2   treated2.txt   treated2.txt   treated
#3   treated3.txt   treated3.txt   treated
#4 untreated1.txt untreated1.txt untreated
#5 untreated2.txt untreated2.txt untreated
#6 untreated3.txt untreated3.txt untreated
######

ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=SATable, directory=directory, design=~condition)
#####
#ddsHTSeq
#class: DESeqDataSet
#dim: 7921 6
#exptData(0):
#assays(1): counts
#rownames(7921): seq_1 seq_2 ... seq_7920 seq_7921
#rowData metadata column names(0):
#colnames(6): treated1.txt treated2.txt ... untreated2.txt
#  untreated3.txt
#colData names(1): condition
#######
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('untreated','treated'))

dds<-DESeq(ddsHTSeq)
res<-results(dds)
res<-res[order(res$padj),]
head(res)

plotMA(dds,ylim=c(-4,4),main='Number of S. aureus Gene Counts')
dev.copy(png,'deseq2_MAplot.png')
dev.off()

mcols(res,use.names=TRUE)

rld<- rlogTransformation(dds, blind=TRUE)
vsd<-varianceStabilizingTransformation(dds, blind=TRUE)

library('RColorBrewer')
library('gplots')
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale='none',
          dendrogram='none', trace='none', margin=c(10,6))
dev.copy(png,'DESeq2_heatmap1')
dev.off()
heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale='none',
          dendrogram='none', trace='none', margin=c(10, 6))
dev.copy(png,'DESeq2_heatmap2')
dev.off()
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale='none',
          dendrogram='none', trace='none', margin=c(10, 6))
dev.copy(png,'DESeq2_heatmap3')
dev.off()

distsRL <- dist(t(assay(rld)))
mat<- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition,HTS_Files, sep=' : '))

#updated in latest vignette (See comment by Michael Love)
#this line was incorrect
#heatmap.2(mat, trace='none', col = rev(hmcol), margin=c(16, 16))
#From the Apr 2015 vignette
hc <- hclust(distsRL)
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace='none',
          col = rev(hmcol), margin=c(13, 13))
dev.copy(png,'deseq2_heatmaps_samplebysample.png')
dev.off()

hc <- hclust(distsRL)
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace='none',
          col = rev(hmcol), margin=c(13, 13))
dev.copy(png,'deseq2_heatmaps_samplebysample.png')
dev.off()

print(plotPCA(rld, intgroup=c('condition')))
dev.copy(png,'deseq2_pca.png')
dev.off()

ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < .1,
             cleaned = results(ddsClean)$padj < .1)
addmargins(tab)
write.csv(as.data.frame(tab),file='sim_condition_treated_results_cleaned_summary_deseq2.csv')
resClean <- results(ddsClean)
write.csv(as.data.frame(resClean),file='sim_condition_treated_results_cleaned_deseq2.csv')

#filtering threashold
attr(res,'filterThreshold')
#     10%
#91.48005
plot(attr(res,'filterNumRej'),type='b', ylab='number of rejections')
dev.copy(png,'deseq2_filtering_treshold.png')
dev.off()

W <- res$stat
maxCooks <- apply(assays(dds)[['cooks']],1,max)
idx <- !is.na(W)
plot(rank(W[idx]), maxCooks[idx], xlab='rank of Wald statistic',
     ylab='maximum Cook's distance per gene',
     ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
     m <- ncol(dds)
     p <- 3
     abline(h=qf(.99, p, m - p))
     dev.copy(png,'deseq2_cooksdist.png')
     dev.off()



use <- res$baseMean > attr(res,'filterThreshold')
table(use)
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c('do not pass'='khaki', 'pass'='powderblue')
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = '', ylab='frequency')
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend('topright', fill=rev(colori), legend=rev(names(colori)))
