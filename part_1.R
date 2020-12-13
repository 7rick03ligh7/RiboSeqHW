library( "DESeq2" )
library(ggplot2)
library('stringr')

rna_data = read.csv('C:/Files/Study/course2/NGS/ribo_seq/rna_data.tsv', sep='\t')
rna_data = DataFrame(rna_data)
rna_meta = read.csv('C:/Files/Study/course2/NGS/ribo_seq/rna_meta.tsv', sep='\t')
rna_meta = DataFrame(rna_meta)

dds = DESeqDataSetFromMatrix(countData=rna_data, 
                             colData=rna_meta,
                             design=~label,
                            tidy = TRUE)

dds = DESeq(dds)
res_rna = results(dds)
head(results(dds, tidy=TRUE))

summary(res_rna)
res_rna = res_rna[order(res_rna$padj),]
head(res_rna)


#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res_rna, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res_rna, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res_rna, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


ribo_data = read.csv('C:/Files/Study/course2/NGS/ribo_seq/ribo_data.tsv', sep='\t')
ribo_data = DataFrame(ribo_data)
ribo_meta = read.csv('C:/Files/Study/course2/NGS/ribo_seq/ribo_meta.tsv', sep='\t')
ribo_meta = DataFrame(ribo_meta)

dds = DESeqDataSetFromMatrix(countData=ribo_data, 
                             colData=ribo_meta,
                             design=~label,
                             tidy = TRUE)

dds = DESeq(dds)
res_ribo = results(dds)
head(results(dds, tidy=TRUE))

summary(res_ribo)
res_ribo = res_ribo[order(res_ribo$padj),]
head(res_ribo)

par(mfrow=c(1,1))
with(res_ribo, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res_ribo, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res_ribo, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

rna = subset(res_rna, padj<.01 & abs(log2FoldChange)>2)
ribo = subset(res_ribo, padj<.01 & abs(log2FoldChange)>2)

dim(rna)
dim(ribo)
