# Expression Matrix
## Task 1
- 对于small RNA-seq
$$CPM/RPM = \frac{Number\ of\ reads\ mapped\ to\ a\ gene}{Total\ number\ of\ mapped\ reads}\times 10^6$$
- 对于poly-A/total RNA-seq和paired end RNA-seq
$$RPKM = \frac{Number\ of\ reads\ mapped\ to\ a\ gene}{Total\ number\ of\ mapped\ reads\times gene\ length\ in\ Kbp}\times 10^6$$
$$FPKM=0.5\times RPKM$$
$$TPM=\frac{RPKM}{\sum (RPKM)}\times 10^6$$
- 以及一些其他的方式

例如edgeR中利用TMM counts；RLE中利用RLE等等

## Task 2
1. E
2. D
3. A

## Task 3
- 判断建库方式
```shell
infer_experiment.py -r GTF/Arabidopsis_thaliana.TAIR10.34.bed -i bam/Shape02.bam
# 输出如下
# Loading SAM/BAM file ...  Total 200000 usable reads were sampled


# This is PairEnd Data
# Fraction of reads failed to determine: 0.0315
# Fraction of reads explained by "1++,1--,2+-,2-+": 0.4769
# Fraction of reads explained by "1+-,1-+,2++,2--": 0.4916
# 由于"1++,1--,2+-,2-+"与"1+-,1-+,2++,2--"的比例几乎相同，故有很大的把握认定这个数据是由非链特异性建库（strand nonspecific）得到的
```
- 计算 read count matrix
```shell
featureCounts -s 0 -p -t exon -g gene_id \
-a GTF/Arabidopsis_thaliana.TAIR10.34.gtf \
-o result/Shape02.featurecounts.exon.txt bam/Shape02.bam
cd result
cut -f 1,7 Shape02.featurecounts.exon.txt | grep -v '^#' > Shape02.count.txt
grep AT1G09530 Shape02.count.txt
# 输出如下
# AT1G09530       86
```

## Task 4
合并数据见 [craete_matrix.R](Rscript/create_matrix.R)

计算CPM和z-scores见 [calculate.R](Rscript/calculate.R)

画图见 [plot.R](Rscript/plot.R)

输出文件在output文件夹中

分析[热图](output/Rplot.png)知，COAD和READ的转录本是最相近的

