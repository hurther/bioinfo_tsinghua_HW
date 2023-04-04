# samtools and bedtools
## Task 1
```shell
samtools flagstat COAD.ACTB.bam
# 输出如下
185650 + 0 in total (QC-passed reads + QC-failed reads)
4923 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
185650 + 0 mapped (100.00% : N/A)
0 + 0 paired in sequencing  # 此行为0表示应为单端测序
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
## Task 2
* secondary alignment代表的是该条read比对到多个位置，该条记录为次优比对，在双端测序中，代表hardclip，可以是read的同一部分有不同匹配区域，也可以是一条read上的不同区域
* 由`Task 1`输出中第二行可知secondary alignment有4923条
## Task 3
* 准备
```shell
# 下载BEDOPS辅助该任务
git clone https://github.com/bedops/bedops.git
cd bedops
make
make install
cp bin/* /usr/local/bin
```
* gff to bed
```shell
awk '$3=="gene"{print $0}' hg38.ACTB.gff > hg38.ATCB.gene.gff
awk '$3=="exon"{print $0}' hg38.ACTB.gff > hg38.ATCB.exon.gff
# gff2bed是BEDOPS中的命令
gff2bed <hg38.ATCB.gene.gff > hg38.ATCB.gene.bed
gff2bed <hg38.ATCB.exon.gff > hg38.ATCB.exon.bed
```
* get intron
```shell
bedtools subtract -a hg38.ATCB.gene.bed -b hg38.ATCB.exon.bed > hg38.ATCB.intron.bed
```
* get bam
```shell
bedtools intersect -abam COAD.ACTB.bam -b hg38.ATCB.intron.bed > result.bam
bedtools intersect -a hg38.ATCB.intron.bed -b COAD.ACTB.bam -c > result.txt
cat result.txt
# 输出如下
chr7    5528185 5528280 ENSG00000075624.17      .       -       HAVANA  gene    .       ID=ENSG00000075624.17;gene_id=ENSG00000075624.17;gene_type=protein_coding;gene_name=ACTB;level=2;hgnc_id=HGNC:132;havana_gene=OTTHUMG00000023268.12    10886
chr7    5529982 5530523 ENSG00000075624.17      .       -       HAVANA  gene    .       ID=ENSG00000075624.17;gene_id=ENSG00000075624.17;gene_type=protein_coding;gene_name=ACTB;level=2;hgnc_id=HGNC:132;havana_gene=OTTHUMG00000023268.12    5101
chr7    5530627 5540675 ENSG00000075624.17      .       -       HAVANA  gene    .       ID=ENSG00000075624.17;gene_id=ENSG00000075624.17;gene_type=protein_coding;gene_name=ACTB;level=2;hgnc_id=HGNC:132;havana_gene=OTTHUMG00000023268.12    268
chr7    5540771 5561851 ENSG00000075624.17      .       -       HAVANA  gene    .       ID=ENSG00000075624.17;gene_id=ENSG00000075624.17;gene_type=protein_coding;gene_name=ACTB;level=2;hgnc_id=HGNC:132;havana_gene=OTTHUMG00000023268.12    290
chr7    5561949 5562389 ENSG00000075624.17      .       -       HAVANA  gene    .       ID=ENSG00000075624.17;gene_id=ENSG00000075624.17;gene_type=protein_coding;gene_name=ACTB;level=2;hgnc_id=HGNC:132;havana_gene=OTTHUMG00000023268.12    177
chr7    5562828 5563713 ENSG00000075624.17      .       -       HAVANA  gene    .       ID=ENSG00000075624.17;gene_id=ENSG00000075624.17;gene_type=protein_coding;gene_name=ACTB;level=2;hgnc_id=HGNC:132;havana_gene=OTTHUMG00000023268.12    173
```
* bam to fastq
```shell
samtools sort -n result.sam > result.sort.sam
bedtools bamtofastq -i result.sort.sam -fq result.sort.fastq
```
## Task 4
```shell
samtools sort COAD.ACTB.bam > COAD.ACTB.sort.bam
samtools index COAD.ACTB.sort.bam
bedtools genomecov -ibam COAD.ACTB.sort.bam -bg -split > COAD.ACTB.coverage.bedgraph
```