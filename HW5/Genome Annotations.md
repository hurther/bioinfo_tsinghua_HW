# Genome Annotations
## Task 1
人类基因组是对24条染色体测序的结果。

根据ENSEMBL的数据显示，人类基因组（GRCh38.p13）有3,096,649,726个碱基对，contig长度为 3.4 GB，染色体长度为 3.1 GB

具体分类为：
* **Coding genes**	19,827 (excl 649 readthrough)
* **Non coding genes**	25,967
  * **Small non coding genes**	4,864
  * **Long non coding genes**	18,882
  * **Misc non coding genes**	2,221
* **Pseudogenes**	15,241

根据NCBI的数据显示，人类基因组（GRCh38.p14）有59,444个基因，基因组大小为 3.1 GB

具体分类为：
* **Protein-coding gene**	20,636
* **Non-coding gene**	22,283
* **Small RNAs**	2,841
* **Pseudogenes**	17,392
* **Other**	21,293

## Task 2
如果根据ENSEMBL的数据计算，非编码RNA的比例大致在$\frac{25967}{19827+25967+15241}\times100\%=42.5\%$

如果根据NCBI的数据计算，非编码RNA的比例大致在$\frac{22283+2841}{59444}\times100\%=42.3\%$

实际上由于人类目前对ncRNA的研究正处于起步阶段，未来这个比例或许还会进一步上升

从已有的分类来看，ncRNA可以分为以下类别：

* snRNA 参与mRNA的剪接
* tRNA  参与mRNA的翻译
* rRNA  最丰富的RNA分子，形成核糖体的框架
* snoRNA 对RNA进行化学修饰
* miRNA 与mRNA结合抑制转录
* siRNA 类似于miRNA，防止翻译
* lncRNA 参与转录的启动子特异性调控、表观遗传调控等