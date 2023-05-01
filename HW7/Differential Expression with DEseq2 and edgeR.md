# Differential Expression with DEseq2 and edgeR
## Task 1
Multiple test correction即多重检验校正。当我们要针对多个特征（例如基因）做假设检验时，由于每次检验都有一定几率会发生第一类错误（Type I errors/false positives），所以每进行一次假设检验，就会降低整体(所有假设检验)的结果的可信度，随着检验次数的增加，至少出现一次错误的可能性也在增加。多重检验校正可以有效地控制假阳性次数。

p value：衡量假阳性率的指标（False positive rate）
q value：衡量错误发现率的指标（False discovery rate，简称FDR）。由于Q value 需要利用公式从P value 校正计算后得到，所以Q value 通常又被称为adjusted p value。

一般情况下：我们可以认为Q value = FDR = adjusted p value，即三者是一个东西，虽然有些定义上的细微区别，但是问题也不大。

P value只为某一次检验负责；这次检验的假阳性率（即这次犯错的概率）
Q value是为所有次数的检验负责；按照这个设定标准，完成多次检验后，阳性显著的结果中错误的比例多大

## Task 2
- edgeR

使用trimmed mean of M values (TMM)进行标准化

TMM normalization首先需要选择一个参考样本，以它为基准进行校正。

默认下，参考样本的选择是通过比较每个样本的CPM (counts per million)的上四分位数与所有样本CPM的平均上四分位数之间的差值，找出差值最小的样本作为参考样本。

接着，在参考样本和非参考样本间两两计算校正因子（normalization factors）。
我们首先需要计算参考样本和非参考样本间的Fold change (M)和平均表达量 (A)

$$M=log_2(\frac{treated\ sample\ count}{control\ sample\ count})$$
$$A=\frac{log_2{(treated\ sample\ count)}+log_2{(control\ sample\ count)}}{2}$$
保留M和A值均为有限值的基因，并过滤极低表达量的基因。

接着，我们对M和A值进行双重截值，截掉M值排在前30%和后30%，A值排在前5%和后5%的基因，计算中间这部分基因M值的权重。

$$w(g)=\frac{D_{treatment}-K_{treatment}}{D_{treatment}\times K_{treatment}}-\frac{D_{control}-K_{control}}{D_{control}\times K_{control}}$$

其中g代表每个基因，D代表每个样本中总的reads数，K代表g的reads数。

接着，我们计算每个样本的标准化因子(C)。

$$TMM=\frac{\sum_{g\in G}w(g)\times M(g)}{\sum_{g\in G}w(g)}$$

$$C=2^{TMM}$$

最后，我们对基因进行校正：

$$\frac{C_{raw}(g)}{C}=\frac{K}{C}$$

- DESeq2

使用RLE (Relative Log Expression)进行标准化

首先，我们为每一个基因创建一个假参考样本(pseudo-reference sample)表达量，这代表理想的基因表达量。该值通过计算每个基因的几何平均值得到。

$$({\prod_{i=1}^nK_{gl}})^\frac{1}{n}=\sqrt[n]{x_1x_2...x_n}$$

接着，对于样本的每一个基因，我们计算一个 sample/reference 的比值。

样本内所有比率的中值即为校正因子(normalization factor), 也就是DESeq2 的size factor。

$$C_j=median\{\frac{K_{gj}}{\prod_{i=1}^nK_{gl}}\}$$

同edgeR一样，我们对基因进行校正即可。

## Task 3
脚本见`DESeq2.R`和`edgeR.R`

输出文本文件见`deseq2.mutant.light.vs.dark.txt`和`edger.mutant.light.vs.dark.txt`

## Task 4
脚本见`Venn_plot.R`

输出图片见`DEseq2.vs.edgeR.Venn.png`

## Task 5
- 关于log2FC大小比较，若以实际数值比较:

    脚本见`heatmap_plot_2.R`

    输出文件见`heatmap_plot_2.pdf`

- 若以log2FC的绝对值比较:

    脚本见`heatmap_plot.R`

    输出文件见`heatmap_plot.pdf`