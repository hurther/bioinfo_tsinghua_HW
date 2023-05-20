# ChIP-seq
## Task 1
input 是 ChIP 实验中的阳性对照，样品经过超声破碎后取出一部分作为 input，input 不进行 ChIP 实验，因而包含样本超声后释放的所有 DNA、蛋白质。

input 的作用包括：
1. 首先，input 可以帮助检测超声破碎的效果，以及目的蛋白在样品中的表达情况。

2. 其次，input 可以帮助判断 ChIP 结果，通过 input 和 IP（抗体免疫沉淀获得的目的蛋白）对比，判断 ChIP 获得的是否为目的蛋白，以及 ChIP 的效果。若 IP 条带较 input 明显或差不多，说明 ChIP 富集效果较高，若 IP 条带较浅，说明 ChIP 富集效果一般。

3. 最后， ChIP 测序数据分析时更离不开 input，input 是寻找结合峰（也就是 Peak Calling）的定盘星。在 ChIP 测序分析中，我们将 input 和 IP 的 reads 平均覆盖深度做归一化处理，采取 input 的 reads 作为背景，使用 MACS 软件对比 input 和 IP 的 reads 的归一化平均覆盖深度（normalized average depth），进行 Peak Calling。

## Task 2
### findPeaks
findPeaks有7种模式（对应参数中的style），分别为：
1. factor（在ChIP-seq中寻找转录因子）
2. histone（在ChIP-seq中寻找组蛋白区域）
3. super（找到super enhancer）
4. groseq（从strand specific GRO-Seq识别转录本）
5. tss（从5'RNA-Seq/CAGE或5'GRO-Seq data识别promotor/TSS）
6. dnase（针对DNase-seq peak finding调节过参数）
7. mC（DNA甲基化分析）

使用方法为：
```bash
findPeaks <tag directory> -style <...> -o auto -i <control tag directory>
```

### findMotifsGenome
使用方法为：
```bash
findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
```
Basic options:
1. -bg \<background position file\> (genomic positions to be used as background, default=automatic) 去除和目标位置重合的背景位置
2. -len <#>[,<#>,<#>...] (motif length, default=8,10,12)
3. -size <#> (fragment size to use for motif finding, default=200)
4. -S <#> (Number of motifs to optimize, default: 25)
5. -mis <#> (global optimization: searches for strings with # mismatches, default: 2)

更多选项设置见[官方文档](http://homer.ucsd.edu/homer/ngs/peakMotifs.html)

## Task 3
以下步骤均在 `/home/test/chip-seq/` 下进行:
```bash
cd /home/test/chip-seq/
```
准备输出目录
```bash
mkdir output
```
### Peak Calling
```bash
makeTagDirectory homework/ip    homework/ip.chrom_part.bam
makeTagDirectory homework/input homework/input.chrom_part.bam
findPeaks homework/ip/ -style factor -o homework_output/part.peak -i homework/input/
```
### Motif Finding
```bash
findMotifsGenome.pl homework_output/part.peak sacCer2 homework_output/part.motif.output -len 8
```
### 数据处理
```bash
grep -v ^# part.peak | awk '{if($11>=8 && $12<1e-8) print}' > peak.result
```
输出结果见`output`文件夹
