# Linux Basic Command Homework
[TOC]
## 准备
* 打开 docker 容器
```docker
docker exec -it zhonghaochen_linux bash
```
* 提取文件（已下载文件至挂载文件夹中）
```shell
cd share
mv test_command.gtf ../linux/
```
* 进入工作目录
```shell
cd ../linux
```
## 任务1
* 行数统计
```shell
wc -l test_command.gtf | awk '{print $1}'
# 输出如下
8
```
* 字数统计
```shell
wc -c test_command.gtf | awk '{print $1}'
# 输出如下
636
```
## 任务2
```shell
grep 'chr_' test_command.gtf | grep 'YDL248W'
# 输出如下
chr_IV  ensembl gene    1802    2953    .       +       .       gene_id "YDL248W"; gene_version "1";
chr_IV  ensembl transcript      802     2953    .       +       .       gene_id "YDL248W"; gene_version "1";
chr_IV  ensembl start_codon     1802    1804    .       +       0       gene_id "YDL248W"; gene_version "1";
```
## 任务3
```shell
sed 's/chr_/chromosome_/g' test_command.gtf | cut -f 1,3,4,5
# 输出如下
chromosome_IV   gene    1802    2953
chromosome_IV   transcript      802     2953
chromosome_IV   exon    1802    2953
chromosome_IV   CDS     1802    950
chromosome_IV   start_codon     1802    1804
chromosome_IV   stop_codon      2951    2953
chromosome_IV   gene    762     3836
chromosome_IV   transcript      3762    836
```
## 任务4
```shell
awk '{t=$2;$2=$3;$3=t;print;}' test_command.gtf | sort -k 4 -k 5 -n > result.gtf
# 输出如下
cat result.gtf
chromosome_IV gene ensembl 762 3836 . + . gene_id "YDL247W-A"; gene_version "1";
chr_IV transcript ensembl 802 2953 . + . gene_id "YDL248W"; gene_version "1";
chromosome_IV CDS ensembl 1802 950 . + 0 gene_id "YDL248W"; gene_version "1";
chr_IV start_codon ensembl 1802 1804 . + 0 gene_id "YDL248W"; gene_version "1";
chr_IV gene ensembl 1802 2953 . + . gene_id "YDL248W"; gene_version "1";
chromosome_IV exon ensembl 1802 2953 . + . gene_id "YDL248W"; gene_version "1";
chromosome_IV stop_codon ensembl 2951 2953 . + 0 gene_id "YDL248W"; gene_version "1";
chr_IV transcript ensembl 3762 836 . + . gene_id "YDL247W-A"; gene_version "1";
```
## 任务5
* 变更前
```shell
ls -lh test_command.gtf
# 输出如下
-rwxrwxrwx 1 test test 636 Feb 25 14:22 test_command.gtf
```
* 变更后
```shell
chmod 774 test_command.gtf | ls -lh test_command.gtf
# 输出如下
-rwxrwxr-- 1 test test 636 Feb 25 14:22 test_command.gtf
```
## 其他
* 删除相关文件并退出 docker
```shell
rm test_command.gtf result.gtf
exit
```