# Bash Script
## bash 脚本
* read.sh
```shell
#!/bin/bash
Dirname=bash_homework
CUR_DIR=`ls $Dirname`

for i in $CUR_DIR
do
	if [ -f $Dirname/$i ];then
		echo $i >> "/home/test/linux/filename.txt"
	elif [ -d $Dirname/$i ];then
		echo $i >> "/home/test/linux/dirname.txt"
	else
		echo "Cannot recognize!"
	fi
done

exit 0
```
运行时在`bash_homework`的上级目录运行即可，我将`bash_homework/`保存在了`/home/test/linux`中，所以运行以下命令即可
```shell
cd linux
./read.sh
```
## 输出结果
* dirname.txt
```shell
a-docker
app
backup
bin
biosoft
c1-RBPanno
datatable
db
download
e-annotation
exRNA
genome
git
highcharts
home
hub29
ibme
l-lwl
map2
mljs
module
mogproject
node_modules
perl5
postar2
postar_app
postar.docker
RBP_map
rout
script
script_backup
software
tcga
test
tmp
tmp_script
var
x-rbp
```
* filename.txt
```shell
a1.txt
a.txt
b1.txt
bam_wig.sh
b.filter_random.pl
c1.txt
chrom.size
c.txt
d1.txt
dir.txt
e1.txt
f1.txt
human_geneExp.txt
if.sh
image
insitiue.txt
mouse_geneExp.txt
name.txt
number.sh
out.bw
random.sh
test3.sh
test4.sh
test.sh
test.txt
wigToBigWig
```