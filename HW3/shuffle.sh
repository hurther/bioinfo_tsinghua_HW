#!/bin/bash
seq="MSTRSVSSSSYRRMFGGPGTASRPSSSRSYVTTSTRTYSLGSALRPSTSRSLYASSPGGVYATRSSAVRL"
var=1
bg=0
L=${#seq}
ed=$(($L-1))
shuffled=""
while [ "$var" -le "10" ]
do
for i in `seq $bg $ed | shuf `;do
    shuffled=$shuffled${seq:$i:1}
done
echo $shuffled > seq$var.txt
shuffled=""
((var++))
done

