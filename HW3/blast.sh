#!/bin/bash
var1=1
var2=2
while [ "$var1" -le "9" ]
do
while [ "$var2" -le "10" ]
do
blastp -query seq$var1.txt -subject seq$var2.txt -out output/seq$var1-$var2.blastp
((var2++))
done
((var1++))
var2=$(($var1+1))
done

