#!/bin/bash

for F in `ls *_noclip_*.txt`
do
	OUT=`echo ${F//_noclip/}`
	cat $F | awk '{if(NR%2==0) print substr($0, 7, length); else print $0;}' > $OUT	
done
