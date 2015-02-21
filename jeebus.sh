#!/bin/bash


#for mgg
#root -l -b -q "Make.C(1,0)"
#root -l -b -q "Make.C(1,1)"

#for mggjj
#root -l -b -q "Make.C(0,1)"
#root -l -b -q "Make.C(0,2)"

#for 2D
for icat in `echo "0 1 2 3"`; do
    for imgg in `echo "0 1 4"`; do
	for imjj in `echo "0 1 4"`; do
	    echo "Running cat $icat mgg $imgg mjj $imjj"
	    ./R2GGBBBiasStudy_2D.exe $icat $imgg $imjj > a.txt
	done
    done
done
