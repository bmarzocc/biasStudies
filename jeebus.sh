#!/bin/bash


#for mgg
#root -l -b -q "Make.C(1,0)"
#root -l -b -q "Make.C(1,1)"

#for mggjj
#root -l -b -q "Make.C(0,1)"
#root -l -b -q "Make.C(0,2)"

#for 2D
searchMass=0
#searchMass=300

withCorr="1 0"

catrange="0 1"
if [ "$searchMass" -eq "0" ]; then
    catrange="0 1 2 3"
fi

for iCorr in `echo ${withCorr}`; do
    for icat in `echo ${catrange}`; do
	for imgg in `echo "0 1 4"`; do
	    for imjj in `echo "0 1 4"`; do
		echo "Running cat $icat mgg $imgg mjj $imjj searchMass $searchMass withCorr $iCorr"
		./R2GGBBBiasStudy_2D.exe $icat $imgg $imjj $searchMass $iCorr > a.txt
	    done
	done
    done
done
