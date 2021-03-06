#!/bin/bash


#for mgg
#root -l -b -q "Make.C(1,0)"
#root -l -b -q "Make.C(1,1)"

#for mggjj
#root -l -b -q "Make.C(0,1)"
#root -l -b -q "Make.C(0,2)"

#for 2D
searchMass="0 270 300"
withCorr="0 1"

catrange_res="0 1"
catrange_nonres="0 1 2 3"

#temporary settings
if [ $# -ge 1 ]; then
    searchMass=$1
fi
if [ $# -ge 2 ]; then
    withCorr=$2
fi
if [ $# -ge 3 ]; then
    catrange_res=$3
    catrange_nonres=$3
fi


for iMass in `echo ${searchMass}`; do
    for iCorr in `echo ${withCorr}`; do
	if [ "${iMass}" -eq "0" ]; then
	    catrange=${catrange_nonres}
	else
	    catrange=${catrange_res}
	fi
	for icat in `echo ${catrange}`; do
	    for imgg in `echo "0 1 4"`; do
		for imjj in `echo "0 1 4"`; do
		    echo "Running cat $icat mgg $imgg mjj $imjj searchMass $iMass withCorr $iCorr"
		    ./R2GGBBBiasStudy_2D.exe $icat $imgg $imjj $iMass $iCorr > output-${icat}-${imgg}-${imjj}-${iMass}-${iCorr}.txt
		    rm output-${icat}-${imgg}-${imjj}-${iMass}-${iCorr}.txt
		    #if [ "$imgg" -eq "4" ] && [ "$imjj" -eq "4" ]; then
		    #	./R2GGBBBiasStudy_2D.exe $icat $imgg $imjj $iMass $iCorr > output-${icat}-${imgg}-${imjj}-${iMass}-${iCorr}.txt
		    #	sleep 100
		    #	rm output*txt
		    #else
		    #	./R2GGBBBiasStudy_2D.exe $icat $imgg $imjj $iMass $iCorr > output-${icat}-${imgg}-${imjj}-${iMass}-${iCorr}.txt &
		    #fi
		done
	    done
	done
    done
done
