#!/bin/bash


#for mgg
#root -l -b -q "Make.C(1,0)"
#root -l -b -q "Make.C(1,1)"

#for mggjj
#root -l -b -q "Make.C(0,1)"
#root -l -b -q "Make.C(0,2)"

#for 2D
searchMass="0 300"
withCorr="0 1"

catrange_m300="0 1"
catrange_m0="0 1 2 3"

#temporary settings
searchMass="0 300"
withCorr="0"

for iMass in `echo ${searchMass}`; do
    for iCorr in `echo ${withCorr}`; do
	if [ "${iMass}" -eq "300" ]; then
	    catrange=${catrange_m300}
	else
	    catrange=${catrange_m0}
	fi
	for icat in `echo ${catrange}`; do
	    for imgg in `echo "0 1 4"`; do
		for imjj in `echo "0 1 4"`; do
		    echo "Running cat $icat mgg $imgg mjj $imjj searchMass $iMass withCorr $iCorr"
		    if [ "$imgg" -eq "4" ] && [ "$imjj" -eq "4" ]; then
			./R2GGBBBiasStudy_2D.exe $icat $imgg $imjj $iMass $iCorr > output-${icat}-${imgg}-${imjj}-${iMass}-${iCorr}.txt
			sleep 100
			rm output*txt
		    else
			./R2GGBBBiasStudy_2D.exe $icat $imgg $imjj $iMass $iCorr > output-${icat}-${imgg}-${imjj}-${iMass}-${iCorr}.txt &
		    fi
		done
	    done
	done
    done
done
