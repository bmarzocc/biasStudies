#!/usr/bin/env python

import sys

searchStr="resMass=0,cat0,withCorr=0"
if len(sys.argv)>1:
    searchStr=sys.argv[1]

file=open('resultsBias2D.txt','r')
printStr="Truth Model & Exp1,Exp1 & Pow1,Pow1 & Ber1,Ber1 & Ber1,Ber2 & Ber1,Ber3 & Ber1,Ber4 & Ber2,Ber1 \\\\ \hline \n"
expexpBiasList=[]
highBiasFlag=False

for line in file:
    if line.find(searchStr)<0:
        continue
    data=line.replace(searchStr,'').replace('Mgg','').replace('Mjj','').replace('N','').replace('\n','').replace('\x00','').split()

    expexpBiasList.append(float(data[1]))

    for i in data[1:]:
        if float(i)>1.0 or highBiasFlag:
            highBiasFlag=True
            break

    for i in data[:-1]:
        printStr+=i+" & "
    printStr+=data[-1]+" \\\\ \n"

aveBias = sum(expexpBiasList)/float(len(expexpBiasList))

print "\nThe table:"
print printStr
print "\nExp1,Exp1 has an average bias of %.1f\\%%.\n" % (100*sum(expexpBiasList)/float(len(expexpBiasList)))
if highBiasFlag:
    print "Note that some entries show a high bias\n"

    
