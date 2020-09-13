#!/bin/sh
#g++ qcone.cpp -o qcone
#cd lrslib-070-32bit 
#make lrs64  (the 32 bits version is still called this way)
#cd ..
#g++ isquad.cpp -std=c++11 -o isquad

CASENO=0
while IFS="," read -r f0 f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26 f27 f28 f29 f30 f31
do
CASENO=$(($CASENO + 1))

#change the first two numbers to n, m for n-var m-aux

/home/nike/projects/rrg-mcrowley/nike/QuadratizationCode/LRS-based/qcone 5 2 ${f0} ${f1} ${f2} ${f3} ${f4} ${f5} ${f6} ${f7} ${f8} ${f9} ${f10} ${f11} ${f12} ${f13} ${f14} ${f15} ${f16} ${f17} ${f18} ${f19} ${f20} ${f21} ${f22} ${f23} ${f24} ${f25} ${f26} ${f27} ${f28} ${f29} ${f30} ${f31} --lrs > temp_${CASENO}.ine
/home/nike/projects/rrg-mcrowley/nike/QuadratizationCode/LRS-based/lrslib-070-32bit/lrs < temp_${CASENO}.ine > temp_${CASENO}.ext
awk '{ if ($1==1) print$0 } END { print "$" }' < temp_${CASENO}.ext > temp_${CASENO}.vtx0

#change the first two numbers to n, m for n-var m-aux

{ printf "5 2 ${f0} ${f1} ${f2} ${f3} ${f4} ${f5} ${f6} ${f7} ${f8} ${f9} ${f10} ${f11} ${f12} ${f13} ${f14} ${f15} ${f16} ${f17} ${f18} ${f19} ${f20} ${f21} ${f22} ${f23} ${f24} ${f25} ${f26} ${f27} ${f28} ${f29} ${f30} ${f31}\n"; cat temp_${CASENO}.vtx0; } >temp_${CASENO}.vtx
/home/nike/projects/rrg-mcrowley/nike/QuadratizationCode/LRS-based/isquad temp_${CASENO} < temp_${CASENO}.vtx
NUMLINES=$(< temp_${CASENO}.gd wc -l)
if [ $NUMLINES -eq 2 ]
then
echo "${f0} ${f1} ${f2} ${f3} ${f4} ${f5} ${f6} ${f7} ${f8} ${f9} ${f10} ${f11} ${f12} ${f13} ${f14} ${f15} ${f16} ${f17} ${f18} ${f19} ${f20} ${f21} ${f22} ${f23} ${f24} ${f25} ${f26} ${f27} ${f28} ${f29} ${f30} ${f31}" >> "noquad.txt"
tail -n 7 temp_${CASENO}.ext >> "noquad.txt"

else
cat temp_${CASENO}.gd >> "hasquad.txt"
tail -n 7 temp_${CASENO}.ext >> "hasquad.txt"

fi

done < input_coeff.txt
