#!/bin/sh
#g++ qcone.cpp -o qcone
#make lrs64 (the 32 bits version is still called this way)
#g++ isquad.cpp -o isquad

rm -f temp.ine
rm -f temp.ext
rm -f temp.vtx0
rm -f temp.vtx
rm -f temp.bd
rm -f temp.gd

while IFS="," read -r f0 f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26 f27 f28 f29 f30 f31
do

/home/nike/qcone 5 1 ${f0} ${f1} ${f2} ${f3} ${f4} ${f5} ${f6} ${f7} ${f8} ${f9} ${f10} ${f11} ${f12} ${f13} ${f14} ${f15} ${f16} ${f17} ${f18} ${f19} ${f20} ${f21} ${f22} ${f23} ${f24} ${f25} ${f26} ${f27} ${f28} ${f29} ${f30} ${f31} --lrs > temp.ine
/home/nike/lrs < temp.ine > temp.ext
awk '{ if ($1==1) print$0 } END { print "$" }' < temp.ext > temp.vtx0
{ printf "5 1 ${f0} ${f1} ${f2} ${f3} ${f4} ${f5} ${f6} ${f7} ${f8} ${f9} ${f10} ${f11} ${f12} ${f13} ${f14} ${f15} ${f16} ${f17} ${f18} ${f19} ${f20} ${f21} ${f22} ${f23} ${f24} ${f25} ${f26} ${f27} ${f28} ${f29} ${f30} ${f31}\n"; cat temp.vtx0; } >temp.vtx
/home/nike/isquad temp < temp.vtx
NUMLINES=$(< temp.gd wc -l)
if [ $NUMLINES -eq 2 ]
then
echo "${f0} ${f1} ${f2} ${f3} ${f4} ${f5} ${f6} ${f7} ${f8} ${f9} ${f10} ${f11} ${f12} ${f13} ${f14} ${f15} ${f16} ${f17} ${f18} ${f19} ${f20} ${f21} ${f22} ${f23} ${f24} ${f25} ${f26} ${f27} ${f28} ${f29} ${f30} ${f31}" >> "noquad.txt"
tail -n 7 temp.ext >> "noquad.txt"

else
cat temp.gd >> "hasquad.txt"
tail -n 7 temp.ext >> "hasquad.txt"

fi

rm -f temp.ine
rm -f temp.ext
rm -f temp.vtx0
rm -f temp.vtx
rm -f temp.bd
rm -f temp.gd

done < input_coeff.txt
