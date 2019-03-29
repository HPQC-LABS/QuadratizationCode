#!/bin/sh
g++ qcone.cpp -o qcone
make lrs
g++ isquad.cpp -o isquad

file="input_coeff_$1.txt"
while IFS="," read -r f1 f2 f3 f4 f5 f6
do
./qcone --n=5 --m=1 --a1234=$f1 --a2345=$f2 --a3451=$f3 --a4512=$f4 --a5123=$f5 --a12345=$f6 --lrs > temp.ine
./lrs < temp.ine > temp.ext
awk '{ if ($1==1) print$0 } END { print "$" }' < temp.ext > temp.vtx0
{ printf "5 1 ${f1} ${f2} ${f3} ${f4} ${f5} ${f6}\n"; cat temp.vtx0; } >temp.vtx
./isquad temp < temp.vtx
NUMLINES=$(< temp.gd wc -l)
if [ $NUMLINES -eq 2 ]
then
echo "${f1} ${f2} ${f3} ${f4} ${f5} ${f6}" >> "noquad_$1.txt"
else
cat temp.gd >> "hasquad_$1.txt"
fi
rm temp.ine
rm temp.ext
rm temp.vtx0
rm temp.vtx
rm temp.bd
rm temp.gd
done <"$file"
