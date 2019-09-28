for i in {1..51} ; do cp quadGruber.sh 5var_all_neg/job_$i/ ; done
for i in {1..51} ; do cd job_$i ; sed -i -e  '8s/job1_5varAllnegWithCubics/job'"$i"'_5varAllnegWithCubics/' 0submit ; cd .. ;  done
for i in {4..51} ; do cd job_$i ; sbatch 0submit ; cd .. ; done

