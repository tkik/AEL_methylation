oarsub -n homer_$bed -l core=16,walltime=24 -I


source activate homer

/mnt/std-pool/homedirs/rtoth/scMeth_HSC/dmr_beds /mnt/std-pool/homedirs/rtoth/scMeth_HSC/dmr_beds_promoter
bed_dirs=( /mnt/std-pool/homedirs/rtoth/scMeth_HSC/dmr_beds_enhancer)


for bed_dir in $bed_dirs
do

if [[ $bed_dir == *"enhancer"* ]]; then
  output_dir=/mnt/std-pool/homedirs/rtoth/scMeth_HSC/res_enhancer
elif [[ $bed_dir == *"promoter"* ]]; then
  output_dir=/mnt/std-pool/homedirs/rtoth/scMeth_HSC/res_promoter
else
  output_dir=/mnt/std-pool/homedirs/rtoth/scMeth_HSC/res_all
fi

#bed_dir=/mnt/std-pool/homedirs/rtoth/scMeth_HSC/dmr_beds
#output_dir=/mnt/std-pool/homedirs/rtoth/scMeth_HSC/res_all
beds=`ls -l $bed_dir | awk '{print $9}'`


#perl ../configureHomer.pl -install hg19
# might need to preparse the genome as well. Then add -preparse option.
for bed in $beds
do
comp="${bed%.*}"
echo $comp
#mkdir $output_dir/$comp
#oarsub -l /nodes=4,/walltime=12:00:00  
findMotifsGenome.pl $bed_dir/$bed hg19  $output_dir/$comp -preparse -len 8,10,12 -size 100 -S 8 -p 8 -cache 6921 -fdr 0
done

done

