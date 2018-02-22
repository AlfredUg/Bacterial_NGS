#taxonomic classification of host free-sequence reads
while read sample; do
  kaiju -v -x -z 20 -t /home/HCV2/Alfred/AFIplus/databases/kaijudb/nodes.dmp -f /home/HCV2/Alfred/AFIplus/databases/kaijudb/kaiju_db.fmi -i /home/HCV2/Alfred/AFIplus/mapping/${sample}_host_removed_R1.fastq -j /home/HCV2/Alfred/AFIplus/mapping/${sample}_host_removed_R2.fastq -o ${sample}_kaiju.out ;
done < /home/HCV2/Alfred/AFIplus/scripts/uniq_prefix.txt

#convert output to krona compatible format and append sample names to kaiju output files
for f in $(ls *.out_names); do
        kaiju2krona -t /home/HCV2/Alfred/AFIplus/databases/kaijudb/nodes.dmp -n /home/HCV2/Alfred/AFIplus/databases/kaijudb/names.dmp -i $s -o ${s}.krona
        sample=$(echo $f | cut -f 1 -d 'k' | sed 's/_/ /g' | sed 's/ /_/')
        sed -i "s/$/\t$sample/" $f;
done
