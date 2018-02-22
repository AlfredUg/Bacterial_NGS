#trim sequences, remove host sequences, concatenate Reverse and forward reads using fq2fa 
while read p; do
        sickle pe -f ${sample}_L001_R1_001.fastq -r ${sample}_L001_R2_001.fastq -t sanger -o ../trimmed/${sample}_qtrim.1.fastq -p ../trimmed/${sample}_qtrim.2.fastq -s ../trimmed/${sample}_qtrim_unpaired.fastq -q 20 -l 20 
        bowtie2 -x ../humangenome/human_DB -1 ../trimmed/${p}_qtrim.1.fastq -2 ../trimmed/${p}_qtrim.2.fastq --threads 20 -S ${p}_mapped_and_unmapped.sam
        samtools view -bS ${p}_mapped_and_unmapped.sam | samtools view -b -f 12 -F 256 | samtools sort -n > ${p}_mapped_unmapped_sorted.bam
        rm -f ${p}_mapped_and_unmapped.sam 
        bedtools bamtofastq -i ${p}_mapped_unmapped_sorted.bam -fq ${p}_host_removed_R1.fastq -fq2 ${p}_host_removed_R2.fastq
        rm -f ${p}_mapped_unmapped_sorted.bam
        fq2fa --merge  ${p}_host_removed_R1.fastq ${p}_host_removed_R2.fastq  ${p}_R12.fa
        
        #metagenome assembly using meta-idba_ud
        rm -rf ${p}
        mkdir ${p}
        idba_ud -r /home/HCV2/Alfred/AFIplus/mapping/${p}_R12.fa --num_threads 20 -o ${p}
        
        #use diamond to obtain potential virulence factors
        diamond blastx -d /home/HCV2/Alfred/AFIplus/vfdb/vfdb_pro_improved -q /home/HCV2/Alfred/AFIplus/assembly/${p}/contig.fa -p 24 -k 1  -a  ${p}
        diamond view -a ${p}.daa -o ${p}.m8       
done<../scripts/uniq_prefix.txt

#write desired fields of the diamond output: accession ids, %ge of matched seqs, E-value and sample. 
while read p;
	do awk '{print $2 "\t" $3 "\t"  $11}' ${p}_diamond/${p}.m8 ;
done<../scripts/prefix.txt > diamond_vfs.txt

#SEARCH FOR DETECTED VIRULENCE FACTOR DESCRIPTION AND ASSOCIATED FACTORS
#first is to get the accession numbers
awk -F"\t" '{print $2}' diamond_vfs.txt > diamond_vfs_accs.txt
#we are selecting the 3rd, 4th and 5th fields (Accession id, full name of vf, corresponding pathogen) and sorting the result by accession ids
grep -f "diamond_vfs_accs.txt" "vfdb.txt" | awk -F"," '{print $1"|"$3"|"$5}' | sort -k 1 -t "|" > sorted_diamond_vfs.txt

#leave out fields 2 and 4
#we can use the -t flag in the sort function to specify the field seperator.
join -t "|" sorted_vfdb.txt sorted_diamond_vfs.txt | awk -F"|" '{$2=$4=""; OFS="\t"; print $0}' > final_res.txt
