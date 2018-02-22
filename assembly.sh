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
        
done<../scripts/uniq_prefix.txt
