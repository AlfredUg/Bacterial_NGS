#this script combines reverse and forward reads (using fq2fa) ahead of assembly 
#(using idba-metagenome assembler).
while read p;
	do
		#concatenate Reverse and forward reads using fq2fa
		fq2fa --merge  ${p}_L001_R1_001.fastq ${p}_L001_R2_001.fastq  ${p}_R12.fa
		rm -rf ${p}
		mkdir ${p}
		idba_ud -r ${p}_R12.fa -o ${p} ;
done </home/HCV2/Alfred/AFIplus/scripts/prefix.txt;