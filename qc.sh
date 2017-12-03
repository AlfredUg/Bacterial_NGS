#manipulating file names
for i in $(ls *.fastq); do echo "$i" | awk -F"_" '{print $1,$2}' | awk '{gsub(" ", "_", $0); print}';  done

#concatenate Reverse and forward reads using fq2fa
while read p ; do fq2fa --merge  ${p}_L001_R1_001.fastq ${p}_L001_R2_001.fastq  ${p}_R12.fa; done < ../prefix.txt 
