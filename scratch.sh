#Get accession ids of each virulence factor
#print the first column of a file (-F"\t" '{print $1}')
#print text between two paranthesis ( -F "[()]")
#remove a given number of first/[to]lastcharacters from each line (cut -c N-)
#add > to the beginning of each line (sed -e 's/^/>/')
awk '/^>/' vfdb_pro.fas | awk '{print $1}' | awk -F "[()]" '{print $2}' | cut -c 4- | sed -e 's/^/>/' > vfdb_pro_acc_ids.txt

#prepare the database to DIAMOND compatible format
#replace the header line (which starts with >) with accession numbers, this
#is beacause, after diamond blast-x, this field is used by krona tools for 
#identification during visualisation.
awk '
/^>/{
getline <"vfdb_pro_accession_ids.txt"
}
1
' vfdb_pro.fas > vfdb_pro_diamondf.fas
#obtain a description file for the VFs
#by printing all lines of the database file begining with >
awk '/^>/' vfdb_pro_diamondf.fas > vfdb_description.txt

#BUILD A TXT FILE EQUIVALENT TO VFDB WITH SELECTED FIELDS
#this will be used later for obtaining meta-data of diamond detected virulent factors
#get the accession numbers
awk '/^>/' vfdb_pro_diamondf | awk '{print $1}' | awk -F "[()]" '{print $2}' | cut -c 4- | sed -e 's/>//'  > vf_acc.txt
#get the vfdb id 
awk '{print $1}' vfdb_description.txt | awk -F"(" '{print $1}' | sed 's/>//' > vf_id.txt
#get vfdb_specific_short_names, 
#the cut bit removes the first "("
#"$2" is the second field of the description file
#the sed bit of this one liner is to remove last charater(which is a ")")
awk '{print $2}' vfdb_description.txt | cut -c 2- | sed 's/.$//' > vf_short.txt
#get the full description of the virulence factors
#This bit ('{$1=$2=""}1') is to exclude certain fields from the file, 1 and 2 in this case (already extracted above)
awk '{$1=$2=""}1' vfdb_description.txt | awk -F"[" '{print $1}' > vf_full.txt
#get the short description of virulence factors
awk '{$1=$2=""}1' vfdb_description.txt | awk -F"[" '{print $2}' | sed 's/]//' > vf_cat.txt
#get the associated bacterial pathogens
awk '{$1=$2=""}1' vfdb_description.txt | awk -F"[" '{print $3}' | sed 's/]//' > vf_pathogen.txt
#create a vfdb equivalent using the above fields.
paste -d "|" vf_*.txt | sort -k 1 -t "|" > sorted_vfdb.txt

#trim sequences
while read sample; do sickle pe -f ${sample}_L001_R1_001.fastq -r ${sample}_L001_R2_001.fastq -t sanger -o ../trimmed/${sample}_qtrim.1.fastq -p ../trimmed/${sample}_qtrim.2.fastq -s ../trimmed/${sample}_qtrim_unpaired.fastq -q 20 -l 20 ;done < ../scripts/uniq_prefix.txt

#remove host sequences
while read p; do
        bowtie2 -x ../humangenome/human_DB -1 ../trimmed/${p}_qtrim.1.fastq -2 ../trimmed/${p}_qtrim.2.fastq --threads 20 -S ${p}_mapped_and_unmapped.sam
        samtools view -bS ${p}_mapped_and_unmapped.sam | samtools view -b -f 12 -F 256 | samtools sort -n > ${p}_mapped_unmapped_sorted.bam
        rm -f ${p}_mapped_and_unmapped.sam 
        bedtools bamtofastq -i ${p}_mapped_unmapped_sorted.bam -fq ${p}_host_removed_R1.fastq -fq2 ${p}_host_removed_R2.fastq
        rm -f ${p}_mapped_unmapped_sorted.bam ;
done<../scripts/uniq_prefix.txt

#concatenate Reverse and forward reads using fq2fa
while read p ; do fq2fa --merge  ${p}_host_removed_R1.fastq ${p}_host_removed_R2.fastq  ${p}_R12.fa; done < ../scripts/prefix.txt 

############# METAGENOME ASSEMBLY USING META-IDBA_UD ########################
while read p;
        do
        rm -rf ${p}
        mkdir ${p}
        idba_ud -r /home/HCV2/Alfred/AFIplus/mapping/${p}_R12.fa --num_threads 20 -o ${p} ;
done </home/HCV2/Alfred/AFIplus/scripts/uniq_prefix.txt;

#use DIAMOND to obtain potential virulence factors
while read p;do
        diamond blastx -d /home/HCV2/Alfred/AFIplus/vfdb/vfdb_pro_improved -q /home/HCV2/Alfred/AFIplus/assembly/${p}/contig.fa -p 24 -k 1  -a  ${p}
        diamond view -a ${p}.daa -o ${p}.m8 ;
done</home/HCV2/Alfred/AFIplus/scripts/uniq_prefix.txt;

#write desired fields of the diamond output: accession ids, %ge of matched seqs, E-value and sample. 
while read p;
	do awk '{print $2 "\t" $3 "\t"  $11}' ${p}_diamond/${p}.m8 ;
done<../scripts/prefix.txt > diamond_vfs.txt

#SEARCH FOR DETECTED VIRULENCE FACTOR DESCRIPTION AND ASSOCIATED FACTORS
#first is to get the accession numbers
awk -F"\t" '{print $1}' diamond_vfs.txt > diamond_vfs_accs.txt
#we are selecting the 3rd, 4th and 5th fields (Accession id, full name of vf, corresponding pathogen) and sorting the result by accession ids
awk -F"\t" '{print $2}' diamond_vfs.txt > diamond_vfs_accs.txt
grep -f "diamond_vfs_accs.txt" "vfdb.txt" | awk -F"," '{print $1"|"$3"|"$5}' | sort -k 1 -t "|" > sorted_diamond_vfs.txt

#leave out fields 2 and 4
#we can use the -t flag in the sort function to specify the field seperator.
join -t "|" sorted_vfdb.txt sorted_diamond_vfs.txt | awk -F"|" '{$2=$4=""; OFS="\t"; print $0}' > final_res.txt

awk -F"\t" '{print $0}' combined_megares.txt | sort -k 2 | awk -F"\t" '{print;}'

#TAXONOMIC PROFILLING
#Using Kaiju
for i in $(ls); do python /home/HCV2/Alfred/AFIplus/metaphinder/MetaPhinder.py -i /home/HCV2/Alfred/AFIplus/assembly/${i}/contig.fa -o ./${i} -d /home/HCV2/Alfred/AFIplus/metaphinder/database/ALL_140821_hr -b /usr/bin ; done
while read sample; do
        kaiju -v -x -z 20 -t /home/HCV2/Alfred/AFIplus/databases/kaijudb/nodes.dmp -f /home/HCV2/Alfred/AFIplus/databases/kaijudb/kaiju_db.fmi -i /home/HCV2/Alfred/AFIplus/mapping/${sample}_host_removed_R1.fastq -j /home/HCV2/Alfred/AFIplus/mapping/${sample}_host_removed_R2.fastq -o ${sample}_kaiju.out ;
done < /home/HCV2/Alfred/AFIplus/scripts/uniq_prefix.txt

#convert kaiju output to krona compatible formats
for s in $(ls *.out_names); do kaiju2krona -t /home/HCV2/Alfred/AFIplus/databases/kaijudb/nodes.dmp -n /home/HCV2/Alfred/AFIplus/databases/kaijudb/names.dmp -i $s -o ${s}.krona ; done

#add sample names to kaiju output
for f in $(ls *.out_names); do
        sample=$(echo $f | cut -f 1 -d 'k' | sed 's/_/ /g' | sed 's/ /_/')
        sed -i "s/$/\t$sample/" $f;
done

#combining the above, we have the following
for f in $(ls *.out_names); do
        kaiju2krona -t /home/HCV2/Alfred/AFIplus/databases/kaijudb/nodes.dmp -n /home/HCV2/Alfred/AFIplus/databases/kaijudb/names.dmp -i $s -o ${s}.krona
        sample=$(echo $f | cut -f 1 -d 'k' | sed 's/_/ /g' | sed 's/ /_/')
        sed -i "s/$/\t$sample/" $f;
done

#mining phage from the assembled contigs
#first is to create a directory for each of the samples
while read p; do mkdir $p; done < uniq_prefix.txt
#then loop over the directories and do metaphinder for contigs from each sample and store
#the results in corresponding directories
for i in $(ls) do;
python /home/HCV2/Alfred/AFIplus/metaphinder/MetaPhinder.py -i /home/HCV2/Alfred/AFIplus/assembly/${i}/contig.fa -o ./${i} -d /home/HCV2/Alfred/AFIplus/metaphinder/database/ALL_140821_hr -b /usr/bin;
done
#add sample names to each record in output as an extra column. This helps us to identify which
#corresponding sample after merging results into a single file.
for f in $(ls); do sed -i "s/$/\t$f/" $f/blast.out; done
for f in $(ls); do sed -i "s/$/\t$f/" $f/output.txt; done

#Then we can merge the files, first the output.txt files
#pick only phage classified records
awk -F"\t" '$2=="phage"{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t" $1"-"$7}' output.txt | sort -k 6 > positivePhage.txt
#to pick out some undesired records
grep -vwE "(#)" blast.out | awk -F"\t" '{print $2"\t"$3"\t"$11"\t" $1"-"$13}' | sort -k 4 > blast.txt
grep -vwE "(negative|not|processed|#contigID)" output.txt
#now merge the two files to only contain records that were classified as phage
join -t "\t" -1 6 -2 4 positivePhage.txt blast.txt
