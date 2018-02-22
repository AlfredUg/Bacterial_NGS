
#write desired fields of the diamond output: accession ids, %ge of matched seqs, E-value and sample. 
while read p; do 
	#use diamond to obtain potential virulence factors
        diamond blastx -d /home/HCV2/Alfred/AFIplus/vfdb/vfdb_pro_improved -q /home/HCV2/Alfred/AFIplus/assembly/${p}/contig.fa -p 24 -k 1  -a  ${p}
        diamond view -a ${p}.daa -o ${p}.m8       
	awk '{print $2 "\t" $3 "\t"  $11}' ${p}_diamond/${p}.m8 ;
done<../scripts/prefix.txt > diamond_vfs.txt

#SEARCH FOR DETECTED VIRULENCE FACTOR DESCRIPTION AND ASSOCIATED FACTORS
#first is to get the accession numbers
awk -F"\t" '{print $2}' diamond_vfs.txt > diamond_vfs_accs.txt
#we are selecting the 3rd, 4th and 5th fields (Accession id, full name of vf, corresponding pathogen) and sorting the result by accession ids
grep -f "diamond_vfs_accs.txt" "vfdb.txt" | awk -F"," '{print $1"|"$3"|"$5}' | sort -k 1 -t "|" > sorted_diamond_vfs.txt

#leave out fields 2 and 4
#we can use the -t flag in the sort function to specify the field seperator.
join -t "|" sorted_vfdb.txt sorted_diamond_vfs.txt | awk -F"|" '{$2=$4=""; OFS="\t"; print $0}' > final_res.txt
