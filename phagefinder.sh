#create directories to accomodate output of the metaphinder
while read p; do mkdir $p; done < uniq_prefix.txt

#then loop over the directories and do metaphinder for contigs from each sample and store the results in corresponding directories
#and also add sample names to each record in output as an extra column. This helps us to identify which corresponding sample after merging results into a single file.
for i in $(ls) do;
  python /home/HCV2/Alfred/AFIplus/metaphinder/MetaPhinder.py -i /home/HCV2/Alfred/AFIplus/assembly/${i}/contig.fa -o ./${i} -d /home/HCV2/Alfred/AFIplus/metaphinder/database/ALL_140821_hr -b /usr/bin;
  sed -i "s/$/\t$f/" $f/blast.out
  sed -i "s/$/\t$f/" $f/output.txt
done

#then we can merge the files, first the output.txt files
#pick only phage classified records
awk -F"\t" '$2=="phage"{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t" $1"-"$7}' output.txt | sort -k 6 > positivePhage.txt

#to pick out some undesired records
grep -vwE "(#)" blast.out | awk -F"\t" '{print $2"\t"$3"\t"$11"\t" $1"-"$13}' | sort -k 4 > blast.txt

#now merge the two files to only contain records that were classified as phage
join -t "\t" -1 6 -2 4 positivePhage.txt blast.txt > merged_res.txt
