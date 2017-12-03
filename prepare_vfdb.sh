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

#BUILD A CSV FILE EQUIVALENT TO VFDB WITH SELECTED FIELDS
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
#go to R and create a csv file equivalent using these fields. (see customVFDB.R)

#SEARCH FOR DETECTED VIRULENCE FACTOR DESCRIPTION AND ASSOCIATED PATHOGENS
#we are selecting the 4th and 5th fields (full name of vf and corresponding pathogen)
grep -f "selected.txt" "vfdb.csv" | awk -F"," '{print $4","$5}' > detected_vfs.txt
#write desired fields of the diamond output ()accession ids, %ge of matched seqs and E-value. 
while read p; do awk '{print $2 "\t" $3 "\t"  $11}' ${p}_diamond/${p}.m8; done<../scripts/prefix.txt 