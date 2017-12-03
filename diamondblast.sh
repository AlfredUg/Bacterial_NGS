#This script  compares contigs to bacterial virulence factor database #

#WORKING_DIR="/home/HCV2/Alfred/AFIplus" #working directory
#SCRIPTS_DIR=$WORKING_DIR/scripts #directory containing scripts
#DB_DIR=$WORKING_DIR/vfDB #the directory of the database to be blasted against.
#THIS_DIR=$WORKING_DIR/AFI_test #we suppose that this script is run in the directory containing the assembly created folders(where contigs are contained)
while read p;
	do
	diamond blastx -d /home/HCV2/Alfred/AFIplus/vfdb/vfdb_pro -q /home/HCV2/Alfred/AFIplus/AFI_test/${p}/contig.fa -p 24 -k 1  -a  ${p}
	diamond view -a ${p}.daa -o ${p}.m8 ;
done</home/HCV2/Alfred/AFIplus/scripts/prefix.txt;