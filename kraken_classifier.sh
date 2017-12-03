while read p;
	do 
		mkdir -p /path/to/${p}_kraken
		kraken --preload --threads 24 --fastq-input --paired --db /path/to/krakendb ${p}_L001_R1.fq ${p}_L001_R2.fq > /path/to/${p}_kraken ;
done < /path/to/prefix.txt ;
