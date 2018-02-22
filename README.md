
### assembly.sh
+ Metagenome assembly using [UDBA-UD](https://github.com/loneknightpy/idba), generates contigs which are blasted against certain databases (e.g the virulence factor database ([VFDB](http://www.mgc.ac.cn/VFs/)).
### prepare_vfdb.sh
+ Manipulates VFDB to [DIAMOND](https://github.com/bbuchfink/diamond) compatible format.
### krakenDB_Creation.py
+ Creates a [kraken](http://ccb.jhu.edu/software/kraken/) compatible bacterial database fetched from [NCBI](https://www.ncbi.nlm.nih.gov/genome/microbes/) 
### kraken_classifier.sh
+ Read based taxonomic assignment, uses database created using `krakenDB_Creation.py`. 

