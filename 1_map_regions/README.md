This step maps peaks from one species to another using sequence alignment. For example, to map human DHS peaks mouse:

1. obtain the human DHS peak file (ENCFF03GCK.tsv) from ENCODE portal
2. obtain the corresponding pairwise alignment file in axt.net file format (eg. hg38.mm10.axt.net) from UCSC genome browser, link https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsMm10/
3. obtain the chromatin size file of mouse mm10 from UCSC genome browser(https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/)
4. split the axt.net file into multiple files according to chromosome
5. modify the parameters and run the base_run.py to align submit of human DHS peak to the mm10 genome
6. extend the human DHS peak and aligned nucleotide base in mouse genome to 200 base pair length region
7. use the output file from base_run.py to run the two_region_hg38tomm10.py to calculate the percentage of idenditcal base pairs between human regions and correspond mouse regions