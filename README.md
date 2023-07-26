### Functional Conservation of DNA Elements (FUNCODE)

Scripts for compute FUNCODE scores from ENCODE DNase-seq, ATAC-seq and Histoen ChIP-seq data.

To reproduce the FUNCODE scores, use the following steps:

1. Generate pairs of mapped DHS elements based on pairwsie sequence alignment (UCSC Blastz)

2. Download and process sequence read alignment files (bam) from ENCODE data portal

3. Integrate data between human and mouse

4. Compute FUNCODE CO-V and CO-B scores for sequence alignable DHS pairs

5. Compute FUNCODE CO-V scores for sequence DHS pairs with nearby gene homology, and calls significantly conserved pairs

For details on each step, see README in the corresponding sub-directories.

Report bugs and provide suggestions by sending email to:

Author and maintainer: Weixiang Fang (wfang9@jh.edu)
