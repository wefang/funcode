This step runs on a high performance computing cluster (HPC) with SGE scheduler

1. 'download_files.R' and 'download_control.R' downloads bam files based on ENCODE metadata tables and saves the fragment coordinates as '.rds' files.

2. 'process.R' counts number of fragments in the input regions, and saves the output as integer binary files

3. 'matrix.R' collects the counts from all samples and compile into a region by experiment numeric matrix

4. 'normalize.R' normalizes raw counts by total library size

