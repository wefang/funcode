This step integrates human and mouse experiments for chromatin accessibiltiy (CA) and histone ChIP-seq data

For chromatin accessibility:
1. 'remove_chromatin_batch_effect.R' removes batch effect between DNase-seq and ATAC-seq data for human and mouse experiments separately

2. 'integrate_chromatin.R' integrates all CA experiments between two species

For Histone ChIP-seq:

1. 'human/mouse_reverse_transfer.R' transfers histone ChIP-seq data to the CA data space

2. 'integrate_human_mouse_histone.R' integrates all CA and transferred histone ChIP-seq data between human and mouse

