# script: input gtf ouput gene regulatory domain
library(GenomicRanges)
library(rtracklayer)
make_regulatory_domain <- function(gtf, basal_up = 5000, basal_down = 1000, domain_max = 1e6, protein_coding = F) {
        gtf_genes = make_genes_gr(gtf, protein_coding = protein_coding)
        gtf_genes_basal = promoters(gtf_genes, upstream = basal_up, downstream = basal_down)
        gtf_genes_basal = restrict(gtf_genes_basal, start = 0)
        # calculate the distance to extend basal regions for, to nearest gene's basal domain, ignoring strand
        pre = precede(gtf_genes, gtf_genes_basal, ignore.strand = T)
        fol = follow(gtf_genes, gtf_genes_basal, ignore.strand = T)
        pre_dist = numeric(length(pre))
        pre_dist[is.na(pre)] = domain_max
        pre_dist[!is.na(pre)] = distance(gtf_genes[!is.na(pre)], gtf_genes_basal[pre[!is.na(pre)]], ignore.strand = T)
        pre_dist = pmin(pre_dist, domain_max)
        fol_dist = numeric(length(fol))
        fol_dist[is.na(fol)] = domain_max
        fol_dist[!is.na(fol)] = distance(gtf_genes[!is.na(fol)], gtf_genes_basal[fol[!is.na(fol)]], ignore.strand = T)
        fol_dist = pmin(fol_dist, domain_max)
        # extending to both sides
        gtf_genes_domain = resize(gtf_genes, width(gtf_genes) + pre_dist, fix = "start", ignore.strand = T)
        gtf_genes_domain = resize(gtf_genes_domain, width(gtf_genes_domain) + fol_dist, fix = "end", ignore.strand = T)
        gtf_genes_domain = restrict(gtf_genes_domain, start = 0)
        # each gene will have at least basal domain
        gtf_genes_domain = punion(gtf_genes_domain, gtf_genes_basal)
        gtf_genes_domain$gene_name = gtf_genes_basal$gene_name
        gtf_genes_domain
}
gene_gr = make_genes_gr(import("./metadata/ENCODEv4_Reference/ENCFF159KBI.gtf.gz"))
saveRDS(gene_gr, file = "./metadata_processed/gene_annotations/ENCODE_Ref_gencode_v29_genes.rds")
gene_reg = make_regulatory_domain(import("./metadata/ENCODEv4_Reference/ENCFF159KBI.gtf.gz"))
saveRDS(gene_reg, file = "./metadata_processed/gene_annotations/ENCODE_Ref_gencode_v29_gene_domains.rds")

gene_gr = make_genes_gr(import("./metadata/ENCODEv4_Reference/gencode.vM21.primary_assembly.annotation_UCSC_names.gtf.gz"))
saveRDS(gene_gr, file = "./metadata_processed/gene_annotations/ENCODE_Ref_gencode_vM21_genes.rds")
gene_reg = make_regulatory_domain(import("./metadata/ENCODEv4_Reference/gencode.vM21.primary_assembly.annotation_UCSC_names.gtf.gz"))
saveRDS(gene_reg, file = "./metadata_processed/gene_annotations/ENCODE_Ref_gencode_vM21_gene_domains.rds")

# protein coding only
gene_gr = make_genes_gr(import("./metadata/ENCODEv4_Reference/ENCFF159KBI.gtf.gz"), protein_coding = T)
saveRDS(gene_gr, file = "./metadata_processed/gene_annotations/ENCODE_Ref_gencode_v29_protein_coding_genes.rds")
gene_reg = make_regulatory_domain(import("./metadata/ENCODEv4_Reference/ENCFF159KBI.gtf.gz"), protein_coding = T)
saveRDS(gene_reg, file = "./metadata_processed/gene_annotations/ENCODE_Ref_gencode_v29_protein_coding_gene_domains.rds")

gene_gr = make_genes_gr(import("./metadata/ENCODEv4_Reference/gencode.vM21.primary_assembly.annotation_UCSC_names.gtf.gz"), protein_coding = T)
saveRDS(gene_gr, file = "./metadata_processed/gene_annotations/ENCODE_Ref_gencode_vM21_protein_coding_genes.rds")
gene_reg = make_regulatory_domain(import("./metadata/ENCODEv4_Reference/gencode.vM21.primary_assembly.annotation_UCSC_names.gtf.gz"), protein_coding = T)
saveRDS(gene_reg, file = "./metadata_processed/gene_annotations/ENCODE_Ref_gencode_vM21_protein_coding_gene_domains.rds")

