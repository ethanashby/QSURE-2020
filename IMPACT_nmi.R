#IMPACT NMI

data(impact)
setDT(impact)

#list of impact genes
impact_genelist <- impact[MS=="Non-hypermutated" & !is.na(Cancer_Type)] %>% select(Hugo_Symbol) %>% unique() %>% unlist() %>% unname()

####calculate exome sizes for these genes

exome_sizes<-data.table::fread("gencode.v19.basic.exome.bed")
###replace X chromosome with 23
exome_sizes$V1<-as.numeric(exome_sizes$V1)
exome_sizes$V1[is.na(exome_sizes$V1)]<-23

exome_sizes$V2<-as.numeric(exome_sizes$V2)
exome_sizes$V3<-as.numeric(exome_sizes$V3)
exome_sizes<-exome_sizes[complete.cases(exome_sizes),]
colnames(exome_sizes)<-c("chrom", "chromStart", "chromEnd")

genelengths<-data.frame(chrom=epigenetic_data$chr[epigenetic_data$gene %in% impact_genelist], chromStart=epigenetic_data$start[epigenetic_data$gene %in% impact_genelist], chromEnd=epigenetic_data$end[epigenetic_data$gene %in% impact_genelist])

setDT(exome_sizes)
setDT(genelengths)
setkey(genelengths)
#if any exon overlaps with gene coordinates, then join them
annotated_exons<-foverlaps(exome_sizes, genelengths, type="any", nomatch=0L)
annotated_exons

genelengths<-data.frame(chrom=epigenetic_data$chr[epigenetic_data$gene %in% impact_genelist], chromStart=epigenetic_data$start[epigenetic_data$gene %in% impact_genelist], chromEnd=epigenetic_data$end[epigenetic_data$gene %in% impact_genelist], gene=epigenetic_data$gene[epigenetic_data$gene %in% impact_genelist])

#full join exon coordinates and gene names
exon_annotations<-full_join(annotated_exons, genelengths, by=c("chromStart")) %>% select(i.chromStart, i.chromEnd, gene)

#genes listed by their exon sizes
exon_sizes<- exon_annotations %>% mutate(length=i.chromEnd-i.chromStart) %>% group_by(gene) %>% summarize(size=sum(length)) %>% as.data.frame()

#genes in impact that we don't have metafeature annotations for... only return those genes with annotations
impact_genelist<-setdiff(impact_genelist, setdiff(impact_genelist, exon_sizes$gene))


#######Calculate variant probs

# TCGA variant frequencies & n_tumors by cancer sites
tcga_v_f <- tcga[
  # filter out all hypermutated tumors
  # & tumors with unknown cancer
  MS == "Non-hypermutated" & !is.na(Cancer_Code),
  ][
    ,
    # number of tumors per cancer type
    n_tumor := length(unique(patient_id)),
    by = Cancer_Code
    ][,
      # variant frequencies by gene, cancer type
      # also save n_tumor per cancer
      .(v_f = length(unique(patient_id)),
        n_tumor = n_tumor[1]),
      by = .(Variant, Hugo_Symbol, Cancer_Code)
      ]


# matrix of P(at least 1 new variant in Gene - g | Cancer k)
# for all genes with rel freq >= 0.01
tcga_newprob_given_cancer <- tcga_v_f[
  # keep only top_genes
  Hugo_Symbol %in% impact_genelist,
  ][,
    # Calculate Good Turing probabilities of
    # at least one new variants per gene & cancer type
    {
      GT_probs <- goodturing_probs(
        counts = v_f,
        m = n_tumor[1]
      )
      .(p_atleast_1new_v = GT_probs['atleast_1new'])
    },
    by = .(Hugo_Symbol, Cancer_Code)
    ] %>%
  # cast the long data frame into wide format -- fill missing values
  # (i.e. no mutation for that (gene, cancer) pair)
  # with ZERO
  # replace fill = 0 by
  # fill = 1 - exp(- 1/(length(unique(tcga$patient_id)) + 1))
  # to exactly recreate Fig 3b in Nat. Comm. paper
  dcast(
    Hugo_Symbol ~ Cancer_Code,
    value.var = "p_atleast_1new_v",
    fill = 1 - exp(- 1/(length(unique(tcga$patient_id)) + 1))
  ) %>%
  # put  Hugo_Symbols as row names
  # and convert the data table into a matrix
  magrittr::set_rownames(.$Hugo_Symbol) %>%
  .[, Hugo_Symbol := NULL] %>%
  data.matrix()


###### Adjust probability by exon size and scale to probability per 10-kb

tcga_newprob_given_cancerper10kb=(tcga_newprob_given_cancer * 1e4/exon_sizes$size[match(rownames(tcga_newprob_given_cancer), exon_sizes$gene)])


# Find the relative frequencies of cancer sites
cancer_npatient <- tcga[
  !is.na(Cancer_Code) & MS == "Non-hypermutated",
  .(cancer_npatient = length(unique(patient_id))),
  by = Cancer_Code
  ][,
    cancer_rf := cancer_npatient/sum(cancer_npatient)
    ]


cancer_prob <- cancer_npatient$cancer_rf %>%
  setNames(cancer_npatient$Cancer_Code)

# calc_minfo has two main arguments:
# (a) prob_mat: matrix (( p_gc )), where
# p_ij = P(Event associated with gene - g | Cancer - c)
# and (b) cancer_prob: vector ( p_c ) where p_c = P(Cancer - c).
# names(cancer_prob) and colnames(prob_mat) should be identical
# The function will return a vector
# of length nrow(prob_mat), for all genes in prob_mat
nmi_IMPACT <- calc_minfo(
  tcga_newprob_given_cancerper10kb,
  cancer_prob,
  binary_minfo = FALSE,
  normalize = TRUE
)

nmi_IMPACT_melt<-nmi_IMPACT %>% melt()
nmi_IMPACT_melt %>% max()