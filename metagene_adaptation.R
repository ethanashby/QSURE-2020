#Script to calculate NMI for metagene groups standardized by metagene size
#remove topgene filter
#instead of grouping by Hugo_Symbol we're grouping by the 'Cluster' column in tcga_clust, the left_join of tcga 
#and the cluster_assignment matrix

library(variantprobs)
library(data.table)
library(magrittr)
library(tidyverse)
source("calc_minfo.R")


#read in data
data(tcga)

#filter by nonhypermutated mutation signature and non-NA Cancer Code
tcga_v_f <- tcga_clust[
  # filter out all hypermutated tumors
  # & tumors with unknown cancer
  MS == "Non-hypermutated" & !is.na(Cancer_Code),
  ][
    ,
    # number of tumors per cancer type
    n_tumor := length(unique(patient_id)),
    by = Cancer_Code
    ][,
      # variant frequencies by CLUSTER, cancer type
      # also save n_tumor per cancer
      .(v_f = length(unique(patient_id)),
        n_tumor = n_tumor[1]),
      by = .(Variant, Cluster, Cancer_Code)
      ]

tcga_newprob_given_cancer <- tcga_v_f[,
    # Calculate Good Turing probabilities of
    # at least one new variants per CLUSTER & cancer type
    {
      GT_probs <- goodturing_probs(
        counts = v_f,
        m = n_tumor[1]
      )
      .(p_atleast_1new_perMB = GT_probs['atleast_1new'])
    },
    by = .(Cluster, Cancer_Code)
    ] %>%
  # cast the long data frame into wide format -- fill missing values
  # (i.e. no mutation for that (CLUSTER, cancer) pair)
  # with ZERO
  # replace fill = 0 by
  # fill = 1 - exp(- 1/(length(unique(tcga$patient_id)) + 1))
  # to exactly recreate Fig 3b in Nat. Comm. paper
  dcast(
    Cluster ~ Cancer_Code,
    value.var = "p_atleast_1new_perMB",
    fill = 0
  ) %>%
  # put  Clusters as row names
  # and convert the data table into a matrix
  magrittr::set_rownames(.$Cluster) %>%
  .[, Cluster := NULL] %>%
  data.matrix()

#cancer probs are agnostic to metagene size... correct for metagene size
tcga_newprob_given_cancer_perMB<-tcga_newprob_given_cancer *1e6 /cluster_lengths$cluster_length

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

nmi <- calc_minfo(
  tcga_newprob_given_cancer_perMB,
  cancer_prob,
  binary_minfo = FALSE,
  normalize = TRUE
)

# see the top 20
nmi %>%
  sort(decreasing = TRUE) %>%
  head(20)


