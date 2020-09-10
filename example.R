library(variantprobs)
library(data.table)
library(magrittr)
source("calc_minfo.R")

data(tcga)

tcga<-setDT(tcga)

# get TCGA genes with relative frequency >= 0.01
top_genes <- tcga[
  MS == "Non-hypermutated",
  ][,
    n_tumor := length(unique(patient_id))
    ][,
      .(g_f = length(unique(patient_id)),
        n_tumor = n_tumor[1]),
      by = Hugo_Symbol
      ][,
        g_rf := g_f/n_tumor
        ][
          g_rf >= 0.01,
          ]$Hugo_Symbol


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
  Hugo_Symbol %in% top_genes,
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
  tcga_newprob_given_cancer,
  cancer_prob,
  binary_minfo = FALSE,
  normalize = TRUE
)

# see the top 20
nmi %>%
  sort(decreasing = TRUE) %>%
  head(20)

