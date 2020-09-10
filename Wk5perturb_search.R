#######
#Build metagene clusters and calc NMIs 
#######

library(dplyr)
library(magrittr)
library(knitr)
library(devtools)
library(magrittr)
library(tidyverse)
library(data.table)
#options(buildtools.check = function(action) TRUE )
#install_github("https://github.com/c7rishi/variantprobs.git")
library(variantprobs)
library(rvest)
library(reshape)
library(boot)
library(knitr)
library(kableExtra)
library(parallel)
library(e1071)
library(drc)
library(spatstat.utils)
library(DescTools)
library(GGally)
library(parallel)
library(lattice)
library(phyclust)
source("calc_minfo.R")

data(tcga)
tcga_genelist<-tcga[MS == "Non-hypermutated" & !is.na(Cancer_Code),]$Hugo_Symbol %>% unique()

#########
#Build covariate table of metafeatures

#read in covariates and remove empty rows
epigenetic_data<-read.csv("covariate-table.csv", header=TRUE)
epigenetic_data<-epigenetic_data[epigenetic_data$gene!="",]
epigenetic_data<-epigenetic_data[epigenetic_data$gene %in% tcga_genelist,]

#create space of epigenetic predictors, Z scale predictors, only complete cases
cluster_space<-cbind(epigenetic_data$expression_CCLE, epigenetic_data$replication_time, epigenetic_data$noncoding_mutation_rate, epigenetic_data$local_GC_content, epigenetic_data$HiC_compartment)
rownames(cluster_space)<-epigenetic_data$gene
cluster_space<-cluster_space[complete.cases(cluster_space),]
cluster_space<-scale(cluster_space)
cluster_space<-as.data.frame(cluster_space)
colnames(cluster_space)<-c("expression", "rep_time", "nc_mut_rate", "loc_GC_cont", "HiC")
#only include clusters in tcga
cluster_space<-cluster_space[rownames(cluster_space) %in% tcga_genelist, ]

###########
#Calculate gene sizes

exome_sizes<-data.table::fread("gencode.v19.basic.exome.bed")
###replace X 
exome_sizes$V1<-as.numeric(exome_sizes$V1)
exome_sizes$V1[is.na(exome_sizes$V1)]<-23

exome_sizes$V2<-as.numeric(exome_sizes$V2)
exome_sizes$V3<-as.numeric(exome_sizes$V3)
exome_sizes<-exome_sizes[complete.cases(exome_sizes),]
colnames(exome_sizes)<-c("chrom", "chromStart", "chromEnd")

genelengths<-data.frame(chrom=epigenetic_data$chr[epigenetic_data$gene %in% tcga_genelist], chromStart=epigenetic_data$start[epigenetic_data$gene %in% tcga_genelist], chromEnd=epigenetic_data$end[epigenetic_data$gene %in% tcga_genelist])

setDT(exome_sizes)
setDT(genelengths)
setkey(genelengths)
#if any exon overlaps with gene coordinates, then join them
annotated_exons<-foverlaps(exome_sizes, genelengths, type="any", nomatch=0L)

genelengths<-data.frame(chrom=epigenetic_data$chr[epigenetic_data$gene %in% tcga_genelist], chromStart=epigenetic_data$start[epigenetic_data$gene %in% tcga_genelist], chromEnd=epigenetic_data$end[epigenetic_data$gene %in% tcga_genelist], gene=epigenetic_data$gene[epigenetic_data$gene %in% tcga_genelist])

#full join exon coordinates and gene names
exon_annotations<-full_join(annotated_exons, genelengths, by=c("chromStart")) %>% select(i.chromStart, i.chromEnd, gene)

#genes listed by their exon sizes
exon_sizes<- exon_annotations %>% mutate(length=i.chromEnd-i.chromStart) %>% group_by(gene) %>% summarize(size=sum(length)) %>% as.data.frame()

#edit exon_sizes to only include genes in cluster_space
exon_sizes<-exon_sizes[exon_sizes$gene %in% rownames(cluster_space),]


############
#NMI Calculator

calc_NMI<-function(cluster_assignments, standardize=TRUE){
  
  clusters<-data.frame(gene=rownames(cluster_space), assignments=cluster_assignments)
  
  cluster_sizes<-setNames(data.frame(cluster= cluster_assignments, lengths=exon_sizes$size[match(rownames(cluster_space), exon_sizes$gene)]), c("cluster", "lengths")) %>% group_by(cluster) %>% summarize(length=sum(lengths))
  
  tcga_clust<-left_join(clusters, tcga[MS == "Non-hypermutated" & !is.na(Cancer_Code),], by=c("gene"="Hugo_Symbol"))
  colnames(tcga_clust)<-c("gene", "Cluster", "patient_id", "Variant", "Cancer_Code", "MS")
  setDT(tcga_clust)
  
  #filter by nonhypermutated mutation signature and non-NA Cancer Code
  tcga_v_f <- tcga_clust[MS == "Non-hypermutated" & !is.na(Cancer_Code)][
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
  
  tcga_newprob_given_cancer <- suppressWarnings(tcga_v_f[,
                                                         # Calculate Good Turing probabilities of
                                                         # at least one new variants per CLUSTER & cancer type
                                                         {
                                                           GT_probs <- goodturing_probs(
                                                             counts = v_f,
                                                             m = n_tumor[1]
                                                           )
                                                           .(p_atleast_1new_per10kb = GT_probs['atleast_1new'])
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
                                                    value.var = "p_atleast_1new_per10kb",
                                                    fill = 1 - exp(- 1/(length(unique(tcga$patient_id)) + 1))
                                                  ) %>%
                                                  # put  Clusters as row names
                                                  # and convert the data table into a matrix
                                                  magrittr::set_rownames(.$Cluster) %>%
                                                  .[, Cluster := NULL] %>%
                                                  data.matrix())
  
  #cancer probs are agnostic to metagene size... correct for metagene size
  if(standardize==TRUE){tcga_newprob_given_cancer_per10kb<-tcga_newprob_given_cancer *1E4/cluster_sizes$length}
  else(tcga_newprob_given_cancer_per10kb=tcga_newprob_given_cancer_per10kb)
  
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
    tcga_newprob_given_cancer_per10kb,
    cancer_prob,
    binary_minfo = FALSE,
    normalize = TRUE
  )
  
  # see the top 20
  return(nmi)
}

###########
#K-means with perturbation probs

kmeans_perturb<-function(x, centers, nstart, iter.max, perturb_prob){
  k.res<-kmeans(x=x, centers=centers, nstart = nstart, iter.max=iter.max)
  random_probs=runif(dim(x)[1], 0, 1)
  tmp<-cbind(k.res$cluster, random_probs)
  assignments<-ifelse(tmp[,2]<perturb_prob, sample.int(centers, sum(tmp[,2]<perturb_prob), replace=TRUE), tmp[,1])
  return(assignments)
}

#rand index calculator
adjRRand<-function(trcl, prcl){unname(RRand(trcl, prcl)[2])}

######
#k=50 clusters

cluster_membership<-data.frame(matrix(nrow=dim(cluster_space)[1], ncol=15))

cluster_membership<-lapply(X=1:15, FUN= function(X) kmeans_perturb(cluster_space, 50, nstart = 1, iter.max=20, perturb_prob=0.5))

#read list into dataframe by column, producing a dataframe of 500 columns and 16000 rows
df50 <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)

rownames(df50)<-rownames(cluster_space)

#calculate rand indices
rand_indices<-unname(unlist(apply(df50, adjRRand, trcl=df50[,1], MARGIN=2)))

dt50<-as.data.table(df50)
nmi_50<-dt50[ , lapply(.SD, "calc_NMI"), .SDcols = 1:15]

