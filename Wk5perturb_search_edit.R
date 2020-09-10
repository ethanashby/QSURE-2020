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
library(e1071)
library(drc)
library(DescTools)
library(lattice)
library(phyclust)
library(Rcpp)
source("calc_minfo.R")

data(tcga)
tcga_genelist<-tcga[MS == "Non-hypermutated" & !is.na(Cancer_Code),]$Hugo_Symbol %>% unique()

#########
#Build covariate table of metafeatures

#read in covariates and remove empty rows
epigenetic_data<-read.table("covariate-table.csv", header=TRUE, sep=",")
colnames(epigenetic_data)[1]<-"gene"
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
exon_annotations<-full_join(annotated_exons, genelengths, by=c("chromStart"))
exon_annotations<-exon_annotations %>% dplyr::select(i.chromStart, i.chromEnd, gene)

#genes listed by their exon sizes
exon_sizes<- exon_annotations %>% mutate(length=i.chromEnd-i.chromStart) %>% group_by(gene) %>% summarize(size=sum(length)) %>% as.data.frame()

#edit exon_sizes to only include genes in cluster_space
exon_sizes<-exon_sizes[exon_sizes$gene %in% rownames(cluster_space),]



###################
#Functions





###########
#TCGA Calculator

calc_GT<-function(cluster_assignments, standardize=TRUE){
  
  cluster_assignments<-df50_0perturb[,1]
  
  clusters<-data.frame(gene=rownames(cluster_space), assignments=cluster_assignments)
  
  cluster_sizes<-setNames(data.frame(cluster= cluster_assignments, lengths=exon_sizes$size[match(rownames(cluster_space), exon_sizes$gene)]), c("cluster", "lengths")) %>% group_by(cluster) %>% summarize(length=sum(lengths))
  
  tcga_clust<-left_join(clusters, tcga[MS == "Non-hypermutated" & !is.na(Cancer_Code),], by=c("gene"="Hugo_Symbol"))
  colnames(tcga_clust)<-c("gene", "Cluster", "patient_id", "Variant", "Cancer_Code", "MS")
  setDT(tcga_clust)
  
  tmp<-left_join(tcga_clust, cluster_sizes, by=c("Cluster"="cluster"))
  
  #filter by nonhypermutated mutation signature and non-NA Cancer Code
  tcga_v_f <- tmp[MS == "Non-hypermutated" & !is.na(Cancer_Code)][
    ,
    # number of tumors per cancer type
    n_tumor := length(unique(patient_id))
    ][,
      # variant frequencies by CLUSTER, cancer type
      # also save n_tumor per cancer
      .(v_f = length(unique(patient_id)),
        n_tumor = n_tumor[1]),
      by = .(Variant, Cluster)
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
                                                         by = .(Cluster)
                                                         ] %>%
                                                  # cast the long data frame into wide format -- fill missing values
                                                  # (i.e. no mutation for that (CLUSTER, cancer) pair)
                                                  # with ZERO
                                                  # replace fill = 0 by
                                                  # fill = 1 - exp(- 1/(length(unique(tcga$patient_id)) + 1))
                                                  # to exactly recreate Fig 3b in Nat. Comm. paper
                                                  # put  Clusters as row names
                                                  # and convert the data table into a matrix
                                                  magrittr::set_rownames(.$Cluster) %>%
                                                  .[, Cluster := NULL] %>%
                                                  data.matrix())
  
  #cancer probs are agnostic to metagene size... correct for metagene size
  GT_probs<-tcga_newprob_given_cancer[order(as.numeric(rownames(tcga_newprob_given_cancer))),]
  if(standardize==TRUE){tcga_newprob_given_cancerper10kb<-GT_probs *1E4/cluster_sizes$length}
  else(tcga_newprob_given_cancer_per10kb=tcga_newprob_given_cancer_per10kb)
  return(tcga_newprob_given_cancerper10kb)
}

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








########################
#perturbation assessment



######
#k=50 clusters 0 perturb

cluster_membership<-data.frame(matrix(nrow=dim(cluster_space)[1], ncol=500))

cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 50, nstart = 1, iter.max=20, perturb_prob=0))

#read list into dataframe by column, producing a dataframe of 500 columns and 16000 rows
df50_0perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)

rownames(df50_0perturb)<-rownames(cluster_space)

#calculate rand indices
rand_indices_0perturb<-unname(unlist(apply(df50_0perturb, adjRRand, trcl=df50_0perturb[,1], MARGIN=2)))

dt50_0perturb<-as.data.table(df50_0perturb)
nmi_50_0perturb<-dt50_0perturb[ , lapply(.SD, "calc_NMI"), .SDcols = 1:500]

nmi_50_0perturb_melt<-nmi_50_0perturb %>% melt()
nmi_50_0perturb_melt$variable<-"0 perturb"
nmi_50_0perturb_melt %>% ggplot()+geom_boxplot(aes(x=variable, y=value))

#########
#k=50 0.1 perturb

cluster_membership<-data.frame(matrix(nrow=dim(cluster_space)[1], ncol=500))

cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 50, nstart = 1, iter.max=20, perturb_prob=0.1))

#read list into dataframe by column, producing a dataframe of 500 columns and 16000 rows
df50_0.1perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)

rownames(df50_0.1perturb)<-rownames(cluster_space)

#calculate rand indices
rand_indices_0.1perturb<-unname(unlist(apply(df50_0.1perturb, adjRRand, trcl=df50_0.1perturb[,1], MARGIN=2)))

dt50_0.1perturb<-as.data.table(df50_0.1perturb)
nmi_50_0.1perturb<-dt50_0.1perturb[ , lapply(.SD, "calc_NMI"), .SDcols = 1:500]

nmi_50_0.1perturb_melt<-nmi_50_0.1perturb %>% melt()

#########
#k=50 0.25 perturb

cluster_membership<-data.frame(matrix(nrow=dim(cluster_space)[1], ncol=500))

cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 50, nstart = 1, iter.max=20, perturb_prob=0.25))

#read list into dataframe by column, producing a dataframe of 500 columns and 16000 rows
df50_0.25perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)

rownames(df50_0.25perturb)<-rownames(cluster_space)

#calculate rand indices
rand_indices_0.25perturb<-unname(unlist(apply(df50_0.25perturb, adjRRand, trcl=df50_0.25perturb[,1], MARGIN=2)))

dt50_0.25perturb<-as.data.table(df50_0.25perturb)
nmi_50_0.25perturb<-dt50_0.25perturb[ , lapply(.SD, "calc_NMI"), .SDcols = 1:500]

nmi_50_0.25perturb_melt<-nmi_50_0.25perturb %>% melt()

#########
#k=50 0.5 perturb

cluster_membership<-data.frame(matrix(nrow=dim(cluster_space)[1], ncol=500))

cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 50, nstart = 1, iter.max=20, perturb_prob=0.5))

#read list into dataframe by column, producing a dataframe of 500 columns and 16000 rows
df50_0.5perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)

rownames(df50_0.5perturb)<-rownames(cluster_space)

#calculate rand indices
rand_indices_0.5perturb<-unname(unlist(apply(df50_0.5perturb, adjRRand, trcl=df50_0.5perturb[,1], MARGIN=2)))

dt50_0.5perturb<-as.data.table(df50_0.5perturb)
nmi_50_0.5perturb<-dt50_0.5perturb[ , lapply(.SD, "calc_NMI"), .SDcols = 1:500]

nmi_50_0.5perturb_melt<-nmi_50_0.5perturb %>% melt()

#########
#k=50 0.75 perturb

cluster_membership<-data.frame(matrix(nrow=dim(cluster_space)[1], ncol=500))

cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 50, nstart = 1, iter.max=20, perturb_prob=0.75))

#read list into dataframe by column, producing a dataframe of 500 columns and 16000 rows
df50_0.75perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)

rownames(df50_0.75perturb)<-rownames(cluster_space)

#calculate rand indices
rand_indices_0.75perturb<-unname(unlist(apply(df50_0.75perturb, adjRRand, trcl=df50_0.75perturb[,1], MARGIN=2)))

dt50_0.75perturb<-as.data.table(df50_0.75perturb)
nmi_50_0.75perturb<-dt50_0.75perturb[ , lapply(.SD, "calc_NMI"), .SDcols = 1:500]

nmi_50_0.75perturb_melt<-nmi_50_0.75perturb %>% melt()

#########
#k=50 1 perturb

cluster_membership<-data.frame(matrix(nrow=dim(cluster_space)[1], ncol=500))

cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 50, nstart = 1, iter.max=20, perturb_prob=1))

#read list into dataframe by column, producing a dataframe of 500 columns and 16000 rows
df50_1perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)

rownames(df50_1perturb)<-rownames(cluster_space)

#calculate rand indices
rand_indices_1perturb<-unname(unlist(apply(df50_1perturb, adjRRand, trcl=df50_1perturb[,1], MARGIN=2)))

dt50_1perturb<-as.data.table(df50_1perturb)
nmi_50_1perturb<-dt50_1perturb[ , lapply(.SD, "calc_NMI"), .SDcols = 1:500]

nmi_50_1perturb_melt<-nmi_50_1perturb %>% melt()

#############
#IMPACT NMI

data(impact)
setDT(impact)

#list of impact genes
impact_genelist <- impact[MS=="Non-hypermutated" & !is.na(Cancer_Type)] %>% dplyr::select(Hugo_Symbol) %>% unique() %>% unlist() %>% unname()

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
exon_annotations<-full_join(annotated_exons, genelengths, by=c("chromStart")) %>% dplyr::select(i.chromStart, i.chromEnd, gene)

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

##########
#plot

nmi_IMPACT_melt$variable<-str_wrap("IMPACT No Pooling", width=8)
nmi_IMPACT_melt<-nmi_IMPACT_melt[c("variable", "value")]
nmi_50_0perturb_melt$variable<-str_wrap("Perturb Prob 0", width=8)
nmi_50_0.1perturb_melt$variable<-str_wrap("Perturb Prob 0.1", width=8)
nmi_50_0.25perturb_melt$variable<-str_wrap("Perturb Prob 0.25", width=8)
nmi_50_0.5perturb_melt$variable<-str_wrap("Perturb Prob 0.50", width=8)
nmi_50_0.75perturb_melt$variable<-str_wrap("Perturb Prob 0.75", width=8)
nmi_50_1perturb_melt$variable<-str_wrap("Perturb Prob 1", width=8)

to_plot<-rbind(nmi_IMPACT_melt, nmi_50_0perturb_melt, nmi_50_0.1perturb_melt, nmi_50_0.25perturb_melt, nmi_50_0.5perturb_melt, nmi_50_0.75perturb_melt, nmi_50_1perturb_melt)

######
#3 boxplot varieties

ggplot(to_plot)+geom_boxplot(aes(x=variable, y=value, color=variable))+theme_bw()+ylab("NMI per 10kb Exonic Length")+labs(color="Cluster Perturbation")+xlab(element_blank())+theme(axis.text.x=element_text(size=24),axis.text.y=element_text(size=16), axis.title=element_text(size=24, face="bold"), legend.position="none")

ggplot(to_plot)+geom_boxplot(aes(x=variable, y=log10(value), color=variable))+theme_bw()+ylab("Log NMI per 10kb Exonic Length")+labs(color="Cluster Perturbation")+xlab(element_blank())+theme(axis.text.x=element_text(size=24),axis.text.y=element_text(size=16), axis.title=element_text(size=24, face="bold"), legend.position="none")

ggplot(to_plot %>% dplyr::filter(variable!="IMPACT\nNo\nPooling"))+geom_boxplot(aes(x=variable, y=value, color=variable))+theme_bw()+ylab("NMI per 10kb Exonic Length")+labs(color="Cluster Perturbation")+xlab(element_blank())+theme(axis.text.x=element_text(size=24),axis.text.y=element_text(size=16), axis.title=element_text(size=24, face="bold"), legend.position="none")

######
#rand plots
rand<-data.frame("Rand"=c(rand_indices_0perturb, rand_indices_0.1perturb, rand_indices_0.25perturb, rand_indices_0.5perturb, rand_indices_0.75perturb, rand_indices_1perturb))
rand$id<-c(rep("Perturb Prob 0", 500), rep("Perturb Prob 0.1", 500), rep("Perturb Prob 0.25", 500), rep("Perturb Prob 0.5", 500), rep("Perturb Prob 0.75", 500), rep("Perturb Prob 1", 500))
  
ggplot(rand)+geom_density(aes(x=Rand, group=id, fill=id), alpha=0.5, adjust=2)+facet_wrap(~id, scales="free")+theme_bw()+theme(legend.position="none")                  




#############
#Metagene Size

melt0perturb<-df50_0perturb %>% melt()
melt0perturb<-melt0perturb %>% group_by(variable, value) %>% summarize(count=n())
melt0perturb$id<-"Perturb Prob 0"

melt0.1perturb<-df50_0.1perturb %>% melt()
melt0.1perturb<-melt0.1perturb %>% group_by(variable, value) %>% summarize(count=n())
melt0.1perturb$id<-"Perturb Prob 0.1"

melt0.25perturb<-df50_0.25perturb %>% melt()
melt0.25perturb<-melt0.25perturb %>% group_by(variable, value) %>% summarize(count=n())
melt0.25perturb$id<-"Perturb Prob 0.25"

melt0.5perturb<-df50_0.5perturb %>% melt()
melt0.5perturb<-melt0.5perturb %>% group_by(variable, value) %>% summarize(count=n())
melt0.5perturb$id<-"Perturb Prob 0.5"

melt0.75perturb<-df50_0.75perturb %>% melt()
melt0.75perturb<-melt0.75perturb %>% group_by(variable, value) %>% summarize(count=n())
melt0.75perturb$id<-"Perturb Prob 0.75"

melt1perturb<-df50_1perturb %>% melt()
melt1perturb<-melt1perturb %>% group_by(variable, value) %>% summarize(count=n())
melt1perturb$id<-"Perturb Prob 1"

cluster_sizes<-rbind(melt0perturb, melt0.1perturb, melt0.25perturb, melt0.5perturb, melt0.75perturb, melt1perturb)

ggplot(cluster_sizes)+geom_histogram(aes(x=count, fill=id), alpha=0.8)+facet_wrap(~id)+theme_bw()+xlab("Metagene Cluster Size")+theme(legend.position="none", strip.text.x=element_text(size=22), axis.text = element_text(size=16), axis.title=element_text(size=22))


###With scaled Good-Turing probs included

GT_50_0perturb<-dt50_0perturb[ , lapply(.SD, "calc_GT"), .SDcols = 1:500]
GT_50_0perturb_melt<-GT_50_0perturb %>% melt()

GT_50_0.1perturb<-dt50_0.1perturb[ , lapply(.SD, "calc_GT"), .SDcols = 1:500]
GT_50_0.1perturb_melt<-GT_50_0.1perturb %>% melt()

GT_50_0.25perturb<-dt50_0.25perturb[ , lapply(.SD, "calc_GT"), .SDcols = 1:500]
GT_50_0.25perturb_melt<-GT_50_0.25perturb %>% melt()

GT_50_0.5perturb<-dt50_0.5perturb[ , lapply(.SD, "calc_GT"), .SDcols = 1:500]
GT_50_0.5perturb_melt<-GT_50_0.5perturb %>% melt()

GT_50_0.75perturb<-dt50_0.75perturb[ , lapply(.SD, "calc_GT"), .SDcols = 1:500]
GT_50_0.75perturb_melt<-GT_50_0.75perturb %>% melt()

GT_50_1perturb<-dt50_1perturb[ , lapply(.SD, "calc_GT"), .SDcols = 1:500]
GT_50_1perturb_melt<-GT_50_1perturb %>% melt()

GT_50_0perturb_melt$clust<-rep(1:50, 500)
GT_50_0.1perturb_melt$clust<-rep(1:50, 500)
GT_50_0.25perturb_melt$clust<-rep(1:50, 500)
GT_50_0.5perturb_melt$clust<-rep(1:50, 500)
GT_50_0.75perturb_melt$clust<-rep(1:50, 500)
GT_50_1perturb_melt$clust<-rep(1:50, 500)

GT_0perturb_plot<-left_join(GT_50_0perturb_melt, melt0perturb, by=c("variable"="variable", "clust"="value"))
GT_0.1perturb_plot<-left_join(GT_50_0.1perturb_melt, melt0.1perturb, by=c("variable"="variable", "clust"="value"))
GT_0.25perturb_plot<-left_join(GT_50_0.25perturb_melt, melt0.25perturb, by=c("variable"="variable", "clust"="value"))
GT_0.5perturb_plot<-left_join(GT_50_0.5perturb_melt, melt0.5perturb, by=c("variable"="variable", "clust"="value"))
GT_0.75perturb_plot<-left_join(GT_50_0.75perturb_melt, melt0.75perturb, by=c("variable"="variable", "clust"="value"))
GT_0.75perturb_plot<-left_join(GT_50_0.75perturb_melt, melt0.75perturb, by=c("variable"="variable", "clust"="value"))
GT_1perturb_plot<-left_join(GT_50_1perturb_melt, melt1perturb, by=c("variable"="variable", "clust"="value"))

GT_size_plot<-rbind(GT_0perturb_plot, GT_0.1perturb_plot, GT_0.25perturb_plot, GT_0.5perturb_plot, GT_0.75perturb_plot, GT_1perturb_plot)

GT_size_plot %>% ggplot()+stat_density_2d(aes(x=value, y=count, fill= ..level..), geom="polygon")+theme_bw()+facet_wrap(~id, scales="free_y")+scale_fill_distiller(palette="Spectral", direction=1)+theme_bw()+xlab("Good-Turing Probability per 10kb")+ylab("Cluster Size")+theme(legend.position="none", strip.text.x=element_text(size=22), axis.text = element_text(size=16), axis.title=element_text(size=22))

GT_size_plot %>% ggplot()+geom_density(aes(x=count))+facet_wrap(~id, scales="free")+theme_bw()+xlab("density")+xlab("Cluster Size")+theme(legend.position="none", strip.text.x=element_text(size=22), axis.text = element_text(size=16), axis.title=element_text(size=22))

GT_size_plot %>% ggplot()+geom_density(aes(x=value))+facet_wrap(~id, scales="free")+theme_bw()+xlab("density")+xlab("Cluster Size")+theme(legend.position="none", strip.text.x=element_text(size=22), axis.text = element_text(size=16), axis.title=element_text(size=22))

####save rdata
saveRDS(to_plot, "NMI_perturb_500iter.rds")
saveRDS(rand, "Rand_perturb_500iter.rds")
saveRDS(cluster_sizes, "Cluster_Sizes.rds")
saveRDS(GT_size_plot, "GT_Sizes.rds")
GT_size_plot<-readRDS("GT_Sizes.rds")

cut(GT_size_plot %>% filter(id=="Perturb Prob 0") %>% dplyr::select(count) %>% unlist(), breaks=seq(0, 500, 50)) %>% table()


#################################
#Clustering varying k high depth
#################################

#####
#k=50 5000 iter
#####

gc()
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 50, nstart = 1, iter.max=100, perturb_prob=0))
                           
#read list into dataframe by column, producing a dataframe of 5000 columns and ~16000 rows
df50_5k <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)

rownames(df50_5k)<-rownames(cluster_space)

#calculate rand indices
rand_indices_50<-unname(unlist(apply(df50_5k, adjRRand, trcl=df50_5k[,1], MARGIN=2)))

dt50_5k<-as.data.table(df50_5k)
nmi_50_5k<-dt50_5k[ , lapply(.SD, "calc_NMI"), .SDcols = 1:500]

nmi_50_5k_melt<-nmi_50_5k %>% melt()

#k=75 5000 iter

cluster_membership<-data.frame(matrix(nrow=dim(cluster_space)[1], ncol=500))

cluster_membership<-lapply(X=1:500, FUN= function(x) kmeans_perturb(cluster_space, 75, nstart = 1, iter.max=100, perturb_prob = 0))

#read list into dataframe by column, producing a dataframe of 5000 columns and ~16000 rows
df75_5k <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)

rownames(df75_5k)<-rownames(cluster_space)

#calculate rand indices
rand_indices_75<-unname(unlist(apply(df75_5k, adjRRand, trcl=df75_5k[,1], MARGIN=2)))

dt75_5k<-as.data.table(df75_5k)
nmi_75_5k<-dt75_5k[ , lapply(.SD, "calc_NMI"), .SDcols = 1:500]

nmi_75_5k_melt<-nmi_75_5k %>% melt()

#k=100 5000 iter

cluster_membership<-data.frame(matrix(nrow=dim(cluster_space)[1], ncol=500))

cluster_membership<-lapply(X=1:500, FUN= function(x) kmeans(cluster_space, 100, nstart = 1, iter.max=100))

#read list into dataframe by column, producing a dataframe of 5000 columns and ~16000 rows
df100_5k <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)

rownames(df100_5k)<-rownames(cluster_space)

#calculate rand indices
rand_indices_100<-unname(unlist(apply(df100_5k, adjRRand, trcl=df100_5k[,1], MARGIN=2)))

dt100_5k<-as.data.table(df100_5k)
nmi_100_5k<-dt100_5k[ , lapply(.SD, "calc_NMI"), .SDcols = 1:500]

nmi_100_5k_melt<-nmi_100_5k %>% melt()

nmi_50_5k_melt$id<-"k=50"
nmi_75_5k_melt$id<-"k=75"
tmp<-rbind(nmi_50_5k_melt, nmi_75_5k_melt)
saveRDS(tmp, "50_75_nmis.rds")

ggplot(tmp)+geom_boxplot(aes(x=id, y=value, color=id))
nmi_100_5k_melt$id<-"k=100"


