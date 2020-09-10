####k=10
#0 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 10, nstart = 1, iter.max=20, perturb_prob=0))
df10_0perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k10_0perturb<-suppressWarnings(apply(df10_0perturb, MARGIN=2, FUN=calc_NMI))
#0.1 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 10, nstart = 1, iter.max=20, perturb_prob=0.1))
df10_0.1perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k10_0.1perturb<-suppressWarnings(apply(df10_0.1perturb, MARGIN=2, FUN=calc_NMI))
#0.25 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 10, nstart = 1, iter.max=20, perturb_prob=0.25))
df10_0.25perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k10_0.25perturb<-suppressWarnings(apply(df10_0.25perturb, MARGIN=2, FUN=calc_NMI))
#0.5 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 10, nstart = 1, iter.max=20, perturb_prob=0.5))
df10_0.5perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k10_0.5perturb<-suppressWarnings(apply(df10_0.5perturb, MARGIN=2, FUN=calc_NMI))
#1 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 10, nstart = 1, iter.max=20, perturb_prob=1))
df10_1perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k10_1perturb<-suppressWarnings(apply(df10_1perturb, MARGIN=2, FUN=calc_NMI))

#####k=20
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 20, nstart = 1, iter.max=20, perturb_prob=0))
df20_0perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k20_0perturb<-suppressWarnings(apply(df20_0perturb, MARGIN=2, FUN=calc_NMI))
#0.1 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 20, nstart = 1, iter.max=20, perturb_prob=0.1))
df20_0.1perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k20_0.1perturb<-suppressWarnings(apply(df20_0.1perturb, MARGIN=2, FUN=calc_NMI))
#0.25 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 20, nstart = 1, iter.max=20, perturb_prob=0.25))
df20_0.25perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k20_0.25perturb<-suppressWarnings(apply(df20_0.25perturb, MARGIN=2, FUN=calc_NMI))
#0.5 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 20, nstart = 1, iter.max=20, perturb_prob=0.5))
df20_0.5perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k20_0.5perturb<-suppressWarnings(apply(df20_0.5perturb, MARGIN=2, FUN=calc_NMI))
#1 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 20, nstart = 1, iter.max=20, perturb_prob=1))
df20_1perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k20_1perturb<-suppressWarnings(apply(df20_1perturb, MARGIN=2, FUN=calc_NMI))

#####k=30
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 20, nstart = 1, iter.max=20, perturb_prob=0))
df30_0perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k30_0perturb<-suppressWarnings(apply(df30_0perturb, MARGIN=2, FUN=calc_NMI))
#0.1 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 20, nstart = 1, iter.max=20, perturb_prob=0.1))
df30_0.1perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k30_0.1perturb<-suppressWarnings(apply(df30_0.1perturb, MARGIN=2, FUN=calc_NMI))
#0.25 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 30, nstart = 1, iter.max=20, perturb_prob=0.25))
df30_0.25perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k30_0.25perturb<-suppressWarnings(apply(df30_0.25perturb, MARGIN=2, FUN=calc_NMI))
#0.5 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 30, nstart = 1, iter.max=20, perturb_prob=0.5))
df30_0.5perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k30_0.5perturb<-suppressWarnings(apply(df30_0.5perturb, MARGIN=2, FUN=calc_NMI))
#1 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 30, nstart = 1, iter.max=20, perturb_prob=1))
df30_1perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k30_1perturb<-suppressWarnings(apply(df30_1perturb, MARGIN=2, FUN=calc_NMI))

#####k=40
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 40, nstart = 1, iter.max=20, perturb_prob=0))
df40_0perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k40_0perturb<-suppressWarnings(apply(df40_0perturb, MARGIN=2, FUN=calc_NMI))
#0.1 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 40, nstart = 1, iter.max=20, perturb_prob=0.1))
df40_0.1perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k40_0.1perturb<-suppressWarnings(apply(df40_0.1perturb, MARGIN=2, FUN=calc_NMI))
#0.25 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 40, nstart = 1, iter.max=20, perturb_prob=0.25))
df40_0.25perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k40_0.25perturb<-suppressWarnings(apply(df40_0.25perturb, MARGIN=2, FUN=calc_NMI))
#0.5 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 40, nstart = 1, iter.max=20, perturb_prob=0.5))
df40_0.5perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k40_0.5perturb<-suppressWarnings(apply(df40_0.5perturb, MARGIN=2, FUN=calc_NMI))
#1 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 40, nstart = 1, iter.max=20, perturb_prob=1))
df40_1perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k40_1perturb<-suppressWarnings(apply(df40_1perturb, MARGIN=2, FUN=calc_NMI))

#####k=50
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 50, nstart = 1, iter.max=20, perturb_prob=0))
df50_0perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k50_0perturb<-suppressWarnings(apply(df50_0perturb, MARGIN=2, FUN=calc_NMI))
#0.1 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 50, nstart = 1, iter.max=20, perturb_prob=0.1))
df50_0.1perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k50_0.1perturb<-suppressWarnings(apply(df50_0.1perturb, MARGIN=2, FUN=calc_NMI))
#0.25 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 50, nstart = 1, iter.max=20, perturb_prob=0.25))
df50_0.25perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k50_0.25perturb<-suppressWarnings(apply(df50_0.25perturb, MARGIN=2, FUN=calc_NMI))
#0.5 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 50, nstart = 1, iter.max=20, perturb_prob=0.5))
df50_0.5perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k50_0.5perturb<-suppressWarnings(apply(df50_0.5perturb, MARGIN=2, FUN=calc_NMI))
#1 perturb
cluster_membership<-lapply(X=1:500, FUN= function(X) kmeans_perturb(cluster_space, 50, nstart = 1, iter.max=20, perturb_prob=1))
df50_1perturb <- data.frame(matrix(unlist(cluster_membership), nrow=dim(cluster_space)[1], byrow=FALSE),stringsAsFactors=FALSE)
rm(cluster_membership)
NMI_k50_1perturb<-suppressWarnings(apply(df50_1perturb, MARGIN=2, FUN=calc_NMI))


######NMI values

to_plot<-rbind(melt(NMI_k10_0perturb) %>% mutate(k="k=10") %>% mutate("Perturb"=0),
               melt(NMI_k10_0.1perturb) %>% mutate(k="k=10") %>% mutate("Perturb"=0.1),
               melt(NMI_k10_0.25perturb) %>% mutate(k="k=10") %>% mutate("Perturb"=0.25),
               melt(NMI_k10_0.5perturb) %>% mutate(k="k=10") %>% mutate("Perturb"=0.5),
               melt(NMI_k10_1perturb) %>% mutate(k="k=10") %>% mutate("Perturb"=1),
               melt(NMI_k20_0perturb) %>% mutate(k="k=20") %>% mutate("Perturb"=0),
               melt(NMI_k20_0.1perturb) %>% mutate(k="k=30") %>% mutate("Perturb"=0.1),
               melt(NMI_k20_0.25perturb) %>% mutate(k="k=30") %>% mutate("Perturb"=0.25),
               melt(NMI_k20_0.5perturb) %>% mutate(k="k=30") %>% mutate("Perturb"=0.5),
               melt(NMI_k20_1perturb) %>% mutate(k="k=30") %>% mutate("Perturb"=1),
               melt(NMI_k30_0perturb) %>% mutate(k="k=30") %>% mutate("Perturb"=0),
               melt(NMI_k30_0.1perturb) %>% mutate(k="k=30") %>% mutate("Perturb"=0.1),
               melt(NMI_k30_0.25perturb) %>% mutate(k="k=30") %>% mutate("Perturb"=0.25),
               melt(NMI_k30_0.5perturb) %>% mutate(k="k=30") %>% mutate("Perturb"=0.5),
               melt(NMI_k30_1perturb) %>% mutate(k="k=30") %>% mutate("Perturb"=1),
               melt(NMI_k40_0perturb) %>% mutate(k="k=40") %>% mutate("Perturb"=0),
               melt(NMI_k40_0.1perturb) %>% mutate(k="k=40") %>% mutate("Perturb"=0.1),
               melt(NMI_k40_0.25perturb) %>% mutate(k="k=40") %>% mutate("Perturb"=0.25),
               melt(NMI_k40_0.5perturb) %>% mutate(k="k=40") %>% mutate("Perturb"=0.5),
               melt(NMI_k40_1perturb) %>% mutate(k="k=40") %>% mutate("Perturb"=1),
               melt(NMI_k50_0perturb) %>% mutate(k="k=50") %>% mutate("Perturb"=0),
               melt(NMI_k50_0.1perturb) %>% mutate(k="k=50") %>% mutate("Perturb"=0.1),
               melt(NMI_k50_0.25perturb) %>% mutate(k="k=50") %>% mutate("Perturb"=0.25),
               melt(NMI_k50_0.5perturb) %>% mutate(k="k=50") %>% mutate("Perturb"=0.5),
               melt(NMI_k50_1perturb) %>% mutate(k="k=50") %>% mutate("Perturb"=1))

########Rand
rand_10_0perturb<-unname(unlist(apply(df10_0perturb, adjRRand, trcl=df10_0perturb[,1], MARGIN=2)))
rand_10_0.1perturb<-unname(unlist(apply(df10_0.1perturb, adjRRand, trcl=df10_0.1perturb[,1], MARGIN=2)))
rand_10_0.25perturb<-unname(unlist(apply(df10_0.25perturb, adjRRand, trcl=df10_0.25perturb[,1], MARGIN=2)))
rand_10_0.5perturb<-unname(unlist(apply(df10_0.5perturb, adjRRand, trcl=df10_0.5perturb[,1], MARGIN=2)))
rand_10_1perturb<-unname(unlist(apply(df10_1perturb, adjRRand, trcl=df10_1perturb[,1], MARGIN=2)))
rand_20_0perturb<-unname(unlist(apply(df20_0perturb, adjRRand, trcl=df20_0perturb[,1], MARGIN=2)))
rand_20_0.1perturb<-unname(unlist(apply(df20_0.1perturb, adjRRand, trcl=df20_0.1perturb[,1], MARGIN=2)))
rand_20_0.25perturb<-unname(unlist(apply(df20_0.25perturb, adjRRand, trcl=df20_0.25perturb[,1], MARGIN=2)))
rand_20_0.5perturb<-unname(unlist(apply(df20_0.5perturb, adjRRand, trcl=df20_0.5perturb[,1], MARGIN=2)))
rand_20_1perturb<-unname(unlist(apply(df20_1perturb, adjRRand, trcl=df20_1perturb[,1], MARGIN=2)))
rand_30_0perturb<-unname(unlist(apply(df30_0perturb, adjRRand, trcl=df30_0perturb[,1], MARGIN=2)))
rand_30_0.1perturb<-unname(unlist(apply(df30_0.1perturb, adjRRand, trcl=df30_0.1perturb[,1], MARGIN=2)))
rand_30_0.25perturb<-unname(unlist(apply(df30_0.25perturb, adjRRand, trcl=df30_0.25perturb[,1], MARGIN=2)))
rand_30_0.5perturb<-unname(unlist(apply(df30_0.5perturb, adjRRand, trcl=df30_0.5perturb[,1], MARGIN=2)))
rand_30_1perturb<-unname(unlist(apply(df30_1perturb, adjRRand, trcl=df30_1perturb[,1], MARGIN=2)))
rand_40_0perturb<-unname(unlist(apply(df40_0perturb, adjRRand, trcl=df40_0perturb[,1], MARGIN=2)))
rand_40_0.1perturb<-unname(unlist(apply(df40_0.1perturb, adjRRand, trcl=df40_0.1perturb[,1], MARGIN=2)))
rand_40_0.25perturb<-unname(unlist(apply(df40_0.25perturb, adjRRand, trcl=df40_0.25perturb[,1], MARGIN=2)))
rand_40_0.5perturb<-unname(unlist(apply(df40_0.5perturb, adjRRand, trcl=df40_0.5perturb[,1], MARGIN=2)))
rand_40_1perturb<-unname(unlist(apply(df40_1perturb, adjRRand, trcl=df40_1perturb[,1], MARGIN=2)))
rand_50_0perturb<-unname(unlist(apply(df50_0perturb, adjRRand, trcl=df50_0perturb[,1], MARGIN=2)))
rand_50_0.1perturb<-unname(unlist(apply(df50_0.1perturb, adjRRand, trcl=df50_0.1perturb[,1], MARGIN=2)))
rand_50_0.25perturb<-unname(unlist(apply(df50_0.25perturb, adjRRand, trcl=df50_0.25perturb[,1], MARGIN=2)))
rand_50_0.5perturb<-unname(unlist(apply(df50_0.5perturb, adjRRand, trcl=df50_0.5perturb[,1], MARGIN=2)))
rand_50_1perturb<-unname(unlist(apply(df50_1perturb, adjRRand, trcl=df50_1perturb[,1], MARGIN=2)))

#########Size distributions

get_cluster_sizes<-function(assignments){
  clusters<-data.frame(gene=rownames(cluster_space), assignments=assignments)
  
  cluster_sizes<-setNames(data.frame(cluster= cluster_assignments, lengths=exon_sizes$size[match(rownames(cluster_space), exon_sizes$gene)]), c("cluster", "lengths")) %>% group_by(cluster) %>% summarize(length=sum(lengths))
  return(cluster_sizes %>% as.data.frame())
}

sizes_10_0perturb<-apply(df10_0perturb, get_cluster_sizes, MARGIN=2)
sizes_10_0.1perturb<-apply(df10_0.1perturb, get_cluster_sizes, MARGIN=2)
sizes_10_0.25perturb<-apply(df10_0.25perturb, get_cluster_sizes, MARGIN=2)
sizes_10_0.5perturb<-apply(df10_0.5perturb, get_cluster_sizes, MARGIN=2)
sizes_10_1perturb<-apply(df10_1perturb, get_cluster_sizes, MARGIN=2)
sizes_20_0perturb<-apply(df20_0perturb, get_cluster_sizes, MARGIN=2)
sizes_20_0.1perturb<-apply(df20_0.1perturb, get_cluster_sizes, MARGIN=2)
sizes_20_0.25perturb<-apply(df20_0.25perturb, get_cluster_sizes, MARGIN=2)
sizes_20_0.5perturb<-apply(df20_0.5perturb, get_cluster_sizes, MARGIN=2)
sizes_20_1perturb<-apply(df20_1perturb, get_cluster_sizes, MARGIN=2)
sizes_30_0perturb<-apply(df30_0perturb, get_cluster_sizes, MARGIN=2)
sizes_30_0.1perturb<-apply(df30_0.1perturb, get_cluster_sizes, MARGIN=2)
sizes_30_0.25perturb<-apply(df30_0.25perturb, get_cluster_sizes, MARGIN=2)
sizes_30_0.5perturb<-apply(df30_0.5perturb, get_cluster_sizes, MARGIN=2)
sizes_30_1perturb<-apply(df30_1perturb, get_cluster_sizes, MARGIN=2)
sizes_40_0perturb<-apply(df40_0perturb, get_cluster_sizes, MARGIN=2)
sizes_40_0.1perturb<-apply(df40_0.1perturb, get_cluster_sizes, MARGIN=2)
sizes_40_0.25perturb<-apply(df40_0.25perturb, get_cluster_sizes, MARGIN=2)
sizes_40_0.5perturb<-apply(df40_0.5perturb, get_cluster_sizes, MARGIN=2)
sizes_40_1perturb<-apply(df40_1perturb, get_cluster_sizes, MARGIN=2)
sizes_50_0perturb<-apply(df50_0perturb, get_cluster_sizes, MARGIN=2)
sizes_50_0.1perturb<-apply(df50_0.1perturb, get_cluster_sizes, MARGIN=2)
sizes_50_0.25perturb<-apply(df50_0.25perturb, get_cluster_sizes, MARGIN=2)
sizes_50_0.5perturb<-apply(df50_0.5perturb, get_cluster_sizes, MARGIN=2)
sizes_50_1perturb<-apply(df50_1perturb, get_cluster_sizes, MARGIN=2)