library("scde")
library("DBI")

setwd("./")

# input the read count table
count_tab <- read.table("summary_star_readcount",head=T,row.names =1)

# drop low express genes and low coverage samples
cd <- clean.counts(count_tab, min.lib.size = 500, min.reads = 10, min.detected = 5)
#> dim(cd)
#[1] 9021   48

# get color group flag  
#            red       green       skyblue       purple
l2cols <- c(rep("coral4",12),  rep("mediumpurple2",7))

# building error model
knn <- knn.error.models(cd, k = ncol(cd)/2, n.cores = 20, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 5, save.model.plots = TRUE, verbose = 1)  # turn on verbosity

# time consuming step, better to save the result
save.image("knn.RData")

# normalize out expected levels of technical and intrinsic biological noise
varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 20, plot = T)
# Controlling for sequencing depth, how many gene is present in each cell
varinfo <- pagoda.subtract.aspect(varinfo, colSums(cd[, rownames(knn)]>0))
save.image("varinfo.RData")


write.table(varinfo$mat,file="tab_varinfo_mat",quote=F,sep="\t")
write.table(varinfo$matw,file="tab_varinfo_matw",quote=F,sep="\t")
write.table(varinfo$arv,file="tab_varinfo_arv",quote=F,sep="\t")
write.table(varinfo$armodes,file="tab_varinfo_arv",quote=F,sep="\t")
write.table(sort(varinfo$arv, decreasing = TRUE)[1:1000],file="tab_top_1000_overdis_genes",sep="\t",quote=F)

# Evaluate overdispersion of pre-defined gene sets
# get GO annotations
library(org.Hs.eg.db)

ids <- unlist(lapply(mget(rownames(cd), org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
rids <- names(ids); names(rids) <- ids

go.env <- eapply(org.Hs.egGO2ALLEGS, function(x) as.character(na.omit(rids[x])))
go.env <- go.env[unlist(lapply(go.env, length))>5]

# alternative choice of some sepecific go terms
# goterm <- read.table("interested_goterm")  # the intersted_goterm is a file containg interested go term list
# goterm <- t(goterm)
# goterms = goterm[1,]
# go.env <- lapply(mget(goterms, org.Hs.egGO2ALLEGS), function(x) as.character(na.omit(rids[x])))

# add description to GO term
library(GO.db)
desc <- unlist(lapply(mget(names(go.env), GOTERM, ifnotfound = NA), function(x) if(is.logical(x)) { return("") } else { slot(x, "Term")}))
names(go.env) <- paste(names(go.env), desc)  # append description to the names
# convert to an environment
go.env <- clean.gos(go.env)
go.env <- list2env(go.env)

# calculate weighted first principal component magnitudes for each GO gene set in the provided environment
pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = 20, verbose = 1)
save.image("pwpca.RData")


# get singnificant go terms with redundance
df <- pagoda.top.aspects(pwpca, return.table = TRUE, plot = TRUE, z.score = 1.96)
write.table(df,file="tab_sig_aspects",sep="\t",quote=F)

# get significant aspects of heterogeneity
tam <- pagoda.top.aspects(pwpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))
#hc <- pagoda.cluster.cells(tam, varinfo, include.aspects = T)
hc <- pagoda.cluster.cells(tam, varinfo)

pdf("plot_cluster.pdf",width = 9, height = 5)
plot(hc)
dev.off()



# combine pathways that are driven by the same sets of genes
tamr <- pagoda.reduce.loading.redundancy(tam, pwpca)

# combine aspects that show similar patterns
tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot = TRUE, cell.clustering = hc, labRow = NA, labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 0)

# review top aspect
pdf("plot_top_aspect.pdf",width = 9, height = 5)
pagoda.view.aspects(tam, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = rbind(l2cols))
pagoda.view.aspects(tamr, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = rbind(l2cols))
pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20),  col.cols = rbind(l2cols))
dev.off()


# specific Goterms
pdf("GO_0060048_cardiac_muscle_contraction.pdf",width =9, height =5)
d <-pagoda.show.pathways(c("GO:0060048 cardiac muscle contraction"), varinfo, go.env, cell.clustering = hc, margins = c(1,5), show.cell.dendrogram = TRUE, showRowLabels = TRUE, showPC = TRUE, n.genes=200,cexRow=1,return.details=1)
dev.off()
write.table(d$rotation,file="tab_GO_0060048_cardiac_muscle_contraction.pdf",sep="\t",quote=F)
