library("scde")
library("DBI")

setwd("./")

# input the read count table
count_tab_all <- read.table("summary_star_readcount",head=T,row.names =1)

# drop low express genes and low coverage samples
#cd <- clean.counts(count_tab, min.lib.size = 500, min.reads = 10, min.detected = 5)
#> dim(cd)
#[1] 9021   48

l2cols <- c(rep("coral4",12),  rep("mediumpurple2",9))
sg <- factor(l2cols,levels = c("red","purple"))
names(sg) <- colnames(count_tab)

# building error model
o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 20, threshold.segmentation = TRUE, save.crossfit.plots = T, save.model.plots = T, verbose = 1)

# discard bad samples
#valid.cells <- o.ifm$corr.a > 0
#o.ifm <- o.ifm[valid.cells, ]

# time consuming step, better to save the result
save.image("o.ifm.RData")

# estimate prior 
prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = T)

# call DE genes
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  sg, n.randomizations  =  100, n.cores  =  20, verbose  =  1)

# write to a table
write.table(ediff[order(abs(ediff$Z), decreasing = TRUE), ], file = "tab_DE_genes.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

# get expression magnitude
o.fpm <- scde.expression.magnitude(o.ifm, counts = cd)
write.table(o.fpm, file = "tab_genes_magnitude.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)


