## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----input--------------------------------------------------------------------
library(rNeighborGWAS)

#simulate "fake_nei" dataset using nei_simu()
set.seed(1)
g <- matrix(sample(c(-1,1),100*1000,replace = TRUE),100,1000)
gmap <- cbind(c(rep(1,nrow(g)/2),rep(2,nrow(g)/2)),c(1:ncol(g)))
x <- runif(nrow(g),1,100)
y <- runif(nrow(g),1,100)
smap <- cbind(x,y)
grouping <- c(rep(1,nrow(g)/2), rep(2,nrow(g)/2))
pheno <- nei_simu(geno=g, smap=smap, scale=44,
                  grouping=grouping, n_causal=50,
                  pveB=0.4, pve=0.8
                  )

fake_nei <- list()
fake_nei[[1]] <- g
fake_nei[[2]] <- gmap
fake_nei[[3]] <- smap
fake_nei[[4]] <- data.frame(pheno,grouping)
names(fake_nei) <- c("geno","gmap","smap","pheno")

fake_nei$geno[1:5,1:10]
head(fake_nei$smap)

## ----PVE----------------------------------------------------------------------
min_s <- min_dist(fake_nei$smap, fake_nei$pheno$grouping)
scale_seq <- c(min_s, quantile(dist(fake_nei$smap),c(0.2*rep(1:5))))

pve_out <- calc_PVEnei(geno=fake_nei$geno, pheno=fake_nei$pheno[,1],
                       smap=fake_nei$smap, scale_seq=scale_seq,
                       addcovar=as.matrix(fake_nei$pheno$grouping),
                       grouping=fake_nei$pheno$grouping
                       )
delta_PVE(pve_out)

## ----GWAS---------------------------------------------------------------------
scale <- 43
gwas_out <- neiGWAS(geno=fake_nei$geno, pheno=fake_nei$pheno[,1],
                    gmap=fake_nei$gmap, smap=fake_nei$smap,
                    scale=scale, addcovar=as.matrix(fake_nei$pheno$grouping),
                    grouping=fake_nei$pheno$grouping
                    )

gaston::manhattan(gwas_out)
gaston::qqplot.pvalues(gwas_out$p)

