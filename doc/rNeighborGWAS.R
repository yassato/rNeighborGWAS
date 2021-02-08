## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----input--------------------------------------------------------------------
set.seed(1234)
library(rNeighborGWAS)

# convert "TTN" genotype data into a rNeighborGWAS format
data("TTN", package="gaston")
x <- gaston::as.bed.matrix(TTN.gen, TTN.fam, TTN.bim)
g <- gaston2neiGWAS(x)

# simulate "fake_nei" dataset using nei_simu()
geno <- g$geno
gmap <- g$gmap
x <- runif(nrow(geno),1,100)
y <- runif(nrow(geno),1,100)
smap <- cbind(x,y)
grouping <- c(rep(1,nrow(geno)/2), rep(2,nrow(geno)/2), 2)
pheno <- nei_simu(geno=geno, smap=smap, scale=43,
                  grouping=grouping, n_causal=50,
                  pveB=0.3, pve=0.6
                  )

fake_nei <- list()
fake_nei[[1]] <- geno
fake_nei[[2]] <- gmap
fake_nei[[3]] <- smap
fake_nei[[4]] <- data.frame(pheno,grouping)
names(fake_nei) <- c("geno","gmap","smap","pheno")

fake_nei$geno[1:5,1:10] # Note: 0 indicates heterozygotes
head(fake_nei$smap)

## ----PVE----------------------------------------------------------------------
scale_seq <- quantile(dist(fake_nei$smap),c(0.2*rep(1:5)))

pve_out <- calc_PVEnei(geno=fake_nei$geno, pheno=fake_nei$pheno[,1],
                       smap=fake_nei$smap, scale_seq=scale_seq,
                       addcovar=as.matrix(fake_nei$pheno$grouping),
                       grouping=fake_nei$pheno$grouping
                       )
delta_PVE(pve_out)

## ----GWAS---------------------------------------------------------------------
scale <- 43.9
gwas_out <- neiGWAS(geno=fake_nei$geno, pheno=fake_nei$pheno[,1],
                    gmap=fake_nei$gmap, smap=fake_nei$smap,
                    scale=scale, addcovar=as.matrix(fake_nei$pheno$grouping),
                    grouping=fake_nei$pheno$grouping
                    )

gaston::manhattan(gwas_out)
gaston::qqplot.pvalues(gwas_out$p)

## ----LMM, eval=FALSE----------------------------------------------------------
#  scale <- 43.9
#  g_nei <- nei_coval(geno=fake_nei$geno, smap=fake_nei$smap,
#                     scale=scale, grouping=fake_nei$pheno$grouping
#                     )
#  
#  gwas_out <- nei_lmm(geno=fake_nei$geno, g_nei=g_nei,
#                      pheno=fake_nei$pheno[,1],
#                      addcovar=as.matrix(fake_nei$pheno$grouping)
#                      )

## ----bin, eval=FALSE----------------------------------------------------------
#  fake_nei$pheno[,1][fake_nei$pheno[,1]>mean(fake_nei$pheno[,1])] <- 1
#  fake_nei$pheno[,1][fake_nei$pheno[,1]!=1] <- 0
#  
#  pve_out <- calc_PVEnei(geno=fake_nei$geno, pheno=fake_nei$pheno[,1],
#                         smap=fake_nei$smap, scale_seq=scale_seq,
#                         addcovar=as.matrix(fake_nei$pheno$grouping),
#                         grouping=fake_nei$pheno$grouping,
#                         response="binary"
#                         )
#  
#  gwas_out <- neiGWAS(geno=fake_nei$geno, pheno=fake_nei$pheno[,1],
#                      gmap=fake_nei$gmap, smap=fake_nei$smap,
#                      scale=scale, addcovar=as.matrix(fake_nei$pheno$grouping),
#                      grouping=fake_nei$pheno$grouping,
#                      response="binary"
#                      )
#  gaston::manhattan(gwas_out)
#  gaston::qqplot.pvalues(gwas_out$p)
#  
#  gwas_out <- nei_lmm(geno=fake_nei$geno, g_nei=g_nei,
#                      pheno=fake_nei$pheno[,1],
#                      addcovar=as.matrix(fake_nei$pheno$grouping),
#                      response="binary"
#                      )

## ----asymmetry, eval=FALSE----------------------------------------------------
#  scale <- 43.9
#  g_nei <- nei_coval(geno=fake_nei$geno, smap=fake_nei$smap,
#                     scale=scale, grouping=fake_nei$pheno$grouping
#                     )
#  
#  gwas_out <- nei_lmm(geno=fake_nei$geno, g_nei=g_nei,
#                      pheno=fake_nei$pheno[,1],
#                      addcovar=as.matrix(fake_nei$pheno$grouping),
#                      asym=TRUE)

