context("neiGWAS")

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

scale <- 43
pheno_bin <- as.numeric(pheno>mean(pheno))

test_that(
  desc = "nei_GWAS_works",
  code = {

    gwas_lmm <- neiGWAS(geno=fake_nei$geno, pheno=fake_nei$pheno[,1],
                        gmap=fake_nei$gmap, smap=fake_nei$smap,
                        scale=scale, addcovar=as.matrix(fake_nei$pheno$grouping),
                        grouping=fake_nei$pheno$grouping)

    gwas_glmm <- neiGWAS(geno=fake_nei$geno, pheno=pheno_bin,
                        gmap=fake_nei$gmap, smap=fake_nei$smap,
                        scale=scale, addcovar=as.matrix(fake_nei$pheno$grouping),
                        grouping=fake_nei$pheno$grouping, response="binary")

    gwas_lm <- neiGWAS(geno=fake_nei$geno, pheno=fake_nei$pheno[,1],
                        gmap=fake_nei$gmap, smap=fake_nei$smap,
                        scale=scale, addcovar=as.matrix(fake_nei$pheno$grouping),
                        grouping=fake_nei$pheno$grouping, model="lm")

    gwas_glm <- neiGWAS(geno=fake_nei$geno, pheno=pheno_bin,
                         gmap=fake_nei$gmap, smap=fake_nei$smap,
                         scale=scale, addcovar=as.matrix(fake_nei$pheno$grouping),
                         grouping=fake_nei$pheno$grouping, response="binary", model="lm")

    expect_true(min(-log10(gwas_lmm[,3]))>=0)

    expect_true(min(-log10(gwas_glmm[,3]))>=0)

    expect_true(min(-log10(gwas_lm[,3]))>=0)

    expect_true(min(-log10(gwas_glm[,3]))>=0)
  }
)
