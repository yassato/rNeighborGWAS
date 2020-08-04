context("nei_lmm")

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
g_nei <- nei_coval(geno=fake_nei$geno, smap=fake_nei$smap, scale=scale, grouping=fake_nei$pheno$grouping)

test_that(
  desc = "neiLMM_neiGWAS_equal",
  code = {
    gwas_lmm <- neiGWAS(geno=fake_nei$geno, pheno=fake_nei$pheno[,1],
                        gmap=fake_nei$gmap, smap=fake_nei$smap,
                        scale=scale, addcovar=as.matrix(fake_nei$pheno$grouping),
                        grouping=fake_nei$pheno$grouping)

    res_lmm <- nei_lmm(geno=fake_nei$geno,g_nei=g_nei,
                       pheno=fake_nei$pheno$pheno,
                       addcovar=as.matrix(fake_nei$pheno$grouping))

    gwas_glmm <- neiGWAS(geno=fake_nei$geno, pheno=pheno_bin,
                         gmap=fake_nei$gmap, smap=fake_nei$smap,
                         scale=scale, addcovar=as.matrix(fake_nei$pheno$grouping),
                         grouping=fake_nei$pheno$grouping, response="binary")

    res_glmm <- nei_lmm(geno=fake_nei$geno,g_nei=g_nei,
                        pheno=pheno_bin,
                        addcovar=as.matrix(fake_nei$pheno$grouping), response="binary")

    expect_equal(-log10(gwas_lmm$p), -log10(res_lmm$p_nei))
    expect_equal(-log10(gwas_glmm$p), -log10(res_glmm$p_nei))
  }
)
