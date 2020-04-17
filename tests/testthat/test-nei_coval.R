context("nei_coval")

set.seed(1)
g <- matrix(sample(c(-1,1),100*1000,replace = TRUE),100,1000)
gmap <- cbind(c(rep(1,nrow(g)/2),rep(2,nrow(g)/2)),c(1:ncol(g)))
x <- runif(nrow(g),1,100)
y <- runif(nrow(g),1,100)
smap <- cbind(x,y)
grouping <- c(rep(1,nrow(g)/2), rep(2,nrow(g)/2))
min_s <- min_dist(smap, grouping)

test_that(
  desc = "nei_coval_range",
  code = {
    g_nei1 <- nei_coval(geno=g, smap=smap, scale=min_s, grouping=grouping)
    g_nei2 <- nei_coval(geno=g, smap=smap, scale=min_s)
    g_nei3 <- nei_coval(geno=g, smap=smap, scale=min_s, grouping=grouping, alpha=min_s/2, kernel="exp")
    g_nei4 <- nei_coval(geno=g, smap=smap, scale=min_s, grouping=grouping, alpha=min_s/2, kernel="gaussian")

    expect_true((min(g_nei1)<=1)&(max(g_nei1)>=-1))
    expect_true((min(g_nei2)<=1)&(max(g_nei2)>=-1))
    expect_true((min(g_nei3)<=1)&(max(g_nei3)>=-1))
    expect_true((min(g_nei4)<=1)&(max(g_nei4)>=-1))
  }
)
