data(karate,package='igraphdata')
exp_mem<-c(1,4,4,5,5,1,1,4,2,1,2,2,1,2,1,3,3,3,1,3,3,5,1,1,3,1,5,4,4,4,5,1,1,4)
exp_mod<- -0.02638067
test_that("proper class argument is provided", {
  expect_error(rSpectral::spectral_igraph_membership('karate'),'.*igraph.*')
  expect_error(rSpectral::spectral_igraph_communities('karate'),'.*igraph.*')
  expect_named(rSpectral::spectral_igraph_membership(karate),
               c('names','membership'))
  expect_named(rSpectral::spectral_igraph_communities(karate),
               c('vcount', 'algorithm', 'membership', 'modularity', 'names'), 
               ignore.order = TRUE)
})

test_that('membership is correct',{
  c<-rSpectral::spectral_igraph_communities(karate)
  expect_equal(c$membership,exp_mem)
  expect_equal(c$modularity,exp_mod,tolerance=1e-5)
  expect_equal(c$algorithm,'spectral')
})