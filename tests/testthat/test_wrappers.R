test_data = generate_test_data()
niter=50
burnin=0
thinin=5


# TEST FOR COMPUTING FUNCTIONS

test_that("rhm_horseshoe works",{
  chain_horseshoe = rhm_horseshoe(y=test_data$y,V_list=test_data$V_list,D_list=test_data$D_list,niter=niter,burnin=burnin,thinin=thinin)
  expect_equal(length(chain_horseshoe$sigma2),(niter-burnin)%/%thinin)
})

test_that("rhm_fusion works",{
  chain_fusion = rhm_fusion(y=test_data$y,V_list=test_data$V_list,D_list=test_data$D_list,niter=niter,burnin=burnin,thinin=thinin)
  expect_equal(length(chain_fusion$sigma2),(niter-burnin)%/%thinin)
})

test_that("rhm_fused works",{
  chain_fused = rhm_fused(y=test_data$y,V_list=test_data$V_list,D_list=test_data$D_list,niter=niter,burnin=burnin,thinin=thinin)
  expect_equal(length(chain_fused$sigma2),(niter-burnin)%/%thinin)
})

test_that("fast_rhm works",{
  fast_result = fast_rhm(y=test_data$y,V_list=test_data$V_list,D_list=test_data$D_list)
  expect_equal(length(fast_result),2)
})


# TEST FOR PREPARATION FUNCTIONS

snp_data = generate_snp_data()
interval = 20

test_that("cut_snp_matrices works", {
  X_list = cut_snp_matrices(snp_data$chr_list,snp_data$position_list,interval)
  expect_equal(length(X_list),length(snp_data$chr_list))
})

X_list = cut_snp_matrices(snp_data$chr_list,snp_data$position,interval)

test_that("get_similarity_matrices works", {
  G_list = get_snp_similarity_matrices(X_list)
  expect_equal(length(G_list),length(snp_data$chr_list))
})

# TEST FOR FINAL WRAPPER

test_that("rhm works", {
  dat = generate_snp_data(by_chr=FALSE)
  result = rhm(dat$y,dat$X,dat$chromosome_list,dat$position_list,block_size=30,niter=niter,thinin=thinin,burnin=burnin)
  expect_equal(length(result$sigma2),(niter-burnin)%/%thinin)
})

