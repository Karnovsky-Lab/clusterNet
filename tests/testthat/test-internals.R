test_that("run_consensus_cluster works",{


  test_output <- run_consensus_cluster(readRDS(test_path("testdata", "test-inputNetwork.rds")),
                                       tau = 0.5, main.seed = 417)

  #tests the subnetwork output
  expect_equal(test_output$dcl, readRDS(test_path("testdata","test-consensus_membership.rds")))
  #really is testing that get_consensus_matrix is working properly
  expect_equal(test_output$D, readRDS(test_path("testdata", "test-D_adjacency_matrix.rds")))
})

test_that("gatherConsensus_Matrix works",{

  expect_equal(gatherConsensusMatrix(readRDS(test_path("testdata","test-consensus_membership.rds"))),
               readRDS(test_path("testdata", "test-subnetwork_results.rds")))
})


