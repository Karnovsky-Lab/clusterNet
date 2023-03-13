test_that("run_consensus_cluster works",{

  expect_equal(run_consensus_cluster(readRDS(test_path("testdata", "test-inputNetwork.rds")),
                                     tau = 0.5, main.seed = 417)$dcl,
               readRDS(test_path("testdata","test-consensus_membership.rds")))
})

test_that("gatherConsensus_Matrix works",{

  expect_equal(gatherConsensusMatrix(readRDS(test_path("testdata","test-consensus_membership.rds"))),
               readRDS(test_path("testdata", "test-subnetwork_results.rds")))
})


