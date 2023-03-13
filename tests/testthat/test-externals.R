test_that("clusterNet works",{

  test_output <- clusterNet(data = readRDS(test_path("testdata", "test-edge_list.rds")), type = "edge_list",
                            metaboliteA = "Metabolite.A", metaboliteB = "Metabolite.B", tau = 0.5, main.seed = 417)

  expect_equal(test_output$summary,readRDS(test_path("testdata", "test-summary.rds")))
  expect_equal(test_output$subnetworks,readRDS(test_path("testdata", "test-subnetworks.rds")))

})
