Sys.setenv(OMP_NUM_THREADS = "1", OMP_THREAD_LIMIT = "1",
           OMP_WAIT_POLICY = "PASSIVE")

library(testthat)
library(MIML)

try(lightgbm::setLGBMthreads(1L), silent = TRUE)

test_check("MIML")
