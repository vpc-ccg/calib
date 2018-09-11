# Simulated datasets tests

This directory includes the following:

- A CSV file for each dataset tested in Calib's paper with the following columns:
  - Comment: dataset name and tool name
  - Timestamp: When this test record finished running
  - User (CPU) time
  - Wall clock (real) time
  - Max RAM use
  - Adjusted Rand Index (ARI) accuracy
  - Parameter list and their values

## Running simulated datasets tests

Please check the testing script available [here](../../slurm_scripts/).