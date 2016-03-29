source('testlib.r')
for (i in 1:100) {
  foam <- readRDS(paste('./saved_states/large_set/foam', i, '.rds', sep=''))
  base <- readRDS(paste('./saved_states/large_set/baseline', i, '-0.1.rds', sep=''))
  print(paste('Iteration: ', i))
  test_wrapper(foam, base, paste(i, '0.1baseNormFalse', sep='-'), FALSE)
  
  # Have to reload foam because it gets changed.
  foam <- readRDS(paste('./saved_states/large_set/foam', i, '.rds', sep=''))
  base2 <- readRDS(paste('./saved_states/large_set/baseline', i, '-0.9.rds', sep=''))
  test_wrapper(foam, base2, paste(i, '0.9baseNormFalse', sep='-'), FALSE)
	
  foam <- readRDS(paste('./saved_states/large_set/foam', i, '.rds', sep=''))
  base1 <- readRDS(paste('./saved_states/large_set/baseline', i, '-0.1.rds', sep=''))
  test_wrapper(foam, base1, paste(i, '0.1baseNormTrue', sep='-'), TRUE)
  
  # Have to reload foam because it gets changed.
  foam <- readRDS(paste('./saved_states/large_set/foam', i, '.rds', sep=''))
  base2 <- readRDS(paste('./saved_states/large_set/baseline', i, '-0.9.rds', sep=''))
  test_wrapper(foam, base2, paste(i, '0.9baseNormTrue', sep='-'), TRUE)

}
