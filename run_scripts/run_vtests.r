source('testlib.r')
for (i in 86:100) {
  foam <- readRDS(paste('./saved_states/large_sub_set/foam', i, '.rds', sep=''))
  base <- readRDS(paste('./saved_states/large_sub_set/baseline', i, '.rds', sep=''))
  print(paste('Iteration: ', i))
  test_wrapper(foam, base, paste(i, '0.1baseNormFalse', sep='-'), FALSE, FALSE)
	
  foam <- readRDS(paste('./saved_states/large_sub_set/foam', i, '.rds', sep=''))
  base1 <- readRDS(paste('./saved_states/large_sub_set/baseline', i, '.rds', sep=''))
  test_wrapper(foam, base1, paste(i, '0.1baseNormTrue', sep='-'), TRUE)
}
