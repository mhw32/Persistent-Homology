source('testlib.r')
for (i in 1:20) {
  foam <- readRDS(paste('./saved_states/test_set/foam', i, '.rds', sep=''))
  base1 <- readRDS(paste('./saved_states/test_set/baseline', i, '-0.1.rds', sep=''))
  test_wrapper(foam, base1, paste(i, '0.1base', sep='-'))
  # Have to reload foam because it gets changed.
  foam <- readRDS(paste('./saved_states/test_set/foam', i, '.rds', sep=''))
  base2 <- readRDS(paste('./saved_states/test_set/baseline', i, '-0.9.rds', sep=''))
  test_wrapper(foam, base2, paste(i, '0.9base', sep='-'))
}