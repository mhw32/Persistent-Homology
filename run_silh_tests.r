source('testlib.r')
for (i in 1:20) {
  for (p in seq(1,5,0.1)) {
    foam <- readRDS(paste('./saved_states/test_set/foam', i, '.rds', sep=''))
    base <- readRDS(paste('./saved_states/test_set/baseline', i, '-0.1.rds', sep=''))
    test_wrapper(foam, base, paste('iter', i, 'tune', p, 'norm', 'false', sep='-'), p, norm=FALSE)
  }
}