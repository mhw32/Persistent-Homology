# Collection of useful tools to share among files.

# Math functions
# Combination
choose <- function(n, k) {
  factorial(n) / (factorial(k)*factorial(n-k))
}
# Permutation
pick <- function(n, k) {
  factorial(n) / factorial(n-k)
}