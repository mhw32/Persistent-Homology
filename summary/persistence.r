# =====================================================================
# Permutation Tests [Too Slow]
# =====================================================================
# 1. NHST test.
# Should check that nhst is working as expected.
# distfunc <- bottleneckDist
# distfunc <- wassersteinDist
# proba <- nhst(X, L, permN, bottleneckDist)
# This works but it is so slow. 
# system.time(bottleneck(X[[1,1]], X[[2,1]], dimension=1))
# 1.335 * 225 * 2 * 1000 (group size 15 --> 15^2, 2 groups, 1000 iterations.)
# Even 100 iterations = 16 hours...
 
# system.time(bottleneck(X[[1,1]], X[[2,1]], dimension=2))
# user  system elapsed
# 0.212  0.024  0.237

# Is the 2-wasserstein faster (dimension 1)? This is so slow.
# Crashes. 160.656 seconds.

# 2-wasserstein faster (dimension 2)
# user  system elapsed
# 5.760   0.014   5.796

# No this is much slower. How do I make these raw functions faster? Otherwise, tests for as big as the Universe will take decades. I don't think this is a viable path.
# ---------------------------------------------------------------------
# In fact, tests 2-4 all require some form of a distance measurement. This is going to be very problematic since everything takes such a long time. Perhaps we must move to landscapes directly. How did you run these? 

# 2. Gaussian kernel (max likelihood h).
# distfunc <- bottleneckDist
# kernel <- gaussianKernel
# If we want to find the best value for parameter h:
# best <- findBestParam(X, L)
# proba <- permutationTest(permN, X, L, kernelStat, distfunc, best)
# ---------------------------------------------------------------------
# 3. Energy kernel.
# distfunc <- bottleneckDist
# proba <- permutationTest(permN, X, L, energyStat, distfunc, 1)
# ---------------------------------------------------------------------
# 4. Rossenbaum kernel.
# No need for permutation tests for still takes the annoying distance func.
# distfunc <- bottleneckDist
# proba <- rosenbaumStat(X, L, distfunc)
 # =====================================================================
