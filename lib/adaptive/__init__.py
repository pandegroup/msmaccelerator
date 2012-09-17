# Adaptive __init__.py

### Adaptive Sampling Methods -------------------------------------------------
#
# This module contains the general classes used to construct adaptive sampling
# routines. While there is a lot of flexibility in writing or adding to these
# classes, we have adopted the following convention that should be very general.
#
# Each class below is designed to produce some set of "weights" for each state,
# in general these weights can be though of as a discrete categorial distribution
# across the MSM space. Then, each state will be sampled with probability
# proportional to its weight. The weights are set by different methods to
# optimize some MSM parameter of interest, whether total states explored, 
# uncertainty in a specific eigenmode, etc.
#
# These routines have been standardized to be used in the following two-step
# fashion:
# 
# (1) Initialize:       sampler = SamplingMethod( args_controlling_algorithm )
# (2) Request weights:  multinomial = sampler.sample( transition_counts_matrix )
#
# The big benefit of this method is that it allows for simple ways to combine 
# methods when one may want to optimize many features at once. A simple linear
# combination of weights yields a new categorical distribution that contains
# information from all included methods, allowing the workflow to proceed in
# a smooth and simple fashion.
#

