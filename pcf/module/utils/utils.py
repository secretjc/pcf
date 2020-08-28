#!/usr/bin/env python

"""

utils.py
"""

####
#### Imports
####
from pyomo.opt import SolverFactory

####
#### Authorship information
####
__author__ = "Yiyang Chang"
__copyright__ = "Copyright 2016, Purdue ISL Robustness Framerwork Project"
__credits__ = ["Yiyang Chang"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Yiyang Chang"
__email__ = "chang256@purdue.edu"
__status__ = "alpha"


class Utils(object):
  @classmethod
  def get_num_of_zeros_or_ones(cls, prefix):
    return sum(1 for c in prefix if c == '1' or c == '0')

  @classmethod
  def get_num_of_unkonwns(cls, prefix):
    return sum(1 for num in prefix.split('|') if num == '?')

  @classmethod
  def get_num_of_ones(cls, prefix):
    return sum(int(num) for num in prefix.split('|')
               if num != '?' and num != '0')

  @classmethod
  def get_pyomo_solver(cls, solver_name, option_d):
    opt = SolverFactory(solver_name, tee=True)
    for k, v in option_d.iteritems():
      opt.options[k] = v
    return opt

  @classmethod
  def unique_permutations(cls, seq):
    """
    Yield only unique permutations of seq in an efficient way.

    A python implementation of Knuth's "Algorithm L", also known from the
    std::next_permutation function of C++, and as the permutation algorithm
    of Narayana Pandita.
    """

    # Precalculate the indices we'll be iterating over for speed
    i_indices = range(len(seq) - 1, -1, -1)
    k_indices = i_indices[1:]

    # The algorithm specifies to start with a sorted version
    seq = sorted(seq)

    while True:
      yield seq

      # Working backwards from the last-but-one index,           k
      # we find the index of the first decrease in value.  0 0 1 0 1 1 1 0
      for k in k_indices:
        if seq[k] < seq[k + 1]:
          break
      else:
        # Introducing the slightly unknown python for-else syntax:
        # else is executed only if the break statement was never reached.
        # If this is the case, seq is weakly decreasing, and we're done.
        return

      # Get item from sequence only once, for speed
      k_val = seq[k]

      # Working backwards starting with the last item,           k     i
      # find the first one greater than the one at k       0 0 1 0 1 1 1 0
      for i in i_indices:
        if k_val < seq[i]:
          break

      # Swap them in the most efficient way
      (seq[k], seq[i]) = (seq[i], seq[k])                #       k     i
                                                         # 0 0 1 1 1 1 0 0

      # Reverse the part after but not                           k
      # including k, also efficiently.                     0 0 1 1 0 0 1 1
      seq[k + 1:] = seq[-1:k:-1]
