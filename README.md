# quickPDI
Efficient calculation of Van Calster's Polytomous Descrimination Index (PDI) when there are many categories
# Let a "case" be the combination of y-value and corresponding probabilities.
# Together is now a data frame that has, for each unique value of case:
# the case and freq, where freq is the repeat-count of z.  That is, col(1) is
# y, cols 2:ncategories+1 are the probs, and ncategories+2 is freq.
#
# For a given column j, Van Calster inspects all sets of k rows such that
# every one of the k row indices is represented (once) in the set.  Each
# such set is assigned a score that is 1 if the probability from the row
# with index j is the largest probability in column j. Otherwise, the
# score is 0. The sum of these scores, when divided by prod(N), is the PDI
# for column j. The overall PDI is the mean of the column-specific PDI's.
#
# The problem is that there are prod(N) such sets.  If every set must be
# individualy inspected, this number can be prohibitively large,
# especially if there are many categories.  For example, if there are 100
# records from each of 3 categories, there will be 1 million sets.  A method
# is needed to efficiently count them.
#
# One can begin to see a solution upon realizing that PDI is a rank-based
# statistic.  Suppose for a moment that there are no ties.  Consider a set
# with score 1. It has this score because the probability from the row
# with index j is the largest probability in column j. Call this
# probability X. If the column is sorted in ascending order, all the
# members of this set that are not from column j will fall below X.
# Suppose (z1, ..., zk) counts the number of row indices 1, ..., k that
# fall below X in the sorted column.  Then the number of sets with scores
# of 1 will be Prod(z') where z' is z excluding element j. Summing these
# over all X such that row index(X)=j, and dividing by prod(N), will give
# the PDI for column j.
#
# Van Calster handles ties as follows:  For each set of cases, if the
# probability (in column j) from the row with index j is largest, and
# letting t be the number of probabilities in column j equal to the
# largest (so that the number of probabilites tied with the largest is
# t-1), the score should be 1/t.
#
# In contrst, Li scores as 1 only sets in which the row which index j is
# undisputedly larges.  Sets in which there are ties get scored as 0.
#
# How to count the Van Calster possibilities efficiently?  After the sort, all tied
# probabilities will be grouped together.  Considering only this tied
# group, let (w1, ... wk) be the count of records with index 1, ..., k within
# this tied group.  The next part of the calculation is best described by
# example.  Suppose that there are 4 categories, and we are just now
# concentrating on the 3rd of these, and that, just after the w have been
# calculated, we have
#
#    w = a b c d
#    z = A B C D
#
# Then, for each of the c cases that are counted in w(3), the number of
# groups in which there are no ties (i.e. t = 1) is
#
# s(1) = ABD
#
# Also for each of the c cases, the number of groups in which there is
# only 1 tie (i.e. t = 2) is:
#
# s(2) = aBD + AbD + ABd
#
# Continuing:
#
# s(3) = abD + aBd + Abd
# s(4) = abd
#
# Note that s(1)+s(2)+s(3)+s(4) = prod(w+z)[-j].  We now have that the
# number of groups in which there are ties (i.e. cases containing
# probabilities in column j that are equal to the probability of the index
# case in column j) is s(2)+s(3)+s(4) = w[j]*(prod(w+z)[-j]-prod(z[-j])
#
# (A more efficient way to calculate the s-values will be shown below.)
#
# This means that the increment to the score that comes from this
# collection of sets is:  c * [1*s1 + 1/2*s2 + ...  + 1/k*s(k)].  This
# increment needs to be added to the appropriate element of SUMS, the
# cumlative total of such increments
#
# FOr Li's method, the score is c*s1.
#
# Following this, update z = z + w, and continue the calculation with the
# next row, just after the end of this tied group.
#
# A more efficient way to calculate the s-values follows:
#
# Define f(t) = (ta + A)(tb + B)(td + D)
#
# Clearly, f(0) = s(1) and f(1) = s(1) + s(2) + s(3) + s(4).  In fact, it
# turns out that:
#
# [1 0 0  0     [s(1)      [f(0)
#  1 1 1  1   *  s(2)   =   f(1)
#  1 2 4  8      s(3)       f(2)
#  1 3 9 27]     s(4)]      f(3)]
#
# If there are k categories, M, the matrix on the left, is k-by-k with row
# i+1 contaning the elements [i^0, i^1, ..., i^(k-1)], for i = 1, ..., k-1.
# That is, M(i+1,j+1) = i^j, for , for i,j = 1, ..., k-1.
#
# M is the transpose of a nonsingular Vandermonde matrix.  This means that, given the
# easily-calculated f-values, one can calculate s=M^(-1)*f.  If t` is the
# row-vector (1 1/2 1/3 1/4), the entire score increment for the tied
# group is c*t`*M^(-1)*f.  But given that we are interested in t`*M^(-1),
# instead of calculating M^(-1) directly, it is faster and more accurate
# to solve for q in M'q = t, and then calculate the score increment as
# c*q`*f. q need be calculated only once, then used later with multiple c
# and f-values.
#
# Noet that, in Li's case, t` = (1 0 0 0).
#
# Suppose that in some group there are no ties for probabilities with row
# index j, that is, a = b = d = 0, and all c probabilities in the group
# have row indexes j.  In this case, f(0)=f(1)=f(2)=f(3)=ABD.  This means
# that the matrix expression above becomes
#
# [1 0 0  0     [s(1)            [1
#  1 1 1  1   *  s(2)   =  ABD *  1
#  1 2 4  8      s(3)             1
#  1 3 9 27]     s(4)]            1]
#
# Because S(1) = ABD, it must be the case that
#
#   [1 1  1     [s(2)       [0
#    2 4  8   *  s(3)    =   0
#    3 9 27]     s(4)]       0]
#
# from which, by elementary row operations, s(2)=s(3)=s(4)=0.
# This means that
#
#
#                [f(0)                   [1      [ABD
#      M^(-1) *   f(1)  = ABD *  M^(-1) * 1   =   0
#                 f(2)                    1       0
#                 f(3)]                   1]      0]
#
# so that the increment to the score is c * ABD, which is what it would
# have been in the original method of counting that ignored ties.  So the
# method that adjusts for ties reduces to the non-tied method when there
# are no ties.
#
# References:
#
# Van Calster B, Vergouwe Y, Looman CWN, Van Belle V, Timmerman D and
# Steyerberg EW.  Assessing the discriminative ability of risk models for
# more than two outcome categories.  European Journal of Epidemiology
# 2012; 27:  761 C 770.
#
# Li, J., Feng, Q., Fine, J.P., Pencina, M.J., Van Calster, B. (2017).
# Nonparametric estimation and inference for polytomous discrimination
# index.  Statistical Methods in Medical Research.  In Press.
