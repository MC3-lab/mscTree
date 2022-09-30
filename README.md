# mscTree

Implementation of the maximal step-computation tree to compute collective reveals relations on Petri nets.
The inputs are a bounded and equal-conflicts Petri net and a collective reveals relation.
The output states if the collective reveals relation is satisfied on the input net.
The net is specified by its two incidence matrices and by its initial marking.
Collective reveals is specified by two lists of transitions, x and y, and by a natural number n. 
The program checks if in all the maximal runs of the net, n occurrences of transitions in x reveals at least 
a transition in y.
It has three possible output values: 
- 0, if n times x collective reveals y
- 1, if there is at least a maximal run with n occurrences of x and none of y
- 2, if there is no maximal run with n occurrences of x.
