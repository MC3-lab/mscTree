# mscTree

Implementation of the maximal step-computation tree to compute collective reveals relations on Petri nets.
The inputs are a bounded and equal-conflicts Petri net and a collective reveals relation.
The output states if the collective reveals relation is satisfied on the input net.
The net is specified by its two incidence matrices and by its initial marking.
Collective reveals is specified by two lists of transitions, x and y, and by a natural number n. 
The program checks if in all the maximal runs of the net, n occurrences of transitions in x reveals at least 
a transition in y.
It has three possible output values: 
- n times x collective reveals y
- there is at least a maximal run with n occurrences of x and none of y
- there is no maximal run with n occurrences of x.

Requirements:
Python3, numpy, itertools, sys, getopt, numba

Current use:
Run gentree.py with python3 and specify the net and the collective reveals relation as follows:
1. The net can be in a text file represented as PTnet(<input_incidence_matrix>, <output_incidence_matrix>, <initial_marking>)
all these 3 elements must be expressed as numpy array (see examples with txt extensions).
If this is the case, the option to be used is -t.
Alternatively, the program can read a pnml file, with the option -p followed by the name of the file
2. The reveals relation can be both in a txt file or in a string directly given on the command line;
its a triple represented as (<list_revealing_events>, <list_revealed_events>, number_of_occurrences).
(TO FIX: also if the file is passed with pnml, the transitions in the reveals relation need to be
passed with their column number in the incidence matrix)

On the command line, you can write: python3 gentree.py -t <path_file_PTnet> -f <name_file_reveals>
or python3 gentree.py -t <name_file_PTnet> "rev_relations".
There is a third possibility, if you want to produce the tree and test more reveals relations on it, 
consisting in using the option -i after the name of the file with the net. In this case, the program 
computes the tree and starts an interactive mode, therefore you can try as many formulas as you like without
computing everytime the tree. When you are done, tipe "quit()" to exit.
For some reason, after the -i you need to add a string, it is not really used in the program, so everything
you type is ok (it needs to be fixed). In this last case, the command would be something like:
python3 gentree.py -t <path_file_PTnet> -i randomCrap.

Example: python3 gentree.py -t examples/net1.txt -f examples/rev1.txt
or
python3 gentree.py -t examples/net1.txt "([0,2],[1], 2)"
