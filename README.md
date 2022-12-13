# mscTree

Implementation of the maximal step-computation tree to compute collective reveals relations on Petri nets.
The inputs are a bounded and equal-conflicts Petri net and a collective reveals relation.
The output states if the collective reveals relation is satisfied on the input net.
The net must be provided in a PNML file; alternatively, a user can directly give a PTnet object,
made of two incidence matrices and an initial marking, all in the form of numpy arrays.
Collective reveals is specified by a triple *(x, y, n)*: the first element is the list *x* of all the revealing transitions,
the second element is the list *y* of revealed transition, and the third element is a natural number *n*.
The transitions must be specified with their id in the PNML file.
The program checks if in all the maximal runs of the net, *n* occurrences of transitions in *x* reveals at least 
a transition in *y*.
It has two possible output values: 
- true, if the relation holds or there is no maximal run with *n* occurrences of *x*; 
- false, if the relation does not hold

Requirements:
Python3, numpy, itertools, sys, getopt

How to use:
Run collrev.py with python3 and specify the net and the collective reveals relation as follows:
1. If the net is in a PNML file, use the option -p followed by the name of the file.
Alternatively, the net can be in a text file represented as PTnet(<input_incidence_matrix>, <output_incidence_matrix>, <initial_marking>).
All these three elements must be expressed as numpy array. 
[//]: # (see examples with txt extensions).
If this is the case, the option to be used is -t.
2. The reveals relation can be either in a file or in a string given as argument on the command line;
in the first case, use the option -f, otherwise no option is needed.
In both cases the relation must be a triple represented as (<list_of_revealing_events>, <list_of_revealed_events>, number_of_occurrences).
3. The program can also be executed in interactive mode; in this case, the user needs to specify a PT system followed by the
option -i. The tool constructs a structure on which all the collective reveals relations can be computed, and asks
which collective reveals relation needs to be checked. Users can check any number of collective reveals relations, and type quit to
end the session. This mode is generally slower than the one in which the formula is given at the beginning, so it is advised to
use only if the user needs to check several relations.
4. If the users need to check a large number or all the reveals relations on the PT system, they can specify the name of the net as described above, 
followed by the option -a. 

Examples:
- python3 collrev.py -p examples/net1.pnml -f examples/rev1.txt
- python3 collrev.py -p examples/net1.pnml "(['id1', 'id2'], ['id3', 'id4'], 5)"
- python3 collrev.py -p examples/net1.pnml -i
- python3 collrev.py -p examples/net1.pnml -a
