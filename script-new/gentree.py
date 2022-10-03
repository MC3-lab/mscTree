#!/usr/bin/python3
import numpy as np
import itertools
import sys
import getopt
from PTnet import *


def find_cardinality(t, col_t, marking):
  """Input: a transition t, the column with the preset of t and a marking
     Output: the number of times in which t can fire in m
  """
  i = 0
  r = marking
  while np.all(r >= col_t):
    r = r - col_t
    i = i + 1
  return i
    
def maxstep2list(maxsteps,numt):
  """Input: the list of maximal steps as tuples of tuples and the number of transitions in the net
     Output: the list of maximal steps as a numpy array where in each position there is the number of 
     occurrences of the transition in the step
  """
  list_steps = []
  if maxsteps != [()]:
    for step in maxsteps:
      lstep = np.zeros(numt, np.uint8)
      for part in step:
        for i in part:
          lstep[i] += 1
      list_steps.append(lstep)
  return list_steps
          
    
def compute_maximal_steps(prem, marking, confl):
  """Input: incidence matrix with preconditions, a marking, the events partitioned in conflicts
     Output: the maximal steps enabled in the marking
  """
  steps = []
  for c in confl:
    t = c.pop()
    c.add(t)
    n = find_cardinality(t, prem[:,t], marking)
    if n > 0:
      step_part = list(itertools.combinations_with_replacement(c, n))
      steps.append(step_part)
  maxsteps = list(itertools.product(*steps))
  ms = maxstep2list(maxsteps, prem.shape[1])
  return ms

def nextmrk_step(m, step, events):
    """Input: the current marking, a step of transitions and the list of events of the net
       Output: the marking obtained from the current marking after the execution of the step
    """
    inpu = np.zeros(len(m), np.uint8)
    outpu = np.zeros(len(m), np.uint8)
    for i in range(0, len(step)):
      inpu = inpu + step[i] * events[i].pre
      outpu = outpu + step[i] * events[i].post
    new_m = m - inpu + outpu
    return new_m

def genMSCT(mat, node, vn, ltr, events, confl_set): 
  """Input: the matrix of preconditions, the root of the tree, the list of nodes that 
     have already been visited, the list of maximal steps (?), the list of events in the net, 
     the sets of conflicts 
     Output: tree of maximal-step computations
  """
  # vn: list of visited nodes in a path 
  for x in vn:
    if np.array_equal(node.mrk, x.mrk) == True:   # Repeated marking
      node.isomrk = x
      return
  ltr = compute_maximal_steps(mat, node.mrk, confl_set)
  cltr = ltr.copy()
  if len(ltr) == 0:    # Deadlock
    node.dead = 1
    return
  else:
    for step in cltr:
      vn.append(node)
      cvn = vn.copy()
      nm = nextmrk_step (node.mrk, step, events)
      newn = Nodo (nm, node.trace + step)
      node.children.append(newn)
      genMSCT (mat, newn, cvn, ltr, events, confl_set)
  return

def findPaths(n, x, leaves = []):
  """Input: a root node n, a list of transitions, a list of leaves (used in the recursice call, initially empty)
     Output: list of leaves such that there is at least an element of x in the path bringing to it
  """
  if n.children == []:
    hasX = False
    i = 0
    while i < len(x) and hasX == False:
      if n.trace[x[i]] > 0:
        leaves.append(n)
        hasX = True
      else:
        i = i + 1
  else:
    for child in n.children:
      leaves = findPaths(child, x, leaves)
  return leaves
  
def elongate_path(n, x, y, root, leaves = []):
  """Input: node of the tree ordered with respect to root, list of revealing events, 
     list of revealed events, starting point of the elongations to consider, list of 
     leaves that have already been added
     Output: list of nodes that should be analysed next
  """
  if n.children == []:
    hasX = False
    ancestor = n.isomrk
    trace = np.zeros(len(n.trace), np.uint8)
    for i in x:
      trace[i] = trace[i] + n.trace[i] - root.trace[i]
    for i in y:
      trace[i] = trace[i] + n.trace[i] - root.trace[i]
    nchild = Nodo(n.mrk, trace)
    nchild.isomrk = ancestor
    if n.dead == 1:
      nchild.dead = 1
    else:
      nchild.dead = 0
    if nchild.dead == 0 and np.all(ancestor.trace >= root.trace) and np.any(ancestor.trace != root.trace):      
      leaves.append(nchild)
    else:
      i = 0
      while i < len(x) and hasX == False:
        if n.trace[x[i]] - root.trace[x[i]] > 0:
          leaves.append(nchild)
          hasX = True
        else:
          i = i + 1
  else:
    for child in n.children:
      leaves = elongate_path(child, x, y, root, leaves)
  return leaves
  
def evaluate_path(path, y):
  """Input: trace of a node, list of potentially revealed events
     Output: True, if an element of y is in the trace; False otherwise
  """
  j = 0
  while j < len(y):
    if path[y[j]] > 0: 
      return True
    else:
      j = j + 1
  return False
  
def update_paths(paths, x, y):
  """Input: list of nodes to extend, list of revealing event, list of revealed events
     Output: updated list of nodes that can still be useful to examine
  """
  new_leaves = []
  for node in paths:
    if node.dead == 0:
      ancestor = node.isomrk
      leaves = elongate_path(ancestor, x, y, ancestor)
      for leaf in leaves:
        leaf.trace = leaf.trace + node.trace
      new_leaves.extend(leaves)
  return new_leaves
  
def collective_reveals(leaves, x, y, n):
  """Input: least of nodes to explore for collective reveals, list of revealing transitions,
     list of revealed transitions, number of occurrences of x for the collective reveals
     Output: 0 if n.x collective reveals y, 1 if it does not, 2 if there is no run with 
     n occurrences of x
  """
  empty = True       
  while leaves != []:
    cl = leaves.copy()
    for leaf in cl:
      m = 0
      for i in x:
        m = m + leaf.trace[i]
      if m >= n:
        empty = False
        b = evaluate_path(leaf.trace, y)
        if b == False:
          return 1
        leaves.remove(leaf)    
    leaves = update_paths(leaves, x, y)
  if empty == True:
    return 2
  else:
    return 0
 
    
def help():
  print("--- Usage ---")
  print("The script needs two arguments:")
  print("1. The name of the file with the net; it can be a ndr file")
  print("  exported from Tina or a txt file where the net is already")
  print("  specified by its two incidence matrices and its initial") 
  print("  marking. In the first case, please use option -n, otherwise,")
  print("  use option -t")
  print("2. The collective reveals relation; it can be given in the")
  print("  command line or in a txt file. In the latter case, please")
  print("  use option -f.")
        
if __name__ == "__main__":
  try:
    opts, args = getopt.getopt(sys.argv[1:], "n:t:f:")
  except getopt.GetoptError:    
    help()
    sys.exit()
  if len(opts) == 0:
    print("Use -t if the net is in a txt file, -n otherwise")
    sys.exit()
  elif len(opts) > 2 or len(args) > 1 or (len(opts) + len(args) != 2) :
    help()
    sys.exit()  
  elif opts[0][0] == "-n":
    print("Sorry, this option is not implemented yet :(") 
    sys.exit()
  else: 
    f = open(opts[0][1])    
    data = ""
    for line in f:
        data += line
    f.close()
    try:
        net = eval(data) 
    except (SyntaxError, NameError):
      stderr.write("%s: syntax error" % argv[0])
      if form != stdin:
        stderr.write(" in %s" % argv[2])
      stderr.write("\n")
      exit()
    if len(opts) == 2:
      f = open(opts[1][1])    
      data = ""
      for line in f:
        data += line
      f.close()
      try:
        rev = eval(data) 
      except (SyntaxError, NameError):
        stderr.write("%s: syntax error" % argv[0])
        exit() 
    else:
      rev = eval(args[0])     
#    net = PTnet(inPNSE,outPNSE,m)
    test = net.conflict_partition()
    eventi = net.eventList()
    trace = np.zeros(len(eventi), np.uint8)
    n = Nodo(net.m0, trace)
    genMSCT(net.prem, n, [], [], eventi, test)
    n.printsubtree(0)
    x = rev[0]
    y = rev[1]
    z = rev[2]
    leaves = findPaths(n, x)
    ans = collective_reveals(leaves, x, y, z)
    if ans == 0:
      print("The set "+ str(x) + " "+ str(z) + "-collective reveals" + str(y)) 
    elif ans == 1:
      print("The set "+ str(x) + " does NOT "+ str(z) + "-collective reveals" + str(y)) 
    else:
      print("There are no run with " + str(z) + "occurrences of" + str(x))  
