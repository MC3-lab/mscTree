#!/usr/bin/python3
import numpy as np
import itertools
import sys
import getopt
import xml.etree.ElementTree as ET
import time
from numba import njit, prange
from PTnet import *
from Nodo import *

# Needed to print the matrix of reveals relations with option -a
np.set_printoptions(threshold=sys.maxsize)

def pnml2gti(fileName, dpl, dtr):
  ns = {'pnml': 'http://www.pnml.org/version-2009/grammar/pnml'}

  tree = ET.parse(fileName)
  root = tree.getroot()
  inm = []

  i = 0
  places = tree.findall(".//{http://www.pnml.org/version-2009/grammar/pnml}place")
  for p in places:
      dpl[p.attrib['id']] = i
      if p.find('pnml:initialMarking', ns) is not None:
          marked = p.find('pnml:initialMarking', ns)
          inm.append(int(marked.find('pnml:text', ns).text))
      else:
          inm.append(0)
      i += 1

  inm = np.array(inm)

  transitions = tree.findall(".//{http://www.pnml.org/version-2009/grammar/pnml}transition")
  j = 0    
  for t in transitions:
      dtr[t.attrib['id']] = j
      j += 1    

  imat = np.zeros((i,j), dtype=int)
  omat = np.zeros((i,j), dtype=int) 

  arcs = root.findall(".//{http://www.pnml.org/version-2009/grammar/pnml}arc")
  for arc in arcs:
      isweighted = arc.findall("./{http://www.pnml.org/version-2009/grammar/pnml}inscription")
      if isweighted:
          weight = int(arc[0][0].text)
      else:
          weight = 1

      if arc.attrib['source'] in dpl.keys():
          # source is a place
          imat[dpl[arc.attrib['source']]][dtr[arc.attrib['target']]] = weight
      else:
          # source is a transition
          omat[dpl[arc.attrib['target']]][dtr[arc.attrib['source']]] = weight

  return imat, omat, inm
  
def pnml2nf(pf, dtr):
  qrn = ([], [], pf[2])

  for t in pf[0]:
    qrn[0].append(dtr[t])
  for t in pf[1]:
    qrn[1].append(dtr[t])
  return qrn

def rmv(l, m):
  for nxt in l:
    if nxt[1] == m:
      l.remove(nxt)
  return l

@njit
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
  
def nextmrk_step(m, step, dstep):
    """Input: the current marking, a step of transitions and the list of events of the net
       Output: the marking obtained from the current marking after the execution of the step
    """
    change = np.dot(dstep, step)
    new_m = m + change
    return new_m
  
def update_mg(adj, invAdj, bd):
  if bd == []:
    return adj
  new_bad = []
  for m in bd:
    pb = invAdj[m]
    for prec in pb:
      adj[prec] = rmv(adj[prec], m)
      if adj[prec] == []:
        new_bad.append(prec)
        adj.pop(prec)
  adj = update_mg(adj, invAdj, new_bad) 
  return adj
  
def redmsmg (net, revd):
  pmarks  = []                    # List of processed markings
  pending = set()                 # Set of pending markings
  adj     = {}
  invAdj  = {}
  bd      = []

  test    = net.conflict_partition()
  dstep   = net.diff_mrk()
  pending.add(tuple(net.m0))

  while len(pending) > 0:
    current = np.asarray(pending.pop())
    adj[tuple(current)] = []
    ltr     = compute_maximal_steps(net.prem, current, test)
    if ltr == []:
      new_dead = False
    else:
      new_dead = True
      for step in ltr:
        flag = 0
        for z in revd:
          if step[z] > 0:
            flag = 1
        if flag == 0:
          new_dead = False
          newm = nextmrk_step(current,step,dstep)
          if not(tuple(newm) in pmarks) and not(np.array_equal(newm, current)):
            pending.add(tuple(newm))
          adj[tuple(current)].append((tuple(step),tuple(newm)))
          if tuple(newm) in invAdj.keys():
            invAdj[tuple(newm)].append(tuple(current))
          else:
            invAdj[tuple(newm)] = [tuple(current)] 
    if new_dead == True:
      bd.append(tuple(current))
    pmarks.append(tuple(current)) 
  # End while
  return update_mg(adj, invAdj, bd)
  
def genTree(adj, node, vn, leaves = []): 
  """Input: the matrix of preconditions, the root of the tree, the list of nodes that 
     have already been visited, the list of maximal steps (?), the list of events in the net, 
     the sets of conflicts 
     Output: tree of maximal-step computations
  """
  # vn: list of visited nodes in a path
  for x in vn:
    if np.array_equal(node.mrk, x.mrk) == True:   # Repeated marking
      node.isomrk = x
      leaves.append(node)
      return leaves
  ltr = adj[tuple(node.mrk)]
  if len(ltr) == 0:    # Deadlock
    node.dead = 1
    leaves.append(node)
    return leaves
  else:
    for step in ltr:
      vn.append(node)
      cvn = vn.copy()
      ntr = list(node.trace)
      for i in range(len(ntr)):
          ntr[i] += step[0][i]
      newn = Nodo (step[1], tuple(ntr))
      node.children.append(newn)
      genTree (adj, newn, cvn, leaves)
  return leaves
  
@njit(parallel = True)
def update_rev(rev, trace, live1):
  for i in prange(len(trace)):
    if trace[i] > 0:
      live1[i] = 1
      for j in prange(len(trace)):
        if trace[j] == 0:
          rev[i][j] = 1 
  return rev, live1

def comp_leaves(x, leaves):
  xleaves = []
  for i in prange(len(leaves)):
    app = False
    j = 0
    while j < len(x) and app == False: 
      if leaves[i].trace[x[j]] > 0:
        xleaves.append(leaves[i])
        app = True
      j += 1
  return xleaves
  
def elongate_path(n, x, y, root, cleaves = []):
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
      cleaves.append(nchild)
    else:
      i = 0
      while i < len(x) and hasX == False:
        if n.trace[x[i]] - root.trace[x[i]] > 0:
          cleaves.append(nchild)
          hasX = True
        else:
          i = i + 1
  else:
    for child in n.children:
      cleaves = elongate_path(child, x, y, root, cleaves)
  return cleaves
  
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
      cleaves = elongate_path(ancestor, x, y, ancestor)
      for leaf in cleaves:
        leaf.trace = leaf.trace + node.trace
      new_leaves.extend(cleaves)
  return new_leaves

def collective_reveals(ileaves, x, y, n):
  """Input: least of nodes to explore for collective reveals, list of revealing transitions,
     list of revealed transitions, number of occurrences of x for the collective reveals
     Output: 0 if n.x collective reveals y, 1 if it does not, 2 if there is no run with 
     n occurrences of x
  """
  cleaves = ileaves.copy()      
  while cleaves != []:
    cl = cleaves.copy()
    for leaf in cl:
      m = 0
      for i in x:
        m = m + leaf.trace[i]
      if m >= n:
        b = evaluate_path(leaf.trace, y)
        if b == False:
          return 1
        cleaves.remove(leaf) 
    if cleaves != []:
      cleaves = update_paths(cleaves, x, y)
  return 0
  
def comp_1rev(net, r, leaves, n):
    x = np.array(r[0])
    y = np.array(r[1])
    z = r[2]
    xleaves = comp_leaves(x, leaves)
    return collective_reveals(xleaves, x, y, z)


if __name__ == "__main__":
  fn = None
  ft = None
  fr = None
  interactive = False
  al = False
  try:
    opts, args = getopt.getopt(sys.argv[1:], "p:t:f:ia")
  except getopt.GetoptError as e:
    print(e)
    sys.exit(2)
  for o,a in opts:
    if   o == "-p":
      fn = a
    elif o == "-t":
      ft = a
    elif o == "-f":
      fr = a
    elif o == "-i":
      interactive = True
    elif o == "-a":
      al = True
      
  if fn: # Input from PNML file
    dpl = {}
    dtr = {}
    im, om, m0 = pnml2gti(fn,dpl,dtr)    
    net = PTnet(im, om, m0)
  elif ft: 
    with open(ft, 'r') as f:
      data = ""
      for line in f:
        data += line
      try:
        net = eval(data) 
      except (SyntaxError, NameError):
        stderr.write("%s: syntax error" % argv[0])
        if form != stdin:
          stderr.write(" in %s" % argv[2])
          stderr.write("\n")
          sys.exit(2)
  else:
      print("Either use -p or -t to specify a PT net")
      sys.exit(2)
      
  if not(args or fr or interactive or al):
      print("Please, either specify a collective reveals relation to verify,")
      print("or use -i for the interactive mode, or -a to compute all the")
      print("reveals relations")
      sys.exit(2)
      
  fc = net.check_free_choice()
  if fc == False:
      print("This net is not equal-conflict")
      sys.exit(2)      
  else:
      if args:
        qr = eval(args[0])
        if fn:
          qr = pnml2nf(qr, dtr) #reveal question
      elif fr:
        with open(fr, 'r') as f:
          for line in f:
            qr = eval(line)
        if fn:
          qr = pnml2nf(qr, dtr) 
      else:  # interactive execution or all reveals execution
        qr = ([], [], 1)
      start = time.time()
      adj = redmsmg(net, qr[1])
      end = time.time()
      print("Time for the marking graph: " + str(end - start))
      trace = np.zeros(net.prem.shape[1], np.uint8)
      n = Nodo(net.m0, trace)
      leaves = genTree(adj, n, [], [])
      if al == True:
        nTr = net.prem.shape[1]
        rev = np.zeros([nTr,nTr], np.uint8)
        live1 = np.zeros(nTr, np.uint8)
        for leaf in leaves:
          rev, live1 = update_rev(rev, leaf.trace, live1)
        print("Matrix of reveals relations:")
        print(rev)
        print("-----------------")
        print("1-live transitions: " + str(live1))  
        if fn:
          print(dtr)
      if args or fr:
        ans = comp_1rev(net, qr, leaves, n)
        if ans == 0:
          print("The set "+ str(qr[0]) + " "+ str(qr[2]) + "-collective reveals " + str(qr[1])) 
        else:
          print("The set "+ str(qr[0]) + " does NOT "+ str(qr[2]) + "-collective reveals " + str(qr[1])) 
