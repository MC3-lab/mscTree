import numpy as np
import itertools

class PTnet:
  def __init__(self, prem, postm, m0):
    self.prem = prem
    self.postm = postm
    self.m0 = m0
    
  def eventList(self):
    events = []
    for i in range(0, self.prem.shape[1]):
      pre = self.prem[:,i]
      post = self.postm[:,i]
      ev = Event(pre, post)
      events.append(ev)
    return events

# Class Event: transition of a PT nets; self-loops allowed
class Event:
    def __init__ (self, pre, post, lab = 'a') :
        self.lab = lab
        self.pre = pre
        self.post = post
        
    def enabled (self, marking) :
        i = 0
        flag = 1
        while i < len(marking) and flag:
            if marking[i] < self.pre[i]:
                flag = 0
            i = i + 1
        return flag

# Class Nodo: a node in the max step tree of a free-choice net
class Nodo:
  def __init__(self, marking, trace):
    self.mrk      = marking
    self.children = []
    self.isomrk   = None  # Ancestor with same marking
    self.trace    = trace
    self.dead     = 0  # 1: is a deadlock

  def printsubtree (self, level):
    i = 0
    while i < level:
        print ("    ", end='')
        i = i+1
#
    print(self.mrk, "  ", self.dead, "  ", self.trace)     
    for c in self.children:
      c.printsubtree(level+1)
  
  def __repr__(self):
      return str(self.mrk) + " " + str(self.trace) + "  " + str(self.dead)
      
def conflict_set(pre):
  """Input: a row of the incidence matrix with preconditions
     Output: a set with all the transitions having the row precondition as input
  """
  confl = set([])
  for i in range (0, pre.shape[0]):
    if pre[i] != 0:
      confl.add(i)
  return confl

def conflict_partition(m1):
  """Input: incidence matrix with preconditions
     Output: list of sets, where each set is a clique of conflict relation
     Since the net is equal-conflict, the sets form a partition of transitions.
  """
  m1_rows, m1_col = m1.shape
  added = set([])
  cp = []
  for t in range (0, m1_col):  
    if t not in added:
      i = 0
      while m1[i][t] == 0 :
        i = i + 1
      t_confl = conflict_set(m1[i])
      added = added.union(t_confl)
      cp.append(t_confl)
  return cp

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
        
  
def compute_maximal_steps(m1, marking, confl):
  """Input: incidence matrix with preconditions, a marking, the events partitioned in conflicts
     Output: the maximal steps enabled in the marking
  """
  steps = []
  for c in confl:
    t = c.pop()
    c.add(t)
    n = find_cardinality(t, m1[:,t], marking)
    if n > 0:
      step_part = list(itertools.combinations_with_replacement(c, n))
      steps.append(step_part)
  maxsteps = list(itertools.product(*steps))
  ms = maxstep2list(maxsteps, m1.shape[1])
  return ms

# Function nextmrk: compute new marking when step s fires at m
def nextmrk (m, s):
  nm = m.copy()
  for i in range(len(m)):
    nm[i] = m[i] - s.pre[i] + s.post[i]
  return nm

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

def listenabled(m, events):
  lse = []
  for t in events:
    if t.enabled(m):
      lse.append(t)
  return lse

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
#    ltr = listenabled(node.mrk, events)
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
  j = 0
  while j < len(y):
    if path[y[j]] > 0: 
      return True
    else:
      j = j + 1
  return False
  
def update_paths(paths, x, y):
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
        

#input_m = np.array(([1,1,0,0,0,0,0],[0,0,1,1,0,0,0],[0,0,0,0,1,0,0],[0,0,0,0,0,1,0],[0,0,0,0,1,0,0],[0,0,0,0,0,1,0],[0,0,0,0,0,0,1]))  
#input_m = np.array(([0,0,0,0,0,0],[1,0,0,0,0,0],[0,1,1,0,0,0],[0,0,0,1,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]))
#output_m = np.array(([0,0,0,0,1,0], [0,0,0,1,0,0], [0,0,0,1,0,0], [1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,1],[0,0,1,0,0,0]))
#marking = np.array([0,1,1,0,0,0,0])
#net = PTnet(input_m, output_m, marking)
#marking2 = np.array([2,0,0,1,0,0,0])
#test = conflict_partition(input_m)
#eventi = net.eventList()
#ms = compute_maximal_steps(input_m, marking, test)
#print(ms)
#print(ms[0][0])

#m = [0,1,1,0,0,0,0]

inPNSE = np.array(([1,0,1,0,0,0,0,0], [0,1,0,0,1,0,0,0],[0,0,0,1,0,0,0,0], [0,0,0,0,0,0,0,1], [0,0,0,0,0,1,1,0],[0,0,0,0,0,1,1,0],[0,0,0,0,0,0,0,0]))
outPNSE = np.array(([0,1,0,0,0,0,0,0], [1,0,0,1,0,0,0,0],[0,0,1,0,0,0,0,0], [0,0,0,0,1,0,1,0], [0,0,0,0,1,0,0,1],[0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0]))
m = np.array(([1,0,0,0,0,2,0]))
net = PTnet(inPNSE,outPNSE,m)
test = conflict_partition(net.prem)
eventi = net.eventList()

#a=Event('a',[0,1,0,0,0,0,0],[0,0,0,1,0,0,0])
#b=Event('b',[0,0,1,0,0,0,0],[0,0,0,0,1,0,0])
#c=Event('c',[0,0,1,0,0,0,0],[0,0,0,0,0,1,1])
#d=Event('d',[0,0,0,1,1,0,0],[0,1,1,0,0,0,0])
#e=Event('e',[0,0,0,0,0,1,0],[1,0,0,0,0,0,0])
#f=Event('f',[0,0,0,0,0,0,1],[0,0,0,0,0,1,0])
#eventi = [a,b,c,d,e,f]
trace = np.zeros(len(eventi), np.uint8)
n = Nodo(net.m0, trace)
#enab = listenabled(m, eventi)
vn = []
genMSCT(net.prem, n, vn, [], eventi, test)
n.printsubtree(0)
x = [0,2]
y = [1]
leaves = findPaths(n, x)
for l in leaves:
  print("ancestor", l.isomrk)
ans = collective_reveals(leaves, x, y, 2)
print(ans)
