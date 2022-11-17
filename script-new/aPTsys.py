class PTsys:
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
    
  def conflict_set(self, i):
    """Input: an index of a row of the incidence matrix with preconditions
       Output: a set with all the transitions having the row precondition as input
    """
    pre = self.prem[i]
    confl = set([])
    for i in range (0, pre.shape[0]):
      if pre[i] != 0:
        confl.add(i)
    return confl
  
  def diff_mrk(self):
    return self.postm - self.prem
  
  def conflict_partition(self):
    """Input: incidence matrix with preconditions
       Output: list of sets, where each set is a clique of conflict relation
       Since the net is equal-conflict, the sets form a partition of transitions.
    """
    m1_rows, m1_col = self.prem.shape
    added = set([])
    cp = []
    for t in range (0, m1_col):  
      if t not in added:
        i = 0
        while self.prem[i][t] == 0 :
          i = i + 1
        t_confl = self.conflict_set(i)
        added = added.union(t_confl)
        cp.append(t_confl)
    return cp
   
  def check_free_choice(self):
    row, col = self.prem.shape
    tcheck = list(range(0, col))
    for i in range (0, row):
      conflict = []
      for j in tcheck:
        if self.prem[i][j] > 0:
          conflict.append(j)
      if len(conflict) > 1: 
        test = conflict.pop()
        for t in conflict:
          if not all(self.prem[:, test] == self.prem[:, t]):
            return False
          tcheck.remove(t)
        tcheck.remove(test)
    return True
     


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

from mcollrev import *

def prepara():
  im, om, mi = pnml2gti('kanban22.pnml', {}, {})
  return PTsys(im, om, mi)

# Function msmg
# Compute the marking graph by maximal steps of net
# Returns a list of transitions
def msmg (net):
  pmarks  = []                    # List of processed markings
  pending = set()                 # Set of pending markings
  t       = []                    # List of transitions
  adj     = {}

  test    = net.conflict_partition()
  dstep   = net.diff_mrk()
  pending.add(tuple(net.m0))

  while len(pending) > 0:
    current = np.asarray(pending.pop())
    adj[tuple(current)] = []
    ltr     = compute_maximal_steps(net.prem, current, test)
    for step in ltr:
      newm = nextmrk_step(current,step,dstep)
      if not(tuple(newm) in pmarks) and not(np.array_equal(newm, current)):
        pending.add(tuple(newm))
      t.append((current,step,newm))
      adj[tuple(current)].append((tuple(step),tuple(newm)))
    pmarks.append(tuple(current))
    
  return t, pmarks, adj
  
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
  ltr = adj[node.mrk]
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
  
if __name__ == "__main__":
  net = prepara()
  inTrace = tuple(np.zeros((net.postm.shape[1]), int))
  node = Nodo(tuple(net.m0), inTrace)
  t, pmarks, adj = msmg(net)
  leaves = genTree(adj, node, [], [])
  print("N. foglie: ", len(leaves))
  comp_1rev(net, ([0,2], [0], 2), leaves, node)
