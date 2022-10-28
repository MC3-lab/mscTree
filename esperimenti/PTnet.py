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
