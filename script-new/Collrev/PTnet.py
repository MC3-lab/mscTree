class PTnet:
  def __init__(self, prem, postm, m0):
    self.prem = prem
    self.postm = postm
    self.m0 = m0
    
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
