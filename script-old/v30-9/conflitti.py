import numpy as np
import itertools

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
  
def conflict_set(pre):
  """Input: a row of the incidence matrix with preconditions
     Output: a set with all the transitions having the row precondition as input
  """
  confl = set([])
  for i in range (0, pre.shape[0]):
    if pre[i] != 0:
      confl.add(i)
  return confl
  
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
  return maxsteps
  
   
m1 = np.array(([1,1,0,0,0,0,0],[0,0,1,1,0,0,0],[0,0,0,0,1,0,0],[0,0,0,0,0,1,0],[0,0,0,0,1,0,0],[0,0,0,0,0,1,0],[0,0,0,0,0,0,1]))  
test = conflict_partition(m1)
marking = np.array([1,1,0,0,0,0,0])
ms = compute_maximal_steps(m1, marking, test)

print(ms)  
