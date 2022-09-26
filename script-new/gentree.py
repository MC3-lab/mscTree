# Class Event: transition of a PT nets; self-loops allowed
class Event:
    def __init__ (self, lab, pre, post) :
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
  def __init__(self, marking):
    self.mrk      = marking
    self.children = []
    self.isomrk   = None  # Ancestor with same marking
    self.label    = None
    self.dead     = 0  # 1: is a deadlock

  def printsubtree (self, level):
    i = 0
    while i < level:
        print ("    ", end='')
        i = i+1
    print (self.label, " ", self.mrk, "  ", self.dead)
    for c in self.children:
      c.printsubtree(level+1)

# Function nextmrk: compute new marking when step s fires at m
def nextmrk (m, s):
  nm = m.copy()
  for i in range(len(m)):
    nm[i] = m[i] - s.pre[i] + s.post[i]
  return nm

def listenabled(m, events):
  lse = []
  for t in events:
    if t.enabled(m):
      lse.append(t)
  return lse

def genMSCT (node, vm, ltr, events):
  if node.mrk in vm:   # Repeated marking
    #update ancestor of node
    return
  else:
    ltr = listenabled(node.mrk, events)
    cltr = ltr.copy()
    if len(ltr) == 0:    # Deadlock
      node.dead = 1
      return
    else:
      for t in cltr:
        vm.append(node.mrk)
        cvm = vm.copy()
        nm = nextmrk (node.mrk, t)
        newn = Nodo (nm)
        newn.label = t.lab
        node.children.append(newn)
        genMSCT (newn, cvm, ltr, events)
    return

m = [0,1,1,0,0,0,0]

a=Event('a',[0,1,0,0,0,0,0],[0,0,0,1,0,0,0])
b=Event('b',[0,0,1,0,0,0,0],[0,0,0,0,1,0,0])
c=Event('c',[0,0,1,0,0,0,0],[0,0,0,0,0,1,1])
d=Event('d',[0,0,0,1,1,0,0],[0,1,1,0,0,0,0])
e=Event('e',[0,0,0,0,0,1,0],[1,0,0,0,0,0,0])
f=Event('f',[0,0,0,0,0,0,1],[0,0,0,0,0,1,0])
eventi = [a,b,c,d,e,f]
n = Nodo(m)
enab = listenabled(m, eventi)
vmrk = []
genMSCT(n, vmrk, [], eventi)
