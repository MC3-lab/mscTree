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
