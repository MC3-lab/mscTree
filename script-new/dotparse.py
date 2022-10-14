# Reads all the lines of a dot file in a list
import re
import numpy as np

dotfile = open('blob.dot', 'r')
righe   = dotfile.readlines()

ptarcs  = []
tparcs  = []
cp      = {} # Condition -> place
et      = {} # Event -> transition
pld     = {}
trd     = {}

nplaces = 0
ntrans  = 0

ptedge  = re.compile('c([0-9]*) -> e([0-9]*)')
tpedge  = re.compile('e([0-9]*) -> c([0-9]*)')
plabel  = re.compile('c([0-9]*) \S*P([0-9]*)')
tlabel  = re.compile('e([0-9]*) \S*T([0-9]*)')

for r in righe:
  pt = ptedge.search(r)
  if pt:
    ptarcs.append((int(pt.groups()[0]),int(pt.groups()[1])))
  else:
    tp = tpedge.search(r)
    if tp:
      tparcs.append((int(tp.groups()[0]),int(tp.groups()[1])))
    else:
      pl = plabel.search(r)
      if pl:
        cp[int(pl.groups()[0])] = int(pl.groups()[1])
        if int(pl.groups()[1]) not in pld.keys():
          pld[int(pl.groups()[1])] = nplaces
          nplaces += 1
      else:
        tl = tlabel.search(r)
        if tl:
          et[int(tl.groups()[0])] = int(tl.groups()[1])
          if int(tl.groups()[1]) not in trd.keys():
            trd[int(tl.groups()[1])] = ntrans
            ntrans += 1

imat = np.zeros((nplaces,ntrans), dtype=int)
omat = np.zeros((nplaces,ntrans), dtype=int)

for arco in ptarcs:
  imat[pld[cp[arco[0]]]][trd[et[arco[1]]]] = 1
for arco in tparcs:
  omat[pld[cp[arco[1]]]][trd[et[arco[0]]]] = 1

for i in range(0,nplaces):
  print(imat[i,:])

print('-------------------------')

for i in range(0,nplaces):
  print(omat[i,:])
