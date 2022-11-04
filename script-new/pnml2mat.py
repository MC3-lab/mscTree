import xml.etree.ElementTree as ET
import sys
import numpy as np

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
  print(inm)

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
