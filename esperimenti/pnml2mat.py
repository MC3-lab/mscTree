import xml.etree.ElementTree as ET
import sys
import numpy as np

def pnml2gti():
  ns = {'pnml': 'http://www.pnml.org/version-2009/grammar/pnml'}

  tree = ET.parse(sys.argv[1])
  root = tree.getroot()
  #print(root.tag)
  #print(root.attrib)
  #print(root[0][1][0].tag)
  dpl = {}
  dtr = {}
  inm = []

  i = 0
  for place in root[0][0].findall('pnml:place', ns):
      dpl[place.attrib['id']] = i
      if place.find('pnml:initialMarking', ns) is not None:
          marked = place.find('pnml:initialMarking', ns)
          inm.append(int(marked.find('pnml:text', ns).text))
#          inm.append(2)
      else:
          inm.append(0)
      i += 1

  inm = np.array(inm)
  print(inm)

  j = 0    
  for transition in root[0][0].findall('pnml:transition', ns):
      dtr[transition.attrib['id']] = j
      j += 1    

  imat = np.zeros((i,j), dtype=int)
  omat = np.zeros((i,j), dtype=int) 

  for arc in root[0][0].findall('pnml:arc', ns)  :
      if arc.attrib['source'] in dpl.keys():
          imat[dpl[arc.attrib['source']]][dtr[arc.attrib['target']]] = int(arc[0][0].text)
      else:
          omat[dpl[arc.attrib['target']]][dtr[arc.attrib['source']]] = int(arc[0][0].text)

  return imat, omat, inm
