import numpy as np
import re
import sys

np.set_printoptions(threshold=sys.maxsize)

places = 0
transitions = 0
marking = []
lines = 0

with open("buf100.ll_net") as f:
    s = f.readline()
    while s != "PL\n":   
      s = f.readline()
    s = f.readline()
    while s != "TR\n":
      places += 1
      if bool(re.search('m1M1$', s)) == True:
        marking.append(1)
      else:
        marking.append(0)
      s = f.readline()
    s = f.readline()  
    while s != "TP\n":
      transitions += 1
      s = f.readline()
    input_net = np.zeros((places,transitions), dtype = int)
    output_net = np.zeros((places,transitions), dtype = int)
    s = f.readline()
    while s != "PT\n": 
      x = s.split("<")
      output_net[int(x[1])-1][int(x[0])-1] = 1
      s = f.readline()
    s = f.readline()
    while bool(re.search('^TX', s)) == False:
      x = s.split(">")
      input_net[int(x[0])-1][int(x[1])-1] = 1
      s = f.readline()

marking = np.array(marking)      
print("PTnet(")
print("np." +repr(input_net)+ ",") 
print("np."+repr(output_net)+",")
print("np."+repr(marking))
print(")")
