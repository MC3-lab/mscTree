import numpy as np
import re
import sys

np.set_printoptions(threshold=sys.maxsize)

places = 0
transitions = 0
marking = []
lines = 0

with open("buf100.ll_net") as f:
  for line in f:
    lines += 1
    if bool(re.search('^"P\d+"', line)) == True:
      places += 1
      if bool(re.search('m1M1$', line)) == True:
        marking.append(1)
      else:
        marking.append(0)
    elif bool(re.search('^"T\d+"', line)) == True:
      transitions += 1
    elif bool(re.search("^TP", line)) == True:
      input_net = np.zeros((places,transitions), dtype = int)
      output_net = np.zeros((places,transitions), dtype = int)
    elif bool(re.search('\d<\d', line)) == True:
      x = line.split("<")
      output_net[int(x[1])-1][int(x[0])-1] = 1
    elif bool(re.search('\d>\d', line)) == True:
      x = line.split(">")
      input_net[int(x[0])-1][int(x[1])-1] = 1

marking = np.array(marking)      
print("PTnet(")
print("np." +repr(input_net)+ ",") 
print("np."+repr(output_net)+",")
print("np."+repr(marking))
print(")")
