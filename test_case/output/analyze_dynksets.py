#! /usr/bin/env python

import numpy as np

import matplotlib.pyplot as plt

#data structure:
# - one key = one setidx
# - one index = one turn
turn    = {}
element = {}
attribute = {}
function = {}
active = {}
 # turn, element, attribute, function, active
data = {}     # value 1, value2, ...

ifile = open("dynksets.dat", 'r')
for l in ifile:
    l = l.strip()
    if l[0] == "#":
        continue
    ls = l.split()
    
    setidx = int(ls[1])

    if not setidx in turn:
        turn     [setidx] = []
        element  [setidx] = []
        attribute[setidx] = []
        function [setidx] = []
        active   [setidx] = []
        data     [setidx] = []
    # if len(turn[setidx]) > 0:
    #     if int(ls[0]) < turn[setidx][-1]:
    #         #break
    #         pass
        
    turn     [setidx].append(  int(ls[0]) )
    element  [setidx].append(      ls[2]  )
    attribute[setidx].append(      ls[3]  )
    function [setidx].append(      ls[4]  )
    active   [setidx].append( bool(ls[5]) )
    
    numdata = int(ls[6])
    assert len(ls) == numdata+7
    
    data[setidx].append( map(float,ls[7:]) )
    
ifile.close()

for k in data:
    turn[k] = np.asarray(turn[k])
    data[k] = np.asarray(data[k])

for k in data:
    print k
    # print turn[k]
    # print data[k][:,0]
    # print (k,element[k][0],attribute[k][0],function[k][0])
    plt.figure()
    plt.title("SET #%i acting on %s:%s by function '%s'" %(k,element[k][0],attribute[k][0],function[k][0]))
    plt.xlabel("turn")
    plt.plot(turn[k],data[k][:,0])

plt.show()
