#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

class DYNKdata:
    #Data for one (element,attribute)
    
    element   = None
    attribute = None
    
    turn    = None # Array of turn numbers
    setIDX  = None # Array of setIDX
    funName = None # Array of funname
    data    = None # Matrix of data
    #How many elements are we setting = how many columns in data
    datalen = None 

    def __init__(self, element, attribute,datalen):
        self.element   = element
        self.attribute = attribute

        self.turn    = []
        self.setIDX  = []
        self.funName = []
        self.data    = [] #First index: turn, second, element instance
        
        self.datalen = datalen
        
    def toNumarray(self):
        self.turn = np.asarray(self.turn)
        self.data = np.asarray(self.data)
        assert self.datalen == self.data.shape[1]
    
dynkData = {}

ifile = open("dynksets.dat", 'r')
lnum = 1
for l in ifile:
    l = l.strip()
    if l[0] == "#":
        continue
    ls = l.split()
    
    TURN = int(ls[0])
    EL_ATT = (ls[1],ls[2])
    SETIDX = int(ls[3])
    FUNCTION = ls[4]
    if len(ls) == 6:
        if lnum ==1:
            print "New file format"
        DATALEN = 1
        DATA = [float(ls[5])]
    elif len(ls) >= 7: 
        if lnum ==1:
            print "Old file format"
        DATALEN = int(ls[5])
        DATA = map(float,ls[6:])
    else:
        print "ERROR: Invalid file format"
        exit(1)
    
    if DATALEN != len(DATA):
        print "ERROR! DATALEN != DATA"
        exit(1)

    if not EL_ATT in dynkData:
        dynkData[EL_ATT] = DYNKdata(EL_ATT[0],EL_ATT[1],DATALEN)

    dynkData[EL_ATT].turn.append(TURN)
    dynkData[EL_ATT].setIDX.append(SETIDX)
    dynkData[EL_ATT].funName.append(FUNCTION)
    dynkData[EL_ATT].data.append(DATA)

    lnum = lnum+1

ifile.close()

if len(dynkData) == 0:
    print "No data"
    exit(0)

for (k,i) in zip(dynkData,xrange(len(dynkData))):
    dynkData[k].toNumarray()
    
    #Separate figure:
    plt.figure(10+i)
    plt.plot(dynkData[k].turn,dynkData[k].data[:,0])
    plt.xlabel("Turn")
    plt.ylabel("Setting")
    plt.title("%s:%s" %(k[0],k[1]))
#    plt.grid(True)
#    plt.xticks(np.arange(dynkData[k].turn[0], dynkData[k].turn[-1], 1.0))

    #Combined figure
    plt.figure(1)
    plt.plot(dynkData[k].turn,dynkData[k].data[:,0],
             label="%s:%s"%(k[0],k[1]) )
    plt.xlabel("Turn")
    plt.ylabel("Setting")
    #plt.xticks(np.arange(dynkData[k].turn[0], dynkData[k].turn[-1], 1.0))

#Final touches
plt.figure(1)
plt.legend(loc=0,frameon=False,ncol=2)
plt.grid(True)

plt.show()
