#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import os

fileDType = np.dtype([('ID', np.int), ('turn', np.int),
                      ('s', np.float),('x', np.float),('y', np.float),('z', np.float),
                      ('xp', np.float),('yp', np.float),('dEE', np.float),
                      ('pType', np.int)])

def readdumpfile(fname):
    fdata = np.loadtxt(fname, dtype=fileDType)
    
    turnIdxs = [[fdata['turn'][0]]] #1st idx: Pass (initially 64 particles)
                                    #2nd idx: Turn number (last element = end)
    
    tPrev   = fdata['turn'][0]
    passIdx = 0

    for i in xrange(1,len(fdata['turn'])):
        t = fdata['turn'][i]
        if t < tPrev:
            turnIdxs[passIdx].append(i)
            passIdx += 1
            turnIdxs.append([])
        if t == tPrev:
            continue
        turnIdxs[passIdx].append(i)
        tPrev = t
    turnIdxs[-1].append(len(fdata['turn']))
    
    return (fdata,turnIdxs)

def plot_particleNum(fdata,turnIdx):
    plt.figure()

    turn = []
    numParticles = []
    
    for i in xrange(len(turnIdxs)):
        if (len(turnIdxs[i])-1) > len(turn):
            for j in xrange(len(turn),len(turnIdxs[i])-1):
                turn.append(j+1)
                numParticles.append(0)
        for j in xrange(len(turnIdxs[i])-1):
            #assumption: One element/turn
            numParticles[j] += turnIdxs[i][j+1]-turnIdxs[i][j]

    plt.plot(turn,numParticles)
    plt.xlabel("Turn")
    plt.ylabel("Particles tracked")
    
def get_turndata(fdata, turn):
    #"Turn" input variable counts from 1
    assert turn > 0

    rdata = []

    for i in xrange(len(turnIdxs)):
        if len(turnIdxs[i])-1 < turn:
            continue
        
        tidx1 = turnIdxs[i][turn-1]
        tidx2 = turnIdxs[i][turn]
        for j in xrange(tidx1,tidx2):
            rdata.append(fdata[j])
    return rdata
    
(fdata,turnIdxs) = readdumpfile("../fort.660") #IP1
#(fdata,turnIdxs) = readdumpfile("../fort.661") #CRAB5

plot_particleNum(fdata,turnIdxs)

corr = []

for t in xrange(1,61):
    print "TURN=",t
    tdata = np.asarray(get_turndata(fdata,t),dtype=fileDType)

    corr.append(np.corrcoef(tdata['z'], tdata['x'])[0,1])    

    
    plt.figure(5)
    
    #zx
    #plt.figure()
    plt.clf()
    plt.title("TURN =" + str(t))
    plt.plot(tdata['z'], tdata['x'],'+')
    plt.xlabel("z [mm]")
    plt.ylabel("x [mm]")
    plt.xlim(-200,200)
    plt.ylim(-0.05,0.05)

    plt.savefig("pngs/zx_%05i.png" % (t))

    #plt.figure()
    plt.clf()
    plt.title("TURN =" + str(t))
    plt.plot(tdata['x'], tdata['y'],'+')
    plt.xlabel("x [mm]")
    plt.ylabel("y [mm]")
    plt.xlim(-0.05,0.05)
    plt.ylim(-0.25,0.25)
    
    plt.savefig("pngs/xy_%05i.png" % (t))

    continue

    fig = plt.figure(10)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(tdata['x'], tdata['y'], tdata['z'])
    ax.set_xlim(-0.05,0.05)
    ax.set_ylim(-0.25,0.25)
    ax.set_zlim(-200,200)
    ax.set_title("TURN =" + str(t))
    plt.xlabel("x [mm]")
    plt.ylabel("y [mm]")
    ax.set_zlabel("z [mm]")
    
    plt.savefig("pngs/xyz_%05i.png" % (t))

fps = 10

movieFileName = "pngs/xy.mp4"
print "Assembling movie '" + movieFileName + "' at", fps, "fps:"
os.system("rm " + movieFileName)
command = "ffmpeg -sameq -r "+str(fps)+" -i " + "pngs/xy_%05d.png " + movieFileName
print "Command = '" + command + "'"
os.system(command)

print "Converting to .gif:"
movieFileName = "pngs/xy.gif"
command = "convert " + "pngs/xy_*.png -layers Optimize -delay " + str(100/fps) + " " + movieFileName
print "Command = '" + command + "'"
os.system(command)


movieFileName = "pngs/zx.mp4"
print "Assembling movie '" + movieFileName + "' at", fps, "fps:"
os.system("rm " + movieFileName)
command = "ffmpeg -sameq -r "+str(fps)+" -i " + "pngs/zx_%05d.png " + movieFileName
print "Command = '" + command + "'"
os.system(command)

print "Converting to .gif:"
movieFileName = "pngs/zx.gif"
command = "convert " + "pngs/zx_*.png -layers Optimize -delay " + str(100/fps) + " " + movieFileName
print "Command = '" + command + "'"
os.system(command)

plt.figure(5)
plt.clf()
plt.plot(corr)
plt.xlabel("Turn")
plt.ylabel("Corr")

plt.show()
