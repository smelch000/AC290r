import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
#import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['legend.fontsize'] = 10

nodes={}  
order=["W","M1","M2","M4","M8","M16","M32","M64","M128","M256"]
def load(fname,spacing=0,append=False):
    if not append: nodes.clear();
    ps=open(fname,'r').readlines()
        
    if fname[-4:].upper()==".DAT":
        """
        if spacing != 0:
            k="M"+str(spacing)
        else:
            i=5
            while(True):
                try:
                    s=int(fname[-i:-4])
                except:
                    break
                k=str(s)
                i=i+1
                          
        """
        for p in ps:
            nodes[p.split()[3]]=1
        for k in nodes:
            x=[]
            y=[]
            z=[]
            for p in ps:
                if p.split()[3] == k:
                    x.append(int(p.split()[0]))
                    y.append(int(p.split()[1]))
                    z.append(int(p.split()[2]))
            nodes[k]=(list(x),list(y),list(z))
    elif fname[-4:].upper()==".XYZ":
        ps=ps[2:]

        for p in ps:
            nodes[p.split()[0]]=1
        for k in nodes:
            x=[]
            y=[]
            z=[]
            for p in ps:
                if p.split()[0] == k:
                    x.append(int(p.split()[1]))
                    y.append(int(p.split()[2]))
                    z.append(int(p.split()[3]))
            nodes[k]=(list(x),list(y),list(z))
    else:
        raise Exception("Input format unknown");

        
def plot3D(dead=False, tipi=0):
    if tipi==0: tipi=nodes.keys()
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    for k in order:
        if k not in tipi: continue
        if (not dead) and k=='X': continue
        x=nodes[k][0]
        y=nodes[k][1]
        z=nodes[k][2]
        ax.plot(x,y,z,'.',label=k)
    ax.legend()
    plt.show()

def slice(X=-1,Y=-1,Z=-1,tipi=0):
    if tipi==0: tipi=nodes.keys()
    for k in order:
        if k not in tipi: continue
        x=nodes[k][0]
        y=nodes[k][1]
        z=nodes[k][2]
        a=[]
        b=[]
        for p in zip(x,y,z):
            if X != -1 and X == p[0]:
                a.append(p[1])
                b.append(p[2])
            elif  Y != -1 and Y == p[1]:
                a.append(p[0])
                b.append(p[2])
            elif  Z != -1 and Z == p[2]:
                a.append(p[0])
                b.append(p[1])
                
        if len(a) > 0: plt.plot(a,b,'.',label=k)
    plt.legend()
    plt.show()
    
