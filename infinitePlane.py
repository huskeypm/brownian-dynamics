#
# Special class for 2D finite plane
# absorbing/reflect BC problem
#
from bdIterator import *
from finitePlane import *
import numpy as np
import sympy as sp


class infinitePlane(finitePlane):             
  def __init__(self,
    boundaryReflect,
    q0,
    nParticles,
    z=0,
    name="traj.png",
    title="Traj",
    box=True
    ):

    ### THIS POART MIGHT BE REDUNTANT THROUGH INHERITANCE 

    self.boundaryAbsorb = 0.
    self.q0 = q0 # origin (x) for particles 
    self.boundaryReflect = boundaryReflect # for the x-coordinate only. We're assuming particle is absorbed/reflected on edges
    self.nParticles=nParticles
    self.absorbed = np.zeros(nParticles)
    self.active   = np.ones(nParticles)
    self.inactive   = np.zeros(nParticles)#,1))
    #absorbedvsT=np.zeros(int(Dt/dt))
    self.lftIdx = []
    self.absCtr=0
    self.name=name
    self.title=title

    #self.forceScale = 1.
    # back boundary is at x=-1
    # abs boundary is at x=1
    # y is periodic 
    #boxMin = np.array([-1,-1])
    #boxMax = np.array([ 1, 1])
    self.boxMin = np.array([self.boundaryAbsorb,self.boundaryAbsorb])
    self.boxMax = np.array([boundaryReflect,boundaryReflect])
    self.boxD = self.boxMax - self.boxMin

    # params.obstacle 
    self.obsMin = params.obsMin                
    self.obsMax = params.obsMax                
    if(box==False):
        self.obsMin = np.array([5,5])
        self.obsMax = np.array([5.0001,5.0001])
    elif(box==True): 
        1
    else:
        spl  = box.split("x") # 6x6
        dims = np.array([np.float(spl[0]),np.float(spl[1])]) # 6x6
        print "Using box with ", dims
        ctr=np.array([5,5]) 
        self.obsMin = ctr-dims/2.       
        self.obsMax = ctr+dims/2.                   



    # params 
    self.kT = 0.59     # energy [kcal/mol]
    self.z = z         # unit charge 


    self.times =[]
    self.absorbedvsT=[]

    ## Initialize everything
    self.init(nParticles)

    self.rgtPeriod = np.zeros(nParticles)
    self.lftPeriod = np.zeros(nParticles)
    self.btmPeriod = np.zeros(nParticles)
    self.topPeriod = np.zeros(nParticles)

  def calc(self,x,xN,taui):

     ## obstacle
     # either reflect or redraw. Easiest to reflect
     # top: (top-x)*(x-bottom) > 0 if within bounds 
     # left: (right-x)*(x-left) > 0 if within bounds 
     #xN[0,:] = [params.obsMax[0]-0.001,0]
     #xN[1,:] = [0, params.obsMax[1]-0.001]
     #xN[2,:] = [params.obsMax[1]-0.001, params.obsMax[1]-0.001]
     inX = (self.obsMax[0]-xN[:,0])*(xN[:,0]-self.obsMin[0]) >= 0
     inY = (self.obsMax[1]-xN[:,1])*(xN[:,1]-self.obsMin[1]) >= 0
     inBox = np.where(inY*inX)
     #print "Collision", xN[inBox,]  
     # replace 'collision with original
     xO = xN.copy()
     xN[inBox] = x[inBox]

     ## store dx before periodic BC are applied
     dxp = xN - x




     ## check bounds 
     rgtPeriodicIdx = np.where(xN[:,0] > self.boxMax[0]) # xN[0,:] = [-1.25,0] 
     lftPeriodicIdx  = np.where(xN[:,0] < self.boxMin[0]) # xN[1,:] = [1.25,0] 
     btmPeriodIdx = np.where(xN[:,1] < self.boxMin[1]) # xN[2,:] = [0,1.25]
     topPeriodIdx = np.where(xN[:,1] > self.boxMax[1]) # xN[3,:] = [0,-1.25]
     self.rgtPeriod[rgtPeriodicIdx]+=1
     self.lftPeriod[lftPeriodicIdx]+=1
     self.btmPeriod[btmPeriodIdx]+=1
     self.topPeriod[topPeriodIdx]+=1
     #print self.rgtPeriod
     #print self.lftPeriod
     #print self.btmPeriod
     #print self.topPeriod



     # periodic (validates)
     xL = xN[lftPeriodicIdx,0] - self.boxMin[0]
     xN[lftPeriodicIdx,0] = xL + self.boxMax[0]
     xR = xN[rgtPeriodicIdx,0] - self.boxMax[0]
     xN[rgtPeriodicIdx,0] = xR + self.boxMin[0]
     xB = xN[btmPeriodIdx,1] - self.boxMin[1]
     xN[btmPeriodIdx,1] = xB + self.boxMax[1]
     xT = xN[topPeriodIdx,1] - self.boxMax[1]
     xN[topPeriodIdx,1] = xT + self.boxMin[1]
     #print "B/T Period", xN

     ## apply PBC
     unwrap =np.array([self.rgtPeriod*self.boxD[0] - self.lftPeriod*self.boxD[0],self.topPeriod*self.boxD[1] - self.btmPeriod*self.boxD[1]])
     #print "unwra", unwrap
     xp = xN.copy(); xp[:,0] += unwrap[0]; xp[:,1] += unwrap[1] # + unwrap
     #print "x",x
     #print "xN",xN
     #print "xp",xp
     #print "SDS"

     ## polot 
     self.plottinghack(xO,xN,inBox)
     self.times.append(taui)
     return (xN,xp,dxp)

  #only prints to screen interactively from within ipython 
  def plottinghack(self,xO,xN,inBox):
    # can quess time step by length of arrays (stupid trick)(
    i = np.shape(np.asarray(self.times))[0]


    mod = 10
    if(np.mod(i,mod)!=0):
      return

    print i

    fig = params.fig
    fig.clear()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal', 'box')
    title = self.title+" %d"%(i)
    plt.title(title)
    #plt.xlim([-1,1]) #[params.boxMin[0],params.boxMax[0]])
    #plt.ylim([-1,1]) #[params.boxMin[1],params.boxMax[1]])
    plt.xlim([self.boxMin[0],self.boxMax[0]])
    plt.ylim([self.boxMin[1],self.boxMax[1]])
    plt.plot(self.obs[:,0],self.obs[:,1],'k-')
    #print self.obs
    #plt.plot(xN[outBox,0],xN[outBox,1],'b.')
    plt.plot(xN[:,0],xN[:,1],'b.')
    showCollisions=1
    if(showCollisions==1 and np.shape(inBox)[0] > 0):
      plt.plot(xO[inBox,0],xO[inBox,1],'r*')
      plt.plot(xN[inBox,0],xN[inBox,1],'r.')
      #fig.canvas.draw()
      #plt.draw()
      #raw_input("Enter key") 

    #plt.gca().set_aspect('equal', 'box')
    fig.canvas.draw()
    plt.draw()


    #plt.show()
    #plt.gcf().savefig(self.name)

  
