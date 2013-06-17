#
# See TODOs
#

import matplotlib
#matplotlib.use('Agg')

import sympy as sp
import numpy as np 
import matplotlib.pylab as plt 


def gaussian(u0=np.random.rand(1),u1=np.random.rand(1)):
    #    // Box-Muller transform
    #    // uniform() returns a uniformly distributed number strictly between 0 and 1
    pi2 = 2*np.pi
    #     double u0 = uniform();
    #    double u1 = uniform();
    lterm = np.sqrt(-2*np.log( u0));
    return lterm*np.cos( pi2*u1);


class params:
  fig = plt.figure()
  name = "test2d"

def drawplot1d(x,xN,i):
  mod = 10
  if(np.mod(i,1)!=0):
    return

  fig = params.fig
  fig.clear()

  ax = fig.add_subplot(111)
  #ax.set_aspect('equal', 'box')
  plt.title("1D sim")
  plt.xlim([-1,10]) #[params.boxMin[0],params.boxMax[0]])
  plt.scatter(x ,np.arange(x.shape[0]),c="b")
  plt.scatter(xN,np.arange(xN.shape[0]),c="r")
  plt.ylabel("Particle")
  plt.xlabel("x")
  plt.draw()

def drawplot2d(x,xN,i):
  mod = 10
  if(np.mod(i,1)!=0):
    return
  fig = params.fig
  fig.clear()
  ax = fig.add_subplot(121)
  ax.set_aspect('equal', 'box')
  plt.title(params.name)
  plt.xlim([0,15]) #[params.boxMin[0],params.boxMax[0]])
  plt.ylim([0,15]) #[params.boxMin[1],params.boxMax[1]])
  #ADD#plt.plot(params.obs[:,0],params.obs[:,1],'k-')
  #plt.plot(xN[outBox,0],xN[outBox,1],'b.')
  plt.plot(xN[:,0],xN[:,1],'b.')
  #ADDshowCollisions=0
  #ADDif(showCollisions==1 and np.shape(inBox)[0] > 0):
  #ADD  plt.plot(xO[inBox,0],xO[inBox,1],'r*')
  #ADD  plt.plot(xN[inBox,0],xN[inBox,1],'r.')
  #ADD  fig.canvas.draw()
  #ADD  plt.draw()
  #break

  #ADDax = fig.add_subplot(122)
  #ax.set_aspect('equal', 'box')
  #print np.arange(i)*dt
  #ADDplt.plot(np.arange(i)*dt,np.cumsum(absorbedvsT[0:i]))
  #ADDplt.xlim([0,dt*np.max(times)])
  #plt.ylim([0,200])
  #plt.xlabel("time") 
  #plt.ylabel("Cume sum") 
  #plt.gca().set_aspect('equal', 'box')
  fig.canvas.draw()
  plt.draw()

# Does nothing 
def noplotter(x,xN,i):
  1

# Can be overridden to add different BC, etc
class xN_conditions():
  def __init__(self):
    1

  # xN - updated x positions
  # taui - current time 
  # returns xN - updated positions
  def calc(self,x,xN,taui):
    #print "OLD"
    xp = xN  # inherited calc functions will replace this with periodic BC mapped
    dxp = xN -x # inherited calc functions will replace this with periodic BC mapped 
    return (xN,xp,dxp)

  # x - posn tau - timestep, D diff const 
  def calcforce(self,x,tau,D):
    return 0.

# dealing with obstables, reflections, etc 
class BDIterator():
 def __init__(self, 
  nParticles=2,
  T=1e4,
  D=0.1,
  dimensions=2,
  tau=1.,
  plot=False
  ):

  self.debug=0
  self.nParticles=nParticles 
  self.T=T
  self.D=D
  self.dimensions=dimensions
  self.tau=tau
  self.xN_conditions = xN_conditions()
  self.init_pos = np.zeros(self.dimensions)
  self.times = self.tau*(np.arange(self.T)+1)

  if(self.dimensions==1):
    drawplot = drawplot1d
  else:
    drawplot = drawplot2d

  if(plot):
    self.plotter=drawplot
  else:
    self.plotter = noplotter

  #T = 1000
 # run BD simulation 
 #@profile
 def run(self):
  tau = self.tau
  D = self.D   
  dimensions= self.dimensions
  times = self.times  
  #x = np.zeros([self.nParticles,self.dimensions])
  x = np.zeros([self.nParticles,self.dimensions])

  self.init_pos = np.array(self.init_pos).flatten()
  x[:] = self.init_pos #

  if (self.dimensions > 1 and np.shape(self.init_pos)[0]<2):
    print "OVeriding init_pos?"
    x[:,0] = self.init_pos #
    x[:,1] = 10 * np.random.rand( self.nParticles ) # need to make this general 
    #x[:] = [2,8]

  self.init_pos = x
  #print "CHEAT"; x[:]=0

  #x = np.repeat(self.nParticles,self.init_pos)         
  dxs = np.zeros([self.T,self.dimensions]) 
  dsqdn = []                    
  dsqdns=np.zeros(self.nParticles)


  #print "CHEAT"; x[:] = 0
 
  sqditau=0
  disp = self.T/100
  #n = np.zeros([self.nParticles,self.dimensions])
  #xN = np.zeros([self.nParticles,self.dimensions])
  #xp = np.zeros([self.nParticles,self.dimensions])
  #dxp= np.zeros([self.nParticles,self.dimensions])
  #gaussian(u0=np.random.rand(self.nParticles*self.dimensions),u1=np.random.rand(self.nParticles*self.dimensions)),

  import time
  prevTime = time.time()
  for i,taui in enumerate(times):
    ## Random number generator
    #n = np.array([np.random.randn(self.nParticles),np.random.randn(self.nParticles)]).T
    if(np.mod(i,disp)==0):
      print "Iter %d/%d (%d sec)"%(i,self.T, (time.time()-prevTime))
      prevTime = time.time()
      sp.cache.clear_cache()


    n = np.reshape(
      np.random.randn(self.nParticles*self.dimensions),
      [self.nParticles,self.dimensions])

    if(self.debug):
      np.random.seed(1)
      ns = np.random.randn(10000)
      #ns = np.random.randn(100)
      k = np.sqrt(2*D*tau)
      n = ns[i]
      #print n

    #print "CHEATING"; n[:,0]= 0; n[:,1]=1.
    

    ## VALIDATED FOR nParticles
    # use numbers from 1D simulation
    #n  = np.zeros([nParticles,2])
    #n[:,0] = nx[i]
    #n[:,1] = ny[i]

    ## propagate
    # random displacement 
    #dx= np.sqrt(D*dimensions*tau)*n
    k = np.sqrt(2*D*tau)
    #print "CHEATIN"; k=1
    dx= k*n
    #del n
    # force
    fdx = self.xN_conditions.calcforce(x,tau,D)
    #print " DFDFSSDF ", fdx
   # quit()

    xN = x  + fdx + dx
    #del fdx
    #xN = x + fdx

    ## Apply conditions
    #print "x",x
    #print "xN before", xN
    (xN,xp,dxp) = self.xN_conditions.calc(x,xN,taui) 
    #xN[:]=axN[:]; xp[:]=axp[:]; dxp[:] = adxp[:]; del axN; del axp; del adxp
    #print "dx",dx
    #print "xN after", xN

    ## store 
    # For more than one 
    #dsqdnsum += np.sum(np.multiply(dx,dx),axis=nParticles-1)
    #dsqditau = np.sum(np.multiply(dx,dx),axis=1)  # MSD for each particle 
    # we grab the displacements before periodic bcs etc are applied (called dxp) 
    dsqditau = np.sum(np.multiply(dxp,dxp),axis=1)  # MSD for each particle 
    # TODO 
    # Save dsqditau*dsqditau for comuting variance
    # compute MSF from XN-X0 with PBC, not dx 
    # To store ENTIRE trajectory dsqdn.append(dsqditau)
    # This stores the MSD per atom at each time step (so later we need to divide  by timesteps) 
    dsqdns += dsqditau 
    

    ## Plot
    self.plotter(x,xN,i)

    ## update 
    x = xN
 
    ###
    #print x 
    #print self.init_pos
    #print "tau", taui
    #print "xp",xp
    xpdiffsqd = np.sum((xp-self.init_pos)**2,axis=1)
    sqditau += xpdiffsqd/taui
    sqdf  = xpdiffsqd  
    #print "xpdiffsqd", xpdiffsqd 
    #print "sqditau", sqditau

  #print "DONE"
  self.dsqdns = dsqdns
  self.sqditau = sqditau
  self.sqdf = sqdf



  # deleting now, so that hopefully memory usage iproves  
  del self.xN_conditions

  #import pdb; pdb.set_trace()
  return(dsqdns)


  # calculate diffusion constant from stored means squared 'dx' displacement 
  # dx = x - x_prev
# Validated for periodic cases too
 def CalcDiffConst(self):
    # from ENTIRE trajectory Dest2d = np.mean( np.asarray(dsqdn).T,axis=1) / (2 * dimensions * tau )
    # divide by self.T since I take a running sum
    Dest2di = (self.dsqdns/np.float(self.T)) / (2 * self.dimensions * self.tau )
    #print "Dest2d each trj",Dest2di
    Dest= np.mean(Dest2di)

    return Dest

 # Calculate diff const from stored mean squared displacements 
 # xfinal - xoriginal
 def CalcDiffConstDisplacements(self):
    #print "SDFDS", self.sqditau
    tFinal = self.times[-1::]
    # divide by self.T since I take a running sum
    Dest2df = (self.sqdf) / (2 * self.dimensions * tFinal)     
    #print "from last iter",Dest2df
    Dest2di = (self.sqditau/np.float(self.T)) / (2 * self.dimensions)     
    #print "running avg",Dest2di
    Dest= np.mean(Dest2di)
    #print "Particle Dest", Dest
    return Dest
  
   
  
