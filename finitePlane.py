#
# Special class for 2D finite plane
# absorbing/reflect BC problem
#
from bdIterator import *
import numpy as np
import sympy as sp

namespace = sp.__dict__.copy()

class params:
  name ="validate"
  #boxMin = np.array([-1,-1])
  #boxMax = np.array([ 1, 1])
  print "WRN remove"
  boxMin = np.array([-1e-1,-1])
  boxMax = np.array([ 1e-1, 1])

  # params.obstacle 
  obsMin = np.array([2,2])
  obsMax = np.array([8,8])

  fig=plt.figure()

  1



class finitePlane(xN_conditions):
  def __init__(self,
    boundaryReflect,
    q0,
    nParticles,
    z=0,
    name="traj.png",
    title="Traj",
    box=True
    ):

    
    self.boundaryAbsorb = 0.
    self.q0 = q0 # origin (x) for particles 
    self.boundaryReflect = boundaryReflect # for the x-coordinate only. We're assuming particle is absorbed/reflected on edges
    self.nParticles=nParticles
    self.absorbed = np.zeros(nParticles)
    self.active   = np.ones(nParticles)
    self.inactive   = np.zeros(nParticles)#,1))
    #absorbedvsT=np.zeros(int(Dt/dt))
    self.absorbedIdx = []
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
        self.obsMax = np.array([5.1,5.1])



    # params 
    self.kT = 0.59     # energy [kcal/mol]
    self.z = z         # unit charge 


    self.times =[]
    self.absorbedvsT=[]

    ## Initialize everything
    self.init(nParticles)

    self.rgtPeriod = np.zeros(nParticles) # always zero
    self.lftPeriod = np.zeros(nParticles) # always zero
    self.btmPeriod = np.zeros(nParticles)
    self.topPeriod = np.zeros(nParticles)

  def plottinghack(self,xO,xN,inBox):
    # can quess time step by length of arrays (stupid trick)(
    i = np.shape(np.asarray(self.absorbedvsT))[0]


    mod = 100
    if(np.mod(i,mod)!=0):
      return

    fig = params.fig
    fig.clear()
    ax = fig.add_subplot(121)
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

    # show flux 
    ax = fig.add_subplot(122)
    absorbvsTime= np.cumsum(np.asarray(self.absorbedvsT))
    # not quite correct, but taking time averaged flux 
    smoothed = np.convolve(absorbvsTime, np.ones(mod)/float(mod),mode="valid")
    #print absorbvsTime
    #print smoothed          
    smoothed = smoothed[::mod]
    #print smoothed          
    #dAbsdT = absorbvsTime[1::] - absorbvsTime[:-1:]
    dAbsdT = smoothed[1::] - smoothed[:-1:]
    #ax.set_aspect('equal', 'box')
    #plt.plot(np.arange(i)*dt,np.cumsum(self.absorbedvsT[0:i]))
    plt.title("Absorption rate vs timeframe") 
    plt.plot(np.arange(np.shape(dAbsdT)[0])*mod,dAbsdT)
    #plt.xlim([0,dt*np.max(times)])
    #plt.ylim([0,200])
    #plt.xlabel("time") 
    #plt.ylabel("Cume sum") 

    #plt.gca().set_aspect('equal', 'box')
    fig.canvas.draw()
    plt.draw()
    plt.gcf().savefig(self.name)

    #raw_input("Enter key") 


  def calcobs(self):
    self.obs = np.array([ [self.obsMin[0],self.obsMin[1]],
                          [self.obsMax[0],self.obsMin[1]],
                          [self.obsMax[0],self.obsMax[1]],
                          [self.obsMin[0],self.obsMax[1]],
                          [self.obsMin[0],self.obsMin[1]]])

  def init(self,num):
    x = np.zeros([num,2])
    x[:,0] = self.q0
    randPos = (np.random.rand(num)-0.5) * (self.boxMax[1]-self.boxMin[1])
    x[:,1] = randPos

    #print "CHEAT: override boudnaryabsorb"; self.boundaryAbsorb = -1.; x[:]=0

    self.calcobs()
  
    self.gradx,self.grady,eqn = self.force()
    # could use this namespace=np.__dict__.copy()
    self.vforce = np.vectorize(self.forceeq,otypes=[np.float,np.float])
    #params.vforce([1,1],[2,2])
  
    return x 

  #@profile
  def calcforce(self,x,dt,D):
     #return np.array([1,2])
     
     ## add in force
     F = np.transpose(self.vforce(x[:,0],x[:,1]))
     #sp.cache.clear_cache() # help w memory leak due to sympy cache
     #F *= self.forceScale
     #print "x",x 
     #print "F",F
     #print "dt",dt
     return (D/self.kT)*F*dt


  def calc(self,x,xN,taui):
     ## Apply obstacle    
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

     ## apply PBC
     xp = xN


     ## check bounds 
     reflectedIdx = np.where(xN[:,0] > self.boxMax[0]) # xN[0,:] = [-1.25,0] 
     absorbedIdx  = np.where(xN[:,0] < self.boxMin[0]) # xN[1,:] = [1.25,0] 
     btmPeriodIdx = np.where(xN[:,1] < self.boxMin[1]) # xN[2,:] = [0,1.25]
     topPeriodIdx = np.where(xN[:,1] > self.boxMax[1]) # xN[3,:] = [0,-1.25]
     self.btmPeriod[btmPeriodIdx]+=1
     self.topPeriod[topPeriodIdx]+=1

     # reflect (validated)
     #print "orig", xN
     xR = xN[reflectedIdx,0] - self.boxMax[0]
     xN[reflectedIdx,0] = self.boxMax[0] - xR
     #print "Ref", xN

     # periodic (validates)
     xB = xN[btmPeriodIdx,1] - self.boxMin[1]
     xN[btmPeriodIdx,1] = xB + self.boxMax[1]
     xT = xN[topPeriodIdx,1] - self.boxMax[1]
     xN[topPeriodIdx,1] = xT + self.boxMin[1]
     #print "B/T Period", xN

     # absorbed (validated) 
     self.absorbed[absorbedIdx] = 1
     nabsorbed = np.sum(self.absorbed)
     self.absorbedvsT.append(nabsorbed)
     self.times.append(taui)
     self.absCtr+=nabsorbed
     #o if(np.mod(i*dt,1)==0):
     #if 1:
     #  print "abs: %f %d/%d" % (taui,nabsorbed,  self.absCtr)

     # replace absorbed and put at 'starting' point
     #print "in",xN 
     xI = self.init( np.sum(self.absorbed) )
     xN[absorbedIdx] = xI
     self.absorbed[absorbedIdx] = 0
     #print "out",xN 

     ## apply PBC
     unwrap =np.array([self.rgtPeriod*self.boxD[0] - self.lftPeriod*self.boxD[0],self.topPeriod*self.boxD[1] - self.btmPeriod*self.boxD[1]])
     #print "unwra", unwrap
     xp = xN.copy(); xp[:,0] += unwrap[0]; xp[:,1] += unwrap[1] # + unwrap
     #print "x",x
     #print "xN",xN
     #print "xp",xp
     #print "SDS"

     self.plottinghack(xO,xN,inBox)
     return (xN,xp,dxp)
  


  def forceeq(self,vx=1,vy=1):
    ##return 1,2               
    vgradx=self.gradx.subs({namespace['xp']:vx,namespace['yp']:vy}) 
    vgrady=self.grady.subs({namespace['xp']:vx,namespace['yp']:vy}) 
    return vgradx,vgrady

  # NOTE: this isn't quite correct, since we are 
  # dealing with a square obstacle, not a sphereical protein 
  # WARNING: Current [um] is the standard unit, but for the ESP
  #          and the homogenization sims, I need to use [A]. 
  #          Here I will just assume [A]
  def force(self):            
    #vxc=np.array([0,0])
    vxc=(self.obsMin+self.obsMax)/2. 
    molRad = np.mean((self.obsMax-self.obsMin)/2.) # see WARNING above [assumed A]
   

    # ec / 4 pi eps = 14.3840 [V Angstoms] 
    # --> eco4pieps = 14.384e3 [mV Angstroms]
    eco4pieps = 14.384e3 # conversion factor [mV Angstroms]
    #eco4pieps = 14.384 # [V Ang]
    ec = 8.854187817e-12 # electric constant [C^2/(Jm)]
    M_TO_ANG = 1e-10
    J_TO_KCAL = 0.000239005736
    ec = ec / (M_TO_ANG * J_TO_KCAL) # ec [C^2/kcal A]
    epsilonExterior = 80. # dielectric constant in exterior []
    ## System-specific parameters

    # ion
    ionC = .150 # ion conc [M]
    vkappa = 0.328 * np.sqrt(ionC) # inverse debye length for monovalent specieds [1/A]
    ikappa  = 1/vkappa # Debye length [A], "Intermoplecular and surface forces, Israelachvili" 

    prefac=self.z * np.exp(vkappa*molRad) * eco4pieps
    prefac/=epsilonExterior*(1+vkappa*molRad)



  
  
    xp=sp.Symbol('xp')
    xr=sp.Symbol('xr')
    yp=sp.Symbol('yp')
    yr=sp.Symbol('yr')
    r=sp.Symbol('r')
    kappa=sp.Symbol('kappa')
    V=sp.Symbol('V')
  
    # coulomb expr = sp.Eq(1/r,V)
    # Debye 
    #print "Debye eqution is wrong!!"
    expr = sp.Eq(prefac*sp.exp(-kappa*r)/r,V)
    exprs = expr.subs({r:(xp-xr)**2+(yp-yr)**2})
    exprs = exprs.subs({xr:vxc[0],yr:vxc[1],kappa:vkappa})    
    
    sol = sp.solve(exprs,V)[0]
    gradx = sol.diff(xp)
    grady = sol.diff(yp)
    namespace['xp']=xp
    namespace['yp']=yp
#    namespace['kappa']=kappa
    return(gradx,grady,sol)
  
  
  
  

#
#  Todo: 
#    add in correct random displacements
#    use correct BD equation
#    compute mean squared displacement and relate to t
#    consider completely periodic system 
#
 
#  
#  
#  for (int k = 0; k < 3; k++)
#    x[k] += (D/kT)*F[k]*dt + sqrt( 2*D*dt)*gaussian();
#  
#  // x - position
#  // D - diffusivity
#  // kT - Boltzmann times absolute temperature
#  // F - force on particle, compute every timestep
#  // dt - time-step size
#  
#  // Make sure that the force does not change more than a few percent between
#  // time steps. If it does, use a smaller timestep.  Also, you do want to 
#  // take small steps near the absorbing boundary.  I usually take it to
#  // be dt < L*L/(18*D), where L is the distance to the boundary, but still
#  // put limits on how small it gets, by setting a smallest possible stepsize
#  // to be dt > d*d/(18*D), where d is a very small length scale of the 
#  // problem.  So, if the length of the box is 1, d might be something like
#  // 10e-6.
#  
#  // I do have some more sophisticated code that can adjust the time step
#  // on the fly; if you need that, let me know and I can help you integrate
#  // it into your software.
#  
#  // If you are computing electrostatic forces from a grid, let me know
#  // and I can send you some code for that as well.






def test():
  dx,dy,eq = force()                 
  def feq(r,vkappa=0.):
    return eq.subs({sp.Symbol('xp'):r,sp.Symbol('yp'):0,sp.Symbol('kappa'):vkappa})
  
  vfeq = np.vectorize(feq,otypes=[np.float])
  p  = np.linspace(0,1,100)
 
  plt.figure()
  plt.plot(p,vfeq(p,0),"k-",label="dilute")
  plt.plot(p,vfeq(p,1),"k.-",label="semi-dilute")
  plt.plot(p,vfeq(p,10),"k.",label="conc")
  plt.ylim([0,10])
  plt.legend()
  plt.show() 






def simulate(s,z=0):

  #x = init(params.nParticles)
  T = 1e2  
  D = s      
  tau = 0.1
  nParticles = 100
  dimensions = 2
  bdSim = BDIterator(nParticles,T=T,D=D,dimensions=dimensions,tau=tau,
    plot=False)
  boundaryReflect = 10 # x 
  q0= 5 # x 
  fp = finitePlane(boundaryReflect,q0,nParticles,z=z)
  fp.obsMin = params.obsMin
  fp.obsMax = params.obsMax
  bdSim.xN_conditions = fp
 
  
  bdSim.run()
 
 
  times = np.asarray(fp.times)
  absorbvsTime= np.asarray(fp.absorbedvsT)   
  return times,absorbvsTime


  

def doit(mode="neutral",boxPosition="-centered"):
   
  # plot
  #if(ion):
  #  plt.ion()
  #  params.fig=plt.figure()
  
  
  # flags
  s = 10  # should use fluct/diss
  f = 0

  # obstacle 
  if(boxPosition=="-centered"):
    params.obsMin = np.array([  4, 6])
    params.obsMax = np.array([  4, 6])
  elif(boxPosition=="-offset"):
    params.obsMin = np.array([  6, 6])
    params.obsMax = np.array([  8, 8])
  else:
    raise RuntimeError("asdf")
  
  # step
  t=0
  Dt = tF-t

  # run 
  s=1
  if(mode=="attractive"):
    z=5e-1     # attractive 
  elif(mode=="repulsive"):
    z=-5e-1     # repulsive   
  elif(mode=="neutral"):
    z=0.        # neutral
  else:
    raise RuntimeError("mode not known") 

  times,abs2 = simulate(s,f,z=z)
  
  
  

def validate():
  
  s=1
  z=0
  times,abs1 = simulate(s,z=z)
  
  s=1
  z=0.1                 
  times,abs2 = simulate(s,z=z)
  
  s=1
  z=-0.1                 
  times,abs2b = simulate(s,z=z)
  
  s=1
  z=0             
  params.obsMin = np.array([  4,4])
  params.obsMax = np.array([  4.1, 4.1])
  times,abs3 = simulate(s,z=z)
    
  plt.figure()
  #plt.plot(times,absorbedvsT/np.float(params.nParticles),'k-')
  cm1 = np.cumsum(abs1)
  cm2 = np.cumsum(abs2)
  cm2b = np.cumsum(abs2b)
  cm3 = np.cumsum(abs3)
  plt.plot(times,cm1/np.max(cm3),'k-',label="-F")
  plt.plot(times,cm2/np.max(cm3),'b--',label="+F(<0)")
  plt.plot(times,cm2b/np.max(cm3),'r--',label="+F(>0)")
  plt.plot(times,cm3/np.max(cm3),'k.-',label="-box") 
  plt.title(params.name+": Absorbed number vs time") 
  plt.legend()
  plt.gcf().savefig(params.name+"absvstime.png")
  
  


   
   





   
   



import sys
#
# Revisions
#       10.08.10 inception
#


if __name__ == "__main__":
  msg="""
Purpose: 

Usage:
  script.py <-interactive> <-run/-validate> <-force [neutral/attractive/repulsive]> 

Notes:
"""
  remap = "none"

  runMode="normal"
  boxPosition="-centered"
  mode = "neutral"

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"

  for i,arg in enumerate(sys.argv):
    if(arg=="-interactive"):
      ion = True
    if(arg=="-validate"):      
      #ion = True
      runMode="validate"
    if(arg=="-offset"):
      boxPosition="-offset"

    if(arg=="-force"):   
      mode = sys.argv[i+1]  

    if(arg=="-name"):   
      params.name = sys.argv[i+1]  





  if(runMode=="validate"):
    validate()
  else:
    doit(mode=mode,boxPosition=boxPosition)


