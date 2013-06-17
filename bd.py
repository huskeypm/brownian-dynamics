import numpy as np
import sympy as sp
import matplotlib.pylab as plt

#
#  Todo: 
#    add in correct random displacements
#    use correct BD equation
#    compute mean squared displacement and relate to t
#    consider completely periodic system 
#
 
#x0 = x.copy()
#
#MSD = (xN - x0)**2
#tau = MSD/(2*D*t)

ion=0
#params.nParticles=4
#tF = 500
tF = 100
dt = 0.1

gradx=1
grady=1

class params:
  fig=0
  nParticles=10
  kappa = 1. # like Debye length
  D = 1.
  kT = 1.
  # back boundary is at x=-1
  # abs boundary is at x=1
  # y is periodic 
  #boxMin = np.array([-1,-1])
  #boxMax = np.array([ 1, 1])
  boxMin = np.array([-1e-1,-1])
  boxMax = np.array([ 1e-1, 1])

  # params.obstacle 
  obsMin = np.array([  0.25,-0.5])
  obsMax = np.array([  0.75, 0.5])


def calcobs():
    params.obs = np.array([ [params.obsMin[0],params.obsMin[1]],[params.obsMax[0],params.obsMin[1]],[params.obsMax[0],params.obsMax[1]],[params.obsMin[0],params.obsMax[1]],
                  [params.obsMin[0],params.obsMin[1]]])

#  // returns a gaussian with 0 mean and unit variance
def gaussian(u0=np.random.rand(1),u1=np.random.rand(1)):
    #    // Box-Muller transform
    #    // uniform() returns a uniformly distributed number strictly between 0 and 1
    pi2 = 2*np.pi 
    #     double u0 = uniform();
    #    double u1 = uniform();
    lterm = np.sqrt(-2*np.log( u0));
    return lterm*np.cos( pi2*u1);

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






namespace = sp.__dict__.copy()
def force():            
  #vxc=np.array([0,0])
  vxc=(params.obsMin+params.obsMax)/2. 
  #print vxc


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
  expr = sp.Eq(sp.exp(-kappa*r)/r,V)
  exprs = expr.subs({r:(xp-xr)**2+(yp-yr)**2})
  exprs = exprs.subs({xr:vxc[0],yr:vxc[1]})    
  
  sol = sp.solve(exprs,V)[0]
  gradx = sol.diff(xp)
  grady = sol.diff(yp)
  namespace['xp']=xp
  namespace['yp']=yp
  namespace['kappa']=kappa
  return(gradx,grady,sol)

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


def forceeq(vx=1,vy=1,vkappa=0.):
  vgradx=params.gradx.subs({namespace['xp']:vx,namespace['yp']:vy,namespace['kappa']:vkappa}) 
  vgrady=params.grady.subs({namespace['xp']:vx,namespace['yp']:vy,namespace['kappa']:vkappa}) 

  return vgradx,vgrady

def init(num):
  randPos = (np.random.rand(num)-0.5) * (params.boxMax[0]-params.boxMin[0])
  x = np.zeros([num,2])
  x[:,0] = 0.
  x[:,1] = randPos


  calcobs()

  params.gradx,params.grady,eqn = force()
  # could use this namespace=np.__dict__.copy()
  params.vforce = np.vectorize(forceeq,otypes=[np.float,np.float])
  #params.vforce([1,1],[2,2])

  return x 

def simulate(x,s=1,f=0,forceScale=1e-2,Dt=100,dt=0.1):
  absorbed = np.zeros(params.nParticles)
  absorbedvsT=np.zeros(int(Dt/dt))
  absorbedIdx = []
  times =  np.arange(Dt/dt)
  
  
  #x[1]=0.5      
  absCtr=0
  for i in times:
     # replace absorbed and put at 'starting' point
     xI = init( np.sum(absorbed) )
     x[absorbedIdx] = xI
     absorbed[absorbedIdx] = 0
  
     ## force calc (vforce is a vectorized smpy function for the gradient of the potential
     params.kappa=1.
     F = np.transpose(params.vforce(x[:,0],x[:,1],params.kappa))
     F *=forceScale
     #print "force",forceTerm
     
  
     ## step
     # uniform n = np.random.rand(params.nParticles*2)-0.5
     #n = np.reshape(n,[params.nParticles,2])
     nx = gaussian(u0=np.random.rand(params.nParticles),u1=np.random.rand(params.nParticles))
     ny = gaussian(u0=np.random.rand(params.nParticles),u1=np.random.rand(params.nParticles))
     n = np.reshape([nx,ny],[params.nParticles,2])

     xN = x + (params.D/params.kT)*F*dt + np.sqrt( 2*params.D*dt)*n                   
     #print "x",x 
     #print "xN",xN
  
     # DBG
     #x[0,1]=0.4
     #xN = x + [ 0.10,0]         
     #xN = x + [0,0.25]         
     #xN = x + [0.1,0]         
     #print "Iter %d" % i
     #print "Xn", xN
     #x=xN
  
     ## check bounds 
     reflectedIdx = np.where(xN[:,0] < params.boxMin[0]) # xN[0,:] = [-1.25,0] 
     absorbedIdx  = np.where(xN[:,0] > params.boxMax[0]) # xN[1,:] = [1.25,0] 
     btmPeriodIdx = np.where(xN[:,1] < params.boxMin[1]) # xN[2,:] = [0,1.25]
     topPeriodIdx = np.where(xN[:,1] > params.boxMax[1]) # xN[3,:] = [0,-1.25]
  
     # reflect (validated)
     xR = xN[reflectedIdx,0] - params.boxMin[0]
     xN[reflectedIdx,0] = params.boxMin[0] - xR
     #print "Ref", xN
  
     # periodic (validates)
     xB = xN[btmPeriodIdx,1] - params.boxMin[1]
     xN[btmPeriodIdx,1] = xB + params.boxMax[1]
     xT = xN[topPeriodIdx,1] - params.boxMax[1]
     xN[topPeriodIdx,1] = xT + params.boxMin[1]
     #print "B/T Period", xN
  
     # absorbed (validated) 
     absorbed[absorbedIdx] = 1
     nabsorbed = np.sum(absorbed)
     absorbedvsT[i] = nabsorbed            
     absCtr+=nabsorbed     
     if(np.mod(i*dt,1)==0):
       print "abs: %f %d/%d" % (i*dt,nabsorbed,  absCtr)
     #print "absorbed", xN[absorbedIdx,:]
     # do statisticss, but only with non-absorbed
  
     # params.obstacle
     # either reflect or redraw. Easiest to reflect
     # top: (top-x)*(x-bottom) > 0 if within bounds 
     # left: (right-x)*(x-left) > 0 if within bounds 
     #xN[0,:] = [params.obsMax[0]-0.001,0]
     #xN[1,:] = [0, params.obsMax[1]-0.001]
     #xN[2,:] = [params.obsMax[1]-0.001, params.obsMax[1]-0.001]
     inX = (params.obsMax[0]-xN[:,0])*(xN[:,0]-params.obsMin[0]) >= 0
     inY = (params.obsMax[1]-xN[:,1])*(xN[:,1]-params.obsMin[1]) >= 0
     inBox = np.where(inY*inX)
     #print "Collision", xN[inBox,]  
     # replace 'collision with original
     xO = xN.copy()
     xN[inBox] = x[inBox] 
     
     # update
     x = xN
  
     # plot 
     if(ion and np.mod(i,10)==0):
       fig = params.fig
       fig.clear()
       ax = fig.add_subplot(121)
       ax.set_aspect('equal', 'box')
       plt.title(params.name)
       plt.xlim([-1,1]) #[params.boxMin[0],params.boxMax[0]])
       plt.ylim([-1,1]) #[params.boxMin[1],params.boxMax[1]])
       plt.plot(params.obs[:,0],params.obs[:,1],'k-')
       #plt.plot(xN[outBox,0],xN[outBox,1],'b.')
       plt.plot(xN[:,0],xN[:,1],'b.')
       showCollisions=0
       if(showCollisions==1 and np.shape(inBox)[0] > 0):
         plt.plot(xO[inBox,0],xO[inBox,1],'r*')
         plt.plot(xN[inBox,0],xN[inBox,1],'r.')
         fig.canvas.draw()
         plt.draw()
         #break

       ax = fig.add_subplot(122)
       #ax.set_aspect('equal', 'box')
       #print np.arange(i)*dt
       plt.plot(np.arange(i)*dt,np.cumsum(absorbedvsT[0:i]))
       plt.xlim([0,dt*np.max(times)])
       #plt.ylim([0,200])
       #plt.xlabel("time") 
       #plt.ylabel("Cume sum") 

       #plt.gca().set_aspect('equal', 'box')
       fig.canvas.draw()
       plt.draw()
       plt.gcf().savefig(params.name+"traj.png")
     #plt.pause(0.5)
  return times,absorbedvsT
  
    


  

def doit(mode="neutral",boxPosition="-centered"):
   
  # plot
  if(ion):
    plt.ion()
    params.fig=plt.figure()
  
  
  # flags
  s = 10  # should use fluct/diss
  f = 0

  # obstacle 
  if(boxPosition=="-centered"):
    params.obsMin = np.array([  0.25, -0.25])
    params.obsMax = np.array([  0.75,  0.25])
  elif(boxPosition=="-offset"):
    params.obsMin = np.array([  0.25,  0.25])
    params.obsMax = np.array([  0.75,  0.75])
  else:
    raise RuntimeError("asdf")

    

  
  # step
  t=0
  Dt = tF-t

  # run 
  s=1
  f=1
  x = init(params.nParticles)
  if(mode=="attractive"):
    forceScale=5e-1     # attractive 
  elif(mode=="repulsive"):
    forceScale=-5e-1     # repulsive   
  elif(mode=="neutral"):
    forceScale=0.        # neutral
  else:
    raise RuntimeError("mode not known") 

  times,abs2 = simulate(x,s,f,forceScale=forceScale,Dt=Dt)
  
  
  

def validate():
  # step
  t=0
  Dt = tF-t
  
  s=1
  f=0
  x = init(params.nParticles)
  times,abs1 = simulate(x,s,f,Dt=Dt)
  
  s=1
  f=1
  x = init(params.nParticles)
  forceScale=5e-2
  times,abs2 = simulate(x,s,f,forceScale=forceScale,Dt=Dt)
  
  s=1
  f=1
  x = init(params.nParticles)
  forceScale=-5e-2
  times,abs2b = simulate(x,s,f,forceScale=forceScale,Dt=Dt)
  
  s=1
  f=0
  x = init(params.nParticles)
  forceScale=1e-3
  params.obsMin = np.array([  0.45,-0.1])
  params.obsMax = np.array([  0.55, 0.1])
  times,abs3 = simulate(x,s,f,forceScale=forceScale,Dt=Dt)
    
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


