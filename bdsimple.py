import numpy as np 

from bdIterator import *
from finiteLine import *
from finitePlane import *
from infinitePlane import *

# ERROR - 1/3-dim sims give diff const too slow/fast [see validate_BDsim(dimensions=1)]
# 4. return package with dp/fti at boundary  

# q0 is starting position, um 
# a is position of reflective buondary
# D is diff const [um^/ms]
# returns time, [ms]
def tauMFPT(q0,a=1, D=1):
  # from 
  # http://jcp.aip.org/resource/1/jcpsa6/v137/i16/p164108_s1?view=fulltext#citeref_c6
  # The first passage time from q0 to b is 
  # tau(q0) = 1/D [ int_q0^b dq [int_a^q dq' ] ]
  # where b is the absorbing boundary and a is the refelcting boundary 
  # assuming particle is absorbed at x=0
  # tau(q) = (1/D) (a*q0 - q0^2 /2)
  return (1/D) * (a*q0 - 0.5 * q0*q0)


# My test for 'non-iterated' BD simulaton of 1 particle 
def validationBD_single_DiffConstFromMSD():
  T=1e3
  D = 0.100 # [um^2/ms]       
  dimensions = 2;         # two dimensional BD1dSimulaton
  tau = 1;               # time interval in milliseconds
  p=np.zeros([T,dimensions])
  #k = np.sqrt(D*dimensions*tau)
  k = np.sqrt(2*D*tau)
  dx =k*np.random.randn(T) 
  dy =k*np.random.randn(T) 
  p[:,0] = np.cumsum( dx )
  p[:,1] = np.cumsum( dy )
  sqd = np.sum( np.multiply(p,p) , axis = 1) 
  #dsqd = np.multiply(dx,dx) + np.multiply(dy,dy)
  ds = p[1::,] - p[:-1:,]  # ds = p_1..N - p_0..(N-1)
  dsqd = np.sum(np.multiply(ds,ds),axis=1) 
  Dest = np.mean( dsqd ) / (2 * dimensions * tau )
  standardError = np.std( dsqd ) / ( 2 * dimensions * tau * np.sqrt(T) )
  actualError = D - Dest

  print "validation_DiffConstFromMSD: " 
  print "D: ", D
  print "Dest: ", Dest

  Dv = validation_DiffConstFromMSD(D=D,dimensions=dimensions,tau=tau,T=T)
  assert(np.abs(Dest-Dv) < D*0.05),"ERROR" 
  
def test(plot=False,z=1,name="traj.png", title="Attractive",box=True):
#def test(plot=False,z=1,name="traj.png", title="Attractive",box=True):
  T=1e2 # time steps 
  D = 0.1 # [um^2/ms]       
  dimensions = 1;         # two dimensional BD1dSimulaton
  tau = 0.1;             # time interval in milliseconds


  nParticles = 1 
  dimensions = 2 
  bdSim = BDIterator(nParticles,T=T,D=D,dimensions=dimensions,tau=tau,
    plot=plot)   
  boundaryReflect = 10 # x 
  q0 = 8
  #z = -10 # repulsive
  bdSim.init_pos = q0# x 
 # bdSim.xN_conditions = finitePlane(boundaryReflect,q0,nParticles,z=z,name=name,title=title,box=box)
  bdSim.xN_conditions = finitePlane(boundaryReflect,q0,nParticles,z=z,name=name,title=title,box=box)

  bdSim.run()
  Dest = bdSim.CalcDiffConst()
  print "Dest ", Dest


  return 
  # replace with new conditions

def testi(plot=False,z=1,name="traj.png", title="none",box=True,T=2.5e3,nParticles=100):
#  T=1e3 # time steps 
  D = 0.1 # [um^2/ms]       
  tau = 0.1;             # time interval in milliseconds
#  nParticles = 150
  dimensions = 2    

  #print "DEBUG"; T=1e1; nParticles=50 
  bdSim = BDIterator(nParticles,T=T,D=D,dimensions=dimensions,tau=tau,
    plot=plot)   
  boundaryReflect = 10 # x 

  q0 = 0
  #z = -10 # repulsive
  bdSim.init_pos = q0# x 
 # bdSim.xN_conditions = finitePlane(boundaryReflect,q0,nParticles,z=z,name=name,title=title,box=box)
  bdSim.xN_conditions = infinitePlane(boundaryReflect,q0,nParticles,z=z,name=name,title=title,box=box)

  bdSim.run()
  Dest = bdSim.CalcDiffConstDisplacements()
  #Dest = bdSim.CalcDiffConst()
  print title+"Dest ", Dest


  #del bdSim # I don't think this is necessary, but not sure where my mem issues are
  return Dest 
  # replace with new conditions

def excludedComparisons():
  # To generate figures
   
  N = 3
  titles=["8x8","6x6", "4x4", "2x2", "0x0"]                    
  c = len(titles)
  Ds=np.zeros([c,N],dtype=np.float)
  
  for i in np.arange(N):
    for j,t in enumerate(titles):
      name = "traji_"+titles[j]+".png"
      Ds[j,i]=testi(plot=False,z=0.0,name=name,title=titles[j], box=titles[j])

  print Ds
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  y = np.mean(Ds,axis=1)
  N = len(y)
  ind = range(N)
  err = np.std(Ds,axis=1)
  ax.bar(ind, y, facecolor='#777777', 
       align='center', yerr=err, ecolor='black')
  ax.set_ylabel('D [um2/s]')
  ax.set_title('Effective diffusion',fontstyle='italic')
  ax.set_xticks(ind)
  group_labels = titles
  ax.set_xticklabels(group_labels)
  fig.autofmt_xdate()

  plt.gcf().savefig("effdiff_size.png") 



def forceComparisons():
  # To generate figures
   
  c = 4
  N = 2
  Ds=np.zeros([c,N],dtype=np.float)
  titles=["Neutral", "Attractive","Repulsive","No box"] 
  
  for i in np.arange(N):
    Ds[1,i]=testi(plot=False,z= 0.2,name="traji_attractive.png",title=titles[1],box="3x3")            
    Ds[0,i]=testi(plot=False,z=0.0,name="traji_neutral.png",title=titles[0],box="3x3")       
    Ds[2,i]=testi(plot=False,z=-0.2,name="traji_repulsive.png",title=titles[2],box="3x3")             
    Ds[3,i]=testi(plot=False,z=0.0,name="traji_neutralmbox.png",title=titles[3],box=False)       

  print Ds
  fig = plt.figure()
  plt.title("Diff from MSD") 
  ax = fig.add_subplot(1,1,1)
  y = np.mean(Ds,axis=1)
  N = len(y)
  ind = range(N)
  err = np.std(Ds,axis=1)
  ax.bar(ind, y, facecolor='#777777', 
       align='center', yerr=err, ecolor='black')
  ax.set_ylabel('D [um2/s]')
  ax.set_title('Eff diff w obs (STD ERRRO WRONG)',fontstyle='italic')
  ax.set_xticks(ind)
  group_labels = titles
  ax.set_xticklabels(group_labels)
  fig.autofmt_xdate()

  plt.gcf().savefig("effdiff_p5.png") 


  
  #nParticles = 10
  #bdSim = BDIterator(nParticles,T=T,D=D,dimensions=dimensions,tau=tau,
  #  plot=plot)   
  #bdSim.init_pos = q0
  #bdSim.xN_conditions = finiteLine(boundaryReflect,nParticles)
#
#  bdSim.run()
#  Dest = bdSim.CalcDiffConst()
# Verified 
def validate_BDsim(dimensions=1):
  T=1e3
  D = 0.10 # [um^2/ms]       
  tau = 0.1;               # time interval in milliseconds
  nParticles = 500
  #D=1; tau=0.1; T=1e4; dimensions=3; nParticles=1000
  debug=0
  #D=.1; tau=0.1; T=1e3; dimensions=2; nParticles=3; debug=1

  bdSim = BDIterator(nParticles,T=T,D=D,dimensions=dimensions,tau=tau)
  bdSim.debug=debug
  bdSim.run()

  Dest = bdSim.CalcDiffConstDisplacements()
  assert(np.abs(Dest-D) < D*0.05),"ERROR D=%f" % Dest
  print Dest

  Destdx = bdSim.CalcDiffConst()
  print Destdx
  assert(np.abs(Destdx-D) < D*0.05),"ERROR D=%f" % Dest


# My test for iterated BD simulation of N particles 
def validationBD_DiffConstFromMSD():
  T=1e5
  D = 0.100 # [um^2/ms]       
  dimensions = 2;         # two dimensional BD1dSimulaton
  tau = 1;               # time interval in milliseconds

  # validated
  nParticles = 10
  bdSim = BDIterator(nParticles,T=T,D=D,dimensions=dimensions,tau=tau)
  bdSim.run()
  Dest = bdSim.CalcDiffConst()


  print "validation_DiffConstFromMSD: (multiple Trajs) " 
  print "D: ", D
  print "Dest: ", Dest

  Dv = validation_DiffConstFromMSD(D=D,dimensions=dimensions,tau=tau,T=T)
  assert(np.abs(Dest-Dv) < D*0.01),"ERROR" 

def validation1dMFPT():

  q0s = np.linspace(1.,9,17)
  D = 0.1
  D = 0.5
  boundaryReflect=10.

  
  nParticles=1e2; T=1e7; tau = 0.01  # 30s CPU 
  nParticles=1e4; T=1e6; tau = 0.1  # 30s CPU 
  nParticles=1e5; T=1e6; tau = 0.5  # production
  nParticles=1e2; T=1e4; tau = 0.01  # production
 

  tas = []
  ts = []

  q0 = 5.

  for q0 in q0s:
    print q0
    ta = tauMFPT(q0,boundaryReflect,D)
    print "Analy", ta                                
    tm,tstd,pts =BD1dSimulaton(D=D,boundaryReflect=boundaryReflect,q0=q0,nParticles=nParticles,T=T,tau=tau)
    completed = np.shape(pts)[0]
    print "Pred.(avg %f/std %f/completed %d)" % (tm,tstd,completed)     
    tas.append(ta)
    ts.append([tm,tstd,completed])

  tas = np.asarray(tas)
  ts = np.asarray(ts)
  print "VERIFY STD ERROR USED CORRECTLY"
  plt.figure()
  plt.plot(q0s,tas,"k-",label="Analytical")
  stdError = ts[:,1] / np.sqrt(ts[:,2])   # sqrt of num particles reaching abs. boundary  
  plt.errorbar(q0s,ts[:,0],yerr=1.96*stdError,fmt='b.',label="1D BD")         
  plt.ylabel("t [ms]") 
  plt.xlabel("q0 [um]") 
  title = "MFPT in 1D: D=%f [um**2/ms], absorbing at x=0" % (D)
  plt.title(title)
  plt.legend(loc=2)
  plt.gcf().savefig("MFPT1D.png") 


  



#def BD1dSimulaton(
#  D = 0.1, # [um^2/ms] 
#  boundaryReflect = 10, # [um]
#  tau=0.1,    # ms 
#  ion = False,# display 
#  T = 1000,   # time steps     
#  nParticles = 20, # particles 
#  q0 = 5.     # starting position 
#  ):
#
#  # Params 
#  boundaryAbsorb=0., # BY OUR SOLUTION OF THE MFPT ANALYTICAL FORMULATION
#  dimensions = 1 
#
#  # Inits 
#  x = np.repeat(q0,nParticles)
#  times = tau*(np.arange(T)+1)
#  absorbed = np.zeros(nParticles)
#  active   = np.ones(nParticles)
#
#  for i,taui in enumerate(times):                     
#    n = np.random.randn(nParticles)
#    n *= active   # will 'freeze' absorbed particles 
#    #dx = np.sqrt(dimensions*D*tau)*n
#    dx = np.sqrt(2*D*tau)*n
#
#    xN = x  + dx
#
#    # Check bounds 
#    #reflectedIdx = np.where(xN[:,0] < params.boxMin[0])
#    reflectedIdx = np.where(xN > boundaryReflect)       
#    absorbedIdx  = np.where(xN < boundaryAbsorb)        
#
#    # reflect
#    xR = xN[reflectedIdx] - boundaryReflect   
#    xN[reflectedIdx] = boundaryReflect  - xR
#
#    # abs
#    # we multiple by active, since if aprticle hasn't been absorbed,
#    # the tau will be stored  
#    absorbed[absorbedIdx] += active[absorbedIdx]* taui      
#    xN[absorbedIdx] = 0
#    active[absorbedIdx] = 0
#
#    #print "\ntime ", i
#    #print "Pos_o ", x 
#    #print "Active ", active 
#    #print "Absorbed ", absorbed 
#    #print "Pos_n ", xN
#
#
#    # plot 
#    if(ion):
#      drawplot(x,xN)   
#
#    # break if no more particles 
#    if (np.sum(active) < 0.01 * nParticles):
#      print "Leaving with %d/%d particles " % (np.sum(active) ,nParticles)
#      break 
#
#    # update 
#    x = xN
#
#  ## TAKE EM
#  absorbedTimes = absorbed[ np.where(absorbed > 0 ) ]
#  if(np.shape(absorbedTimes)[0] ==0):
#    return (-1,-1,0)
#  mfpt = np.mean(  absorbedTimes )
#  mfpt_std = np.std( absorbedTimes ) 
#
#  return (mfpt, mfpt_std,absorbedTimes)

# http://labs.physics.berkeley.edu/mediawiki/index.php/Simulating_Brownian_Motion
def validation_DiffConstFromMSD(
  D = 4.2902e-013,  
  dimensions = 2,         # two dimensional BD1dSimulaton
  tau = .1,               # time interval in seconds
  T = 1000
  ):

  # single particle 
  #T=100 
  p=np.zeros([T,dimensions])
  nx =np.random.randn(T) 
  ny =np.random.randn(T) 
  #k = np.sqrt(D*dimensions*tau)
  k = np.sqrt(2*D*tau)
  dx =k*nx
  dy =k*ny
  p[:,0] = np.cumsum( dx )
  p[:,1] = np.cumsum( dy )
  sqd = np.sum( np.multiply(p,p) , axis = 1) 
  dsqd = np.multiply(dx,dx) + np.multiply(dy,dy)
  Dest = np.mean( dsqd ) / (2 * dimensions * tau )
  standardError = np.std( dsqd ) / ( 2 * dimensions * tau * np.sqrt(T) )
  actualError = D - Dest

  print "D: ", D
  print "Single Particle Dest: ", Dest

  return Dest

# Validation using BD simulator w periodic BC
#@profile
def validate_MSD_Periodic():

  D = 0.1
  z = 0 
  T = 1e4
  nParticles = 150

  #D=.1; tau=0.1; T=1e3; dimensions=2; nParticles=2; debug=1
  #print "DEBUG" D=.1; tau=0.1; T=1e2; dimensions=2; z=1; nParticles=20; debug=0
  Dest = testi(plot=False,z=z,name="traj.png",box=False, T=T, nParticles=nParticles)


#
#  bdSim = BDIterator(nParticles,T=T,D=D,dimensions=dimensions,tau=tau,
#    plot=False) 
#  q0 = 0
  #z = -10 # repulsive
#  bdSim.init_pos = q0# x 
#  boundaryReflect=10
#  #bdSim.xN_conditions = finitePlane(boundaryReflect,q0,nParticles,z=0,name="name",title="title",box=False)
#
#  print "DEBUG"; T=1e1; nParticles=20; z=1
#  bdSim.xN_conditions = infinitePlane(boundaryReflect,q0,nParticles,z=z,box=False)
#  bdSim.debug = debug
#
#  bdSim.run()

#  Dest = bdSim.CalcDiffConstDisplacements()
  print "Validation Dest ", Dest
  assert(np.abs(Dest-D) < D*0.10),"ERROR" 
#  Destdx = bdSim.CalcDiffConst()
#  print "Validation Destdx ", Destdx
#  assert(np.abs(Destdx-D) < D*0.05),"ERROR" 



import sys
#
# Revisions
#       10.08.10 inception
#


if __name__ == "__main__":
  msg="""
Purpose: 

Usage:
  script.py -mfpt/-msd

Notes:
"""
  remap = "none"


  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"

  for i,arg in enumerate(sys.argv):
    if(arg=="-validate"):   
      # VALIDATED
      #validate_MSD_Periodic()
      validation_DiffConstFromMSD()        
      validationBD_single_DiffConstFromMSD()        
      validationBD_DiffConstFromMSD()        
      validate_BDsim(dimensions=1)
      print "CODE IS VALIDATED - commit me!"
    if(arg=="-test"):
      test()
    if(arg=="-forceComparisons"):
      forceComparisons()
    if(arg=="-excludedComparisons"):
      excludedComparisons()
    if(arg=="-mfpt"):  
      # UNVALIDATED
      validation1dMFPT()

    if(arg=="-misc"):
      validate_MSD_Periodic()
      #validate_BDsim(dimensions=1)()
    if(arg=="-misc2"):
      validate_BDsim(dimensions=2)
      #validate_BDsim(dimensions=1)()


  #raise RuntimeError(msg) 






