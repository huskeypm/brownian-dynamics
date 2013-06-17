#
# Special class for 1D finite line
# absorbing/reflect BC problem
#

from bdIterator import *

class finiteLine(xN_conditions):
    def __init__(self,boundaryReflect,nParticles):
      self.boundaryAbsorb = 0.
      self.boundaryReflect = boundaryReflect
      self.nParticles=nParticles
      self.absorbed = np.zeros(nParticles)
      self.active   = np.ones(nParticles)
      self.inactive   = np.zeros(nParticles)#,1))

    def calc(self,x,xN,taui):
      ## store dx before periodic BC are applied
      dxp = xN - x

      ## apply PBC
      xp = xN

      # check bounds (use [0] I think because how xN is shaped)  
      reflectedIdx = np.where(xN > self.boundaryReflect)[0]
      absorbedIdx  = np.where(xN < self.boundaryAbsorb)[0]

      #print self.boundaryAbsorb  
      #print xN
      #print reflectedIdx
      #print absorbedIdx
      #print self.absorbed
      #print self.active       
      #print taui

      # reflect 
      xR = xN[reflectedIdx] - self.boundaryReflect
      xN[reflectedIdx] = self.boundaryReflect  - xR 

      # absorb
      # We multiply by 'active', which is unity only
      # if particle hasn't been absorbed yet
      # otherwise we ignore tau (e.g. ==0)
      self.absorbed[absorbedIdx] += self.active[absorbedIdx]* taui 

      self.active[absorbedIdx] = 0  
      #self.inactive[absorbedIdx,0] = 1  
      self.inactive[absorbedIdx] = 1  
      inactIdx =np.where(self.inactive==1)[0]
      #xN[(np.where(self.inactive)==1)[0]] = self.boundaryAbsorb # resets ALL previous absorbed to boundary 
      xN[inactIdx] = self.boundaryAbsorb # resets ALL previous absorbed to boundary 


      #if (np.sum(self.active) < 0.01 * self.nParticles):
      # print "Leaving with %d/%d particles " % (np.sum(active) ,nParticles)
      # break  
      
      
      
      return(xN,xp,dxp)

  
