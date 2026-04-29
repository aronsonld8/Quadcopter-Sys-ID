from integrator import Integrator
from propeller import Propeller

class Quadcopter:
  def __init__(self, int_instance: Integrator):
    self.mot = []
    self.Vel_CM   = np.empty(3)
    self.att_body = np.empty(3)
    
    self.mot.append(mot_params['r1'], mot_params['d1'], int_instance))
    self.mot.append(mot_params['r2'], mot_params['d2'], int_instance))
    self.mot.append(mot_params['r3'], mot_params['d3'], int_instance))
    self.mot.append(mot_params['r4'], mot_params['d4'], int_instance))

    self.Body_Area = [body_params['Ax'], body_params['Ay'], body_params['Az']]
    self.Cd        = [body_params['Cx'], body_params['Cy'], body_params['Cz']]
  
  def SetVelocity(self, motOmega, Vel_CM, att_body):
    self.Vel_CM = Vel_CM
    self.att_body = att_body
    for i in range(4):
      self.mot[i].SetVelocity(motOmega[i], Vel_CM, att_body)

  def GetForce(self):
    Fp = np.zeros(3)
    Fd = 0.5*rho_air*self.Vel_CM*np.abs(self.Vel_CM)
    #prop
    for i in range(4):
      Fp += self.mot[i].GetForce()
    #drag
    for i in range(3):
      Fd[i] *= self.Body_Area[i]*self.Cd[i]
    return Fp - Fd

  def GetTorque(self):
    Tau = np.zeros(3)
    #prop
    for i in range(4):
      Tau += self.mot[i].GetTorque()
      Tau += np.cross(self.mot[i].r,self.mot[i].GetForce())
    return Tau
    
