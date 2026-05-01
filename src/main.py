import parameters
from integrator import Integrator
from propeller import Propeller
from quadcopter import Quadcopter


def set_params():
  #for blade
  blade_params['Area']= np.pi*blade_params['radius']**2
  blade_params['e'] = blade_params['ef']*blade_params['radius']
  blade_params['Ib'] = blade_params['mb'] * blade_params['cg']**2   # blade moment of inertia around hinge (estimate)
  blade_params['Ib'] = blade_params['Ib'] + 1/12*blade_params['mb']*((blade_params['radius']/2)**2) # add moment of lenghty rod
  blade_params['Mb'] = blade_params['mb'] * acc_g * blade_params['cg']            # moment of blade around hinge (estimate)
  blade_params['theta_0'] = blade_params['pitch']*np.pi/180 # pitch angle theta0 of the propeller [deg] (datasheet)
  blade_params['theta_1'] = blade_params['twist']*np.pi/180 # blade twist (measured)


def main():
  set_params()

  #integrator
  n_r =21
  n_th=20
  Integrator Solver(n_r, n_th)

  #Quadcopter
  Quadcopter Quad(Solver)

  #get state variables from CSV
  '''
  TODO
  '''
  
  #calculate prop forces
  for i in range(n):
    Quad.SetVelocity(state.getMotspeed(i),state.getVel(i), state.getAtt(i))
    state.F_hat[i] = Quad.GetForce()
    state.Tau_hat[i] = Quad.GetTorque()


  #append file




if __name__ == "__main__":
  main()
