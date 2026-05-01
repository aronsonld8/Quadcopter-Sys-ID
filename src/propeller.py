import parameters
from integrator import Integrator

#ref frame transforms
def FLU_to_FRD(vector_FLU):
  return [vector_FLU[0], -vector_FLU[1], -vector_FLU[2]]
def FRD_to_FLU(vector_FRD):
  return [vector_FRD[0], -vector_FRD[1], -vector_FRD[2]]

class Propeller:
    def __init__(self, r, dir, int_instance: Integrator):
        #variables
        self._solver = int_instance
        self.r   = r
        self.dir = dir
        self.vel = np.empty(3) #velocity relative to propeller in prop frame (forward, right, down)
        self.att = np.empty(3) #angular velocity relative to propeller in prop frame (forward, right, down)
        self.alpha_s=0
        self.state = {
            "a0": 0, #conning angle
            "a1s": 0, #lateral flapping angle
            "b1s": 0, #longitudinal flapping angle
            "omega":0,  #rotor angular velocity
            "mu" :0,  #adv ratio
            "Vz": 0,  #verical velocity
            "vi": 0  #induced velocity
        }
        self.Force = np.empty(3)
        self.Torque = np.empty(3)
   
  '''
  BEM calculation functions.
  reverse engineered from UZH NeuroBEM matlab and cpp code
  '''
    def _BEM_calc_a0(self):
        Omega = np.ones(7)
        mu = np.ones(7)
        for i in range(1,7):
            Omega[i] = self.state['omega']*Omega[i-1]
            mu[i] = self.state['mu']*mu[i-1]
        p = self.att[0]
        q = self.att[1]
        vind = self.state['vi']
        alpha = self.alpha_s
  #WARNING!! DO NOT TOUCH LOOK OR ATTEMPT TO UNDERSTAND THIS:
        self.state['a0'] = (17860.69768 - 9.085245727e-8 * Omega[3] * mu[1] * q +
           5.075429983e-14 * Omega[5] * mu[2] * vind -
           1.201149295e-15 * Omega[5] * mu[3] * p -
           3.287355978e-15 * Omega[6] * alpha * mu[3] -
           8.012930184e-15 * Omega[6] * alpha * mu[1] +
           1.011494156e-15 * Omega[6] * alpha * mu[5] -
           8.680285403 * Omega[1] * mu[1] * p -
           17.36057081 * Omega[2] * alpha * mu[1] +
           1.517241217e-15 * Omega[5] * mu[5] * p -
           5.958332683e-15 * Omega[5] * mu[1] * p -
           1.561670740e-14 * Omega[5] * mu[4] * vind -
           3.121900413e-12 * Omega[4] * mu[4] + 1.237136043e-13 * Omega[5] * vind -
           3.157428176e-15 * Omega[6] * mu[2] + 6.274197178e-16 * Omega[6] * mu[6] -
           2.031302206e-16 * Omega[6] * mu[4] + 268.0341334 * Omega[1] * vind -
           3.589529587 * Omega[2] * mu[2] - 1.902162223e-15 * Omega[6] +
           8.243768249e-12 * Omega[4] - 4.121166794 * Omega[2]) / (-4.377059028e8 - 1.147609987e-7 * Omega[4] * mu[2] +
           7.650732730e-8 * Omega[4] * mu[4] + 1.135865755e-14 * Omega[6] * mu[4] -
           2.999395522e-14 * Omega[6] - 2.020271598e-7 * Omega[4] -
           64.98399099 * Omega[2])

    def _BEM_calc_a1s(self):
        Omega = np.ones(7)
        mu = np.ones(7)
        for i in range(1,7):
            Omega[i] = self.state['omega']*Omega[i-1]
            mu[i] = self.state['mu']*mu[i-1]
        p = self.att[0]
        q = self.att[1]
        vind = self.state['vi']
        alpha = self.alpha_s
  #WARNING!! DO NOT TOUCH LOOK OR ATTEMPT TO UNDERSTAND THIS:
        self.state['a1s'] = 4.543463022e-14 * (0.4564908655 * Omega[5] * mu[3] + 0.7417976561 * Omega[5] * mu[1] +
           4.491002632e7 * Omega[2] * q + 0.9648437489 * Omega[4] * p +
           1.559168274e10 * Omega[1] * mu[1] - 58779.37391 * Omega[3] * mu[3] +
           1.398833618e6 * Omega[3] * mu[1] + 6.498797591e6 * Omega[2] * p +
           3.024957887e14 * q - 15.43924657 * Omega[4] * mu[3] * vind -
           3.578285830e6 * Omega[2] * mu[2] * p +
           6.499534432e7 * Omega[2] * mu[1] * vind -
           25.08877564 * Omega[4] * mu[1] * vind -
           1.039925259e8 * Omega[2] * mu[3] * vind + Omega[5] * alpha * mu[4] +
           1.625000001 * Omega[5] * alpha * mu[2] +
           0.5937500002 * Omega[4] * mu[2] * p +
           6.735595891e6 * Omega[3] * alpha * mu[4] -
           4.209748462e6 * Omega[3] * alpha * mu[2]) * Omega[1] / (-4.377059028e8 - 1.147609987e-7 * Omega[4] * mu[2] +
           7.650732720e-8 * Omega[4] * mu[4] + 1.135865755e-14 * Omega[6] * mu[4] -
           2.999395520e-14 * Omega[6] - 2.020271610e-7 * Omega[4] -
           64.98399099 * Omega[2])
      
    def _BEM_calc_b1s(self):
        Omega = np.ones(7)
        mu = np.ones(7)
        for i in range(1,7):
            Omega[i] = self.state['omega']*Omega[i-1]
            mu[i] = self.state['mu']*mu[i-1]
        p = self.att[0]
        q = self.att[1]
        vind = self.state['vi']
        alpha = self.alpha_s
      #WARNING!! DO NOT TOUCH LOOK OR ATTEMPT TO UNDERSTAND THIS:
        self.state['b1s'] = 3.034482349e-15 * Omega[1] * (14.44639117 * Omega[4] * q + 0.6202900305 * Omega[5] * mu[5] +
           1.177725817e17 * mu[1] * vind - 1.208793256 * Omega[5] * mu[3] -
           1.157259683 * Omega[5] * mu[1] + 9.730505319e7 * Omega[2] * q -
           3.482171723e15 * Omega[1] * mu[1] - 3086.424696 * Omega[3] * mu[3] -
           5.169754829e8 * Omega[3] * mu[1] - 6.724278493e8 * Omega[2] * p -
           5.988003267e7 * Omega[2] * mu[2] * q +
           1.500000041 * Omega[4] * mu[4] * p - 8.890086873 * Omega[4] * mu[2] * q -
           7.628130118e15 * Omega[1] * alpha * mu[2] - 4.529202255e15 * p -
           15.43924687 * Omega[4] * mu[3] * vind +
           1.748510208e10 * Omega[2] * mu[1] * vind +
           75.26632912 * Omega[4] * mu[1] * vind +
           0.9999999999 * Omega[5] * alpha * mu[4] -
           4.875000134 * Omega[5] * alpha * mu[2] -
           3.625000097 * Omega[4] * mu[2] * p -
           1.132510062e9 * Omega[3] * alpha * mu[2]) / (-4.377059028e8 - 1.147609987e-7 * Omega[4] * mu[2] +
           7.650732720e-8 * Omega[4] * mu[4] + 1.135865756e-14 * Omega[6] * mu[4] -
           2.999395520e-14 * Omega[6] - 2.020271601e-7 * Omega[4] -
           64.98399099 * Omega[2])


    
      

    def _BEM_T_calc(self):
        self._solver.update_state(self.state)
        return self._solver.getThrust()
    def _BEM_H_calc(self):
        return self._solver.getHforce()
    def _BEM_Q_calc(self):
        return self._solver.getQtorque()
    
    def _BEM_VRS_correction(self,Vz, Vh):
      k0=1;
      k1=-1.125
      k2=-1.372
      k3=-1.718
      k4=-0.655
      r = -Vz/Vh
      v2 = Vh*(k0 + k1*(r) + k2*(r)**2 + k3*(r)**3 + k4*(r)**4)
      return v2
    def _BEM_vi_costfunc(self, _vi):
      R = blade_params['radius']
      Vhor = self.state.mu*(self.state.omega*blade_params['radius'])
      Vver = self.vel[2] - _vi
      self.state=['vi'] = _v1
      T_BEM = BEM_T_calc()
      T_M =  2 * rho_air * A * _vi * np.sqrt(Vhor**2 + Vver**2)
      res = T_BEM - T_M
      return res
    def _BEM_calc_vi(self):
      v1 = fsolve(_BEM_vi_costfunc,vi_init,xtol=etol, maxfev=500)[0]
      v_ratio = self.vel[2]/v1
      if v_ratio >= 0 and v_ratio <= 2:
        vz = -self.vel[2]
        self.vel[2] = 0;
        vh=fsolve(_BEM_vi_costfunc,vi_init,xtol=etol, maxfev=500)[0]
        v2 = _BEM_VRS_correction(vz,vh)
        v1 = max(v1,v2)
        self.vel[2] = -vz
      self.state=['vi'] = v1
    def _update(self):
        _BEM_calc_vi()
        _BEM_calc_a0()
        _BEM_calc_a1s()
        _BEM_calc_b1s()
        
        T  = _BEM_T_calc()
        H  = _BEM_H_calc()
        Q  = _BEM_Q_calc()
  
        sgn2 = self.dir
        sgn1 = -self.dir
        self.Force = [
          -(H + np.sin(self.state['a1s'])*T),
          sgn1*np.sin(self.state['b1s'])*T,
          -T*np.cos(self.state['a0'])
        self.Torque = [ 
          sgn1*blade_params['k_beta']*self.state['b1s'],
          blade_params['k_beta']*self.state['a1s'],
          sgn2*Q]

    #public functions
    def SetVelocity(self, motOmega, Vel_CM, att_body):
        R = blade_params['radius']
        self.omega = motOmega
        self.vel = FLU_to_FRD(Vel_CM + np.cross(att_body,self.r))
        self.att = FLU_to_FRD(att_body)
        self.alpha_s = np.atan2(self.vel[2], Vhor)
        Vhor = np.sqrt(self.vel[0]**2 + self.vel[1]**2)

        #reset state
        self.state['a0']  = 0
        self.state['a1s'] = 0
        self.state['b1s'] = 0
        self.state['omega'] = motOmega
        self.state['mu']    = Vhor/(motOmega*R)
        self.state['Vz']    = self.vel[2]
        self.state['vi']    = 0
        _update()

    def GetForce(self):
        return FRD_to_FLU(self.Force)
    def GetTorque(self):
        return FRD_to_FLU(self.Torque)


