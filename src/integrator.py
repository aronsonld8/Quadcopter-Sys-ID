class Integrator:
  def __init__(n_r, n_th, blade_params):
    self.n_r = n_r
    self.n_th = n_th

    #blade param constants
    self.R = blade_params['radius']
    self.K = blade_params['K']
    self.ci = blade_params['chord_inner']
    self.co = blade_params['chord_outer']
    self.cd = blade_params['drag_coefficient']
    self.cl = blade_params['lift_coefficient']
    self.theta_0 = blade_params['theta_0']
    self.theta_1 = blade_params['theta_1']

    #prop state variables
    self.a0  = 0
    self.a1s = 0
    self.b1s = 0
    self.omega = 0
    self.mu = 0
    self.Vz = 0
    self.vi = 0

    #get points,weights
    points_r, weights_r = np.polynomial.legendre.leggauss(n_r)
    self.th_points = np.linspace(0,2*np.pi,n_th, endpoint=False)
    self.th_weights = (2*np.pi)/(n_th)*np.ones(n_th)
    
    # Map r from [-1, 1] to [0, R]
    r_div2 = 0.5*blade_params['radius']
    self.r_points  = r_div2 * points_r + r_div2
    self.r_weights = r_div2 * weights_r


  def _diffElement(r, Psi, type='T'):
      r_ratio = r/self.R
      sPsi = np.sin(Psi)
      cPsi = np.cos(Psi)

      beta = self.a0 -  self.a1s*cPsi - self.b1s*sPsi
      U_T = self.omega * (r + self.R*self.mu*sPsi)
      corr = self.K*r_ratio*cPsi
      U_P = self.Vz - self.vi*(1+corr) - r*self.omega*(self.a1s*sPsi + self.b1s*cPsi) - self.Vz*self.beta*cPsi
      phi = np.atan2(U_P,U_T)
      f = 1.5 * (1 - r_ratio) / (r_ratio*np.abs(np.sin(phi)) + 1.0e-8)
      F = 2/np.pi * np.arccos(np.exp(-f))
      U2 = U_T**2 + U_P**2
      alpha = self.theta_0 + (self.theta_1 * r_ratio) + phi;
      _cl = self.cl*np.sin(alpha)*np.cos(alpha);
      _cd = self.cd * np.sin(alpha)**2;
      _c = self.ci + (r_ratio * (self.co - self.ci))
      #print(f"alpha = {alpha/np.pi}pi, cl = {_cl}, cd = {_cd}")
      dL = U2 * _cl * _c
      dD = U2 * _cd * _c
      temp = [0,0,0]
      if type == 'T':
        temp = [np.cos(phi), np.sin(phi), 1]
      elif type == 'H':
        temp = [-np.sin(phi), np.cos(phi), sPsi]
      elif type == 'Q':
        temp = [-np.sin(phi), np.cos(phi), r]
      return F*(temp[0]*dL + temp[1]*dD)*temp[2]

    
        
