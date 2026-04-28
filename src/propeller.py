class Propeller:
    def __init__(self, r, dir):
        #variables
        self.r   = r
        self.dir = dir
        self.omega =0  #rotor angular velocity
        self.vi =0  #induced velocity
        self.vel = np.empty(3) #velocity relative to propeller in prop frame (forward, right, down)
        self.att = np.empty(3) #angular velocity relative to propeller in prop frame (forward, right, down)
        self.mu =0
        self.alpha_s=0
        self.state = {
            "T" : 0, #Thrust
            "H" : 0, #Hforce
            "Q" : 0, #Torque
            "a0": 0, #conning angle
            "a1s": 0, #lateral flapping angle
            "b1s": 0 #longitudinal flapping angle
        }
        self.Force = np.empty(3)
        self.Torque = np.empty(3)
   
  '''
  BEM calculation functions.
  reverse engineered from UZH NeuroBEM matlab and cpp code
  '''
  def _BEM_calc_a0(self):
    R = blade_params['radius']
    a = blade_params['a']
    c = blade_params['mean_chord']
    e = blade_params['e']
    th0 = blade_params['theta_0']
    th1 = blade_params['theta_1']
    k_beta = blade_params['k_beta']
    Mb = blade_params['Mb']
    Ib = blade_params['Ib']
    K = blade_params['K']
    Omega = self.omega
    v1 = self.vi
    mu = self.mu
    alpha_s = self.alpha_s
    pdot = self.att[0]
    qdot = self.att[1]
    #WARNING!! DO NOT TOUCH LOOK AT OR ATTEMPT TO UNDERSTAND THIS:
    a0 = -(((135*R**13*a**3*c**3*acc_g**2*mu**6*rho_air**3*th0 + 90*R**13*a**3*c**3*acc_g**2*mu**6*rho_air**3*th1 - 810*R**12*a**3*c**3*e*acc_g**2*mu**6*rho_air**3*th0 - 630*R**12*a**3*c**3*e*acc_g**2*mu**6*rho_air**3*th1 +
           2025*R**11*a**3*c**3*e**2*acc_g**2*mu**6*rho_air**3*th0 + 1890*R**11*a**3*c**3*e**2*acc_g**2*mu**6*rho_air**3*th1 - 2700*R**10*a**3*c**3*e**3*acc_g**2*mu**6*rho_air**3*th0 - 3150*R**10*a**3*c**3*e**3*acc_g**2*mu**6*rho_air**3*th1 +
           2025*R**9*a**3*c**3*e**4*acc_g**2*mu**6*rho_air**3*th0 + 3150*R**9*a**3*c**3*e**4*acc_g**2*mu**6*rho_air**3*th1 - 810*R**8*a**3*c**3*e**5*acc_g**2*mu**6*rho_air**3*th0 - 1890*R**8*a**3*c**3*e**5*acc_g**2*mu**6*rho_air**3*th1 +
           135*R**7*a**3*c**3*e**6*acc_g**2*mu**6*rho_air**3*th0 + 630*R**7*a**3*c**3*e**6*acc_g**2*mu**6*rho_air**3*th1 - 90*R**6*a**3*c**3*e**7*acc_g**2*mu**6*rho_air**3*th1 + 180*R**13*a**3*alpha_s*c**3*acc_g**2*mu**5*rho_air**3 -
           1530*R**12*a**3*alpha_s*c**3*e*acc_g**2*mu**5*rho_air**3 + 5400*R**11*a**3*alpha_s*c**3*e**2*acc_g**2*mu**5*rho_air**3 - 10350*R**10*a**3*alpha_s*c**3*e**3*acc_g**2*mu**5*rho_air**3 + 11700*R**9*a**3*alpha_s*c**3*e**4*acc_g**2*mu**5*rho_air**3 -
           7830*R**8*a**3*alpha_s*c**3*e**5*acc_g**2*mu**5*rho_air**3 + 2880*R**7*a**3*alpha_s*c**3*e**6*acc_g**2*mu**5*rho_air**3 - 450*R**6*a**3*alpha_s*c**3*e**7*acc_g**2*mu**5*rho_air**3 + 135*R**13*a**3*c**3*acc_g**2*mu**4*rho_air**3*th0 +
           108*R**13*a**3*c**3*acc_g**2*mu**4*rho_air**3*th1 - 1440*R**12*a**3*c**3*e*acc_g**2*mu**4*rho_air**3*th0 - 1242*R**12*a**3*c**3*e*acc_g**2*mu**4*rho_air**3*th1 + 5490*R**11*a**3*c**3*e**2*acc_g**2*mu**4*rho_air**3*th0 +
           5508*R**11*a**3*c**3*e**2*acc_g**2*mu**4*rho_air**3*th1 - 10260*R**10*a**3*c**3*e**3*acc_g**2*mu**4*rho_air**3*th0 - 12852*R**10*a**3*c**3*e**3*acc_g**2*mu**4*rho_air**3*th1 + 9900*R**9*a**3*c**3*e**4*acc_g**2*mu**4*rho_air**3*th0 +
           17388*R**9*a**3*c**3*e**4*acc_g**2*mu**4*rho_air**3*th1 - 3960*R**8*a**3*c**3*e**5*acc_g**2*mu**4*rho_air**3*th0 - 13608*R**8*a**3*c**3*e**5*acc_g**2*mu**4*rho_air**3*th1 - 810*R**7*a**3*c**3*e**6*acc_g**2*mu**4*rho_air**3*th0 +
           5292*R**7*a**3*c**3*e**6*acc_g**2*mu**4*rho_air**3*th1 + 1260*R**6*a**3*c**3*e**7*acc_g**2*mu**4*rho_air**3*th0 - 108*R**6*a**3*c**3*e**7*acc_g**2*mu**4*rho_air**3*th1 - 315*R**5*a**3*c**3*e**8*acc_g**2*mu**4*rho_air**3*th0 -
           648*R**5*a**3*c**3*e**8*acc_g**2*mu**4*rho_air**3*th1 + 162*R**4*a**3*c**3*e**9*acc_g**2*mu**4*rho_air**3*th1 - 1080*R**12*a**3*alpha_s*c**3*e*acc_g**2*mu**3*rho_air**3 + 7200*R**11*a**3*alpha_s*c**3*e**2*acc_g**2*mu**3*rho_air**3 -
           20160*R**10*a**3*alpha_s*c**3*e**3*acc_g**2*mu**3*rho_air**3 + 30240*R**9*a**3*alpha_s*c**3*e**4*acc_g**2*mu**3*rho_air**3 - 25200*R**8*a**3*alpha_s*c**3*e**5*acc_g**2*mu**3*rho_air**3 + 10080*R**7*a**3*alpha_s*c**3*e**6*acc_g**2*mu**3*rho_air**3 -
           1440*R**5*a**3*alpha_s*c**3*e**8*acc_g**2*mu**3*rho_air**3 + 360*R**4*a**3*alpha_s*c**3*e**9*acc_g**2*mu**3*rho_air**3 - 540*R**13*a**3*c**3*acc_g**2*mu**2*rho_air**3*th0 - 360*R**13*a**3*c**3*acc_g**2*mu**2*rho_air**3*th1 +
           2520*R**12*a**3*c**3*e*acc_g**2*mu**2*rho_air**3*th0 +  1920*R**12*a**3*c**3*e*acc_g**2*mu**2*rho_air**3*th1 - 3420*R**11*a**3*c**3*e**2*acc_g**2*mu**2*rho_air**3*th0 - 2920*R**11*a**3*c**3*e**2*acc_g**2*mu**2*rho_air**3*th1 -
           1440*R**10*a**3*c**3*e**3*acc_g**2*mu**2*rho_air**3*th0 - 2640*R**10*a**3*c**3*e**3*acc_g**2*mu**2*rho_air**3*th1 + 7560*R**9*a**3*c**3*e**4*acc_g**2*mu**2*rho_air**3*th0 + 14640*R**9*a**3*c**3*e**4*acc_g**2*mu**2*rho_air**3*th1 -
           5040*R**8*a**3*c**3*e**5*acc_g**2*mu**2*rho_air**3*th0 - 20160*R**8*a**3*c**3*e**5*acc_g**2*mu**2*rho_air**3*th1 - 2520*R**7*a**3*c**3*e**6*acc_g**2*mu**2*rho_air**3*th0 + 11760*R**7*a**3*c**3*e**6*acc_g**2*mu**2*rho_air**3*th1 +
           4320*R**6*a**3*c**3*e**7*acc_g**2*mu**2*rho_air**3*th0 - 480*R**6*a**3*c**3*e**7*acc_g**2*mu**2*rho_air**3*th1 - 1260*R**5*a**3*c**3*e**8*acc_g**2*mu**2*rho_air**3*th0 - 2760*R**5*a**3*c**3*e**8*acc_g**2*mu**2*rho_air**3*th1 -
           360*R**4*a**3*c**3*e**9*acc_g**2*mu**2*rho_air**3*th0 + 960*R**4*a**3*c**3*e**9*acc_g**2*mu**2*rho_air**3*th1 + 180*R**3*a**3*c**3*e**10*acc_g**2*mu**2*rho_air**3*th0 + 120*R**3*a**3*c**3*e**10*acc_g**2*mu**2*rho_air**3*th1 -
           80*R**2*a**3*c**3*e**11*acc_g**2*mu**2*rho_air**3*th1 - 720*R**13*a**3*alpha_s*c**3*acc_g**2*mu*rho_air**3 + 4920*R**12*a**3*alpha_s*c**3*e*acc_g**2*mu*rho_air**3 - 13760*R**11*a**3*alpha_s*c**3*e**2*acc_g**2*mu*rho_air**3 +
           19320*R**10*a**3*alpha_s*c**3*e**3*acc_g**2*mu*rho_air**3 - 12000*R**9*a**3*alpha_s*c**3*e**4*acc_g**2*mu*rho_air**3 - 1680*R**8*a**3*alpha_s*c**3*e**5*acc_g**2*mu*rho_air**3 + 6720*R**7*a**3*alpha_s*c**3*e**6*acc_g**2*mu*rho_air**3 -
           2640*R**6*a**3*alpha_s*c**3*e**7*acc_g**2*mu*rho_air**3 - 720*R**5*a**3*alpha_s*c**3*e**8*acc_g**2*mu*rho_air**3 + 600*R**4*a**3*alpha_s*c**3*e**9*acc_g**2*mu*rho_air**3 - 40*R**2*a**3*alpha_s*c**3*e**11*acc_g**2*mu*rho_air**3 -
           540*R**13*a**3*c**3*acc_g**2*rho_air**3*th0 - 432*R**13*a**3*c**3*acc_g**2*rho_air**3*th1 + 3600*R**12*a**3*c**3*e*acc_g**2*rho_air**3*th0 + 3384*R**12*a**3*c**3*e*acc_g**2*rho_air**3*th1 - 9840*R**11*a**3*c**3*e**2*acc_g**2*rho_air**3*th0 -
           11280*R**11*a**3*c**3*e**2*acc_g**2*rho_air**3*th1 + 13760*R**10*a**3*c**3*e**3*acc_g**2*rho_air**3*th0 + 20448*R**10*a**3*c**3*e**3*acc_g**2*rho_air**3*th1 - 9660*R**9*a**3*c**3*e**4*acc_g**2*rho_air**3*th0 - 20960*R**9*a**3*c**3*e**4*acc_g**2*rho_air**3*th1 +
           2400*R**8*a**3*c**3*e**5*acc_g**2*rho_air**3*th0 + 10584*R**8*a**3*c**3*e**5*acc_g**2*rho_air**3*th1 - 288*R**7*a**3*c**3*e**6*acc_g**2*rho_air**3*th1 + 960*R**6*a**3*c**3*e**7*acc_g**2*rho_air**3*th0 - 1920*R**6*a**3*c**3*e**7*acc_g**2*rho_air**3*th1 -
           660*R**5*a**3*c**3*e**8*acc_g**2*rho_air**3*th0 + 144*R**5*a**3*c**3*e**8*acc_g**2*rho_air**3*th1 - 240*R**4*a**3*c**3*e**9*acc_g**2*rho_air**3*th0 + 360*R**4*a**3*c**3*e**9*acc_g**2*rho_air**3*th1 + 240*R**3*a**3*c**3*e**10*acc_g**2*rho_air**3*th0 +
           48*R**3*a**3*c**3*e**10*acc_g**2*rho_air**3*th1 - 96*R**2*a**3*c**3*e**11*acc_g**2*rho_air**3*th1 - 20*R*a**3*c**3*e**12*acc_g**2*rho_air**3*th0 + 8*a**3*c**3*e**13*acc_g**2*rho_air**3*th1)*Omega**6 + (90*R**13*a**3*c**3*acc_g**2*mu**5*pdot*rho_air**3 -
           495*R**12*a**3*c**3*e*acc_g**2*mu**5*pdot*rho_air**3 + 1080*R**11*a**3*c**3*e**2*acc_g**2*mu**5*pdot*rho_air**3 - 1125*R**10*a**3*c**3*e**3*acc_g**2*mu**5*pdot*rho_air**3 + 450*R**9*a**3*c**3*e**4*acc_g**2*mu**5*pdot*rho_air**3 +
           135*R**8*a**3*c**3*e**5*acc_g**2*mu**5*pdot*rho_air**3 - 180*R**7*a**3*c**3*e**6*acc_g**2*mu**5*pdot*rho_air**3 + 45*R**6*a**3*c**3*e**7*acc_g**2*mu**5*pdot*rho_air**3 - 270*R**12*a**3*c**3*e*acc_g**2*mu**3*pdot*rho_air**3 -
           180*R**12*a**3*c**3*acc_g**2*mu**4*rho_air**3*v1 + 1440*R**11*a**3*c**3*e**2*acc_g**2*mu**3*pdot*rho_air**3 + 1530*R**11*a**3*c**3*e*acc_g**2*mu**4*rho_air**3*v1 - 3060*R**10*a**3*c**3*e**3*acc_g**2*mu**3*pdot*rho_air**3 -
           5400*R**10*a**3*c**3*e**2*acc_g**2*mu**4*rho_air**3*v1 + 3240*R**9*a**3*c**3*e**4*acc_g**2*mu**3*pdot*rho_air**3 + 10350*R**9*a**3*c**3*e**3*acc_g**2*mu**4*rho_air**3*v1 - 1800*R**8*a**3*c**3*e**5*acc_g**2*mu**3*pdot*rho_air**3 -
           11700*R**8*a**3*c**3*e**4*acc_g**2*mu**4*rho_air**3*v1 + 720*R**7*a**3*c**3*e**6*acc_g**2*mu**3*pdot*rho_air**3 + 7830*R**7*a**3*c**3*e**5*acc_g**2*mu**4*rho_air**3*v1 - 540*R**6*a**3*c**3*e**7*acc_g**2*mu**3*pdot*rho_air**3 -
           2880*R**6*a**3*c**3*e**6*acc_g**2*mu**4*rho_air**3*v1 + 360*R**5*a**3*c**3*e**8*acc_g**2*mu**3*pdot*rho_air**3 + 450*R**5*a**3*c**3*e**7*acc_g**2*mu**4*rho_air**3*v1 - 90*R**4*a**3*c**3*e**9*acc_g**2*mu**3*pdot*rho_air**3 -
           360*R**13*a**3*c**3*acc_g**2*mu*pdot*rho_air**3 + 1920*R**12*a**3*c**3*e*acc_g**2*mu*pdot*rho_air**3 - 3640*R**11*a**3*c**3*e**2*acc_g**2*mu*pdot*rho_air**3 + 1080*R**11*a**3*c**3*e*acc_g**2*mu**2*rho_air**3*v1 + 1800*R**10*a**3*c**3*e**3*acc_g**2*mu*pdot*rho_air**3 -
           7200*R**10*a**3*c**3*e**2*acc_g**2*mu**2*rho_air**3*v1 + 3600*R**9*a**3*c**3*e**4*acc_g**2*mu*pdot*rho_air**3 + 20160*R**9*a**3*c**3*e**3*acc_g**2*mu**2*rho_air**3*v1 - 6720*R**8*a**3*c**3*e**5*acc_g**2*mu*pdot*rho_air**3 - 30240*R**8*a**3*c**3*e**4*acc_g**2*mu**2*rho_air**3*v1 +
           5040*R**7*a**3*c**3*e**6*acc_g**2*mu*pdot*rho_air**3 + 25200*R**7*a**3*c**3*e**5*acc_g**2*mu**2*rho_air**3*v1 - 2160*R**6*a**3*c**3*e**7*acc_g**2*mu*pdot*rho_air**3 - 10080*R**6*a**3*c**3*e**6*acc_g**2*mu**2*rho_air**3*v1 + 600*R**5*a**3*c**3*e**8*acc_g**2*mu*pdot*rho_air**3 +
           1440*R**4*a**3*c**3*e**8*acc_g**2*mu**2*rho_air**3*v1 - 120*R**3*a**3*c**3*e**10*acc_g**2*mu*pdot*rho_air**3 - 360*R**3*a**3*c**3*e**9*acc_g**2*mu**2*rho_air**3*v1 + 40*R**2*a**3*c**3*e**11*acc_g**2*mu*pdot*rho_air**3 + 720*R**12*a**3*c**3*acc_g**2*rho_air**3*v1 -
           4920*R**11*a**3*c**3*e*acc_g**2*rho_air**3*v1 + 13760*R**10*a**3*c**3*e**2*acc_g**2*rho_air**3*v1 - 19320*R**9*a**3*c**3*e**3*acc_g**2*rho_air**3*v1 + 12000*R**8*a**3*c**3*e**4*acc_g**2*rho_air**3*v1 + 1680*R**7*a**3*c**3*e**5*acc_g**2*rho_air**3*v1 -
           6720*R**6*a**3*c**3*e**6*acc_g**2*rho_air**3*v1 + 2640*R**5*a**3*c**3*e**7*acc_g**2*rho_air**3*v1 + 720*R**4*a**3*c**3*e**8*acc_g**2*rho_air**3*v1 - 600*R**3*a**3*c**3*e**9*acc_g**2*rho_air**3*v1 + 40*R*a**3*c**3*e**11*acc_g**2*rho_air**3*v1)*Omega**5 +
           (-1080*Mb*R**9*a**2*c**2*acc_g**2*mu**4*rho_air**2 + 4320*Mb*R**8*a**2*c**2*e*acc_g**2*mu**4*rho_air**2 - 6480*Mb*R**7*a**2*c**2*e**2*acc_g**2*mu**4*rho_air**2 + 4320*Mb*R**6*a**2*c**2*e**3*acc_g**2*mu**4*rho_air**2 - 1080*Mb*R**5*a**2*c**2*e**4*acc_g**2*mu**4*rho_air**2 +
            4320*Mb*R**9*a**2*c**2*acc_g**2*rho_air**2 - 23040*Mb*R**8*a**2*c**2*e*acc_g**2*rho_air**2 + 48000*Mb*R**7*a**2*c**2*e**2*acc_g**2*rho_air**2 - 46080*Mb*R**6*a**2*c**2*e**3*acc_g**2*rho_air**2 + 14400*Mb*R**5*a**2*c**2*e**4*acc_g**2*rho_air**2 + 7680*Mb*R**4*a**2*c**2*e**5*acc_g**2*rho_air**2 -
            5760*Mb*R**3*a**2*c**2*e**6*acc_g**2*rho_air**2 + 480*Mb*R*a**2*c**2*e**8*acc_g**2*rho_air**2)*Omega**4 +
             (4320*K*Mb*R**7*a**2*c**2*e**2*acc_g*mu*rho_air**2*v1 - 20160*K*Mb*R**6*a**2*c**2*e**3*acc_g*mu*rho_air**2*v1 + 36000*K*Mb*R**5*a**2*c**2*e**4*acc_g*mu*rho_air**2*v1 - 28800*K*Mb*R**4*a**2*c**2*e**5*acc_g*mu*rho_air**2*v1 + 7200*K*Mb*R**3*a**2*c**2*e**6*acc_g*mu*rho_air**2*v1 +
              2880*K*Mb*R**2*a**2*c**2*e**7*acc_g*mu*rho_air**2*v1 - 1440*K*Mb*R*a**2*c**2*e**8*acc_g*mu*rho_air**2*v1 +
              4320*K*R**7*a**2*c**2*e*acc_g**2*k_beta*mu*rho_air**2*v1 - 20160*K*R**6*a**2*c**2*e**2*acc_g**2*k_beta*mu*rho_air**2*v1 + 36000*K*R**5*a**2*c**2*e**3*acc_g**2*k_beta*mu*rho_air**2*v1 - 28800*K*R**4*a**2*c**2*e**4*acc_g**2*k_beta*mu*rho_air**2*v1 +
              7200*K*R**3*a**2*c**2*e**5*acc_g**2*k_beta*mu*rho_air**2*v1 + 2880*K*R**2*a**2*c**2*e**6*acc_g**2*k_beta*mu*rho_air**2*v1 - 1440*K*R*a**2*c**2*e**7*acc_g**2*k_beta*mu*rho_air**2*v1 - 4320*Mb*R**8*a**2*c**2*e**2*acc_g*mu*qdot*rho_air**2 +
              14400*Mb*R**7*a**2*c**2*e**3*acc_g*mu*qdot*rho_air**2 - 15840*Mb*R**6*a**2*c**2*e**4*acc_g*mu*qdot*rho_air**2 + 5760*Mb*R**5*a**2*c**2*e**5*acc_g*mu*qdot*rho_air**2 - 1440*Mb*R**4*a**2*c**2*e**6*acc_g*mu*qdot*rho_air**2 + 2880*Mb*R**3*a**2*c**2*e**7*acc_g*mu*qdot*rho_air**2 -
              1440*Mb*R**2*a**2*c**2*e**8*acc_g*mu*qdot*rho_air**2 - 4320*R**8*a**2*c**2*e*acc_g**2*k_beta*mu*qdot*rho_air**2 + 14400*R**7*a**2*c**2*e**2*acc_g**2*k_beta*mu*qdot*rho_air**2 - 15840*R**6*a**2*c**2*e**3*acc_g**2*k_beta*mu*qdot*rho_air**2 +
              5760*R**5*a**2*c**2*e**4*acc_g**2*k_beta*mu*qdot*rho_air**2 - 1440*R**4*a**2*c**2*e**5*acc_g**2*k_beta*mu*qdot*rho_air**2 + 2880*R**3*a**2*c**2*e**6*acc_g**2*k_beta*mu*qdot*rho_air**2 - 1440*R**2*a**2*c**2*e**7*acc_g**2*k_beta*mu*qdot*rho_air**2)*Omega**3 +
              (-34560*Mb**2*R**5*a*c*e**2*mu**2*rho_air*th0 - 23040*Mb**2*R**5*a*c*e**2*mu**2*rho_air*th1 + 69120*Mb**2*R**4*a*c*e**3*mu**2*rho_air*th0 + 69120*Mb**2*R**4*a*c*e**3*mu**2*rho_air*th1 - 34560*Mb**2*R**3*a*c*e**4*mu**2*rho_air*th0 - 69120*Mb**2*R**3*a*c*e**4*mu**2*rho_air*th1 +
               23040*Mb**2*R**2*a*c*e**5*mu**2*rho_air*th1 - 69120*Mb*R**5*a*c*e*acc_g*k_beta*mu**2*rho_air*th0 - 46080*Mb*R**5*a*c*e*acc_g*k_beta*mu**2*rho_air*th1 + 138240*Mb*R**4*a*c*e**2*acc_g*k_beta*mu**2*rho_air*th0 + 138240*Mb*R**4*a*c*e**2*acc_g*k_beta*mu**2*rho_air*th1 -
               69120*Mb*R**3*a*c*e**3*acc_g*k_beta*mu**2*rho_air*th0 - 138240*Mb*R**3*a*c*e**3*acc_g*k_beta*mu**2*rho_air*th1 + 46080*Mb*R**2*a*c*e**4*acc_g*k_beta*mu**2*rho_air*th1 - 34560*R**5*a*c*acc_g**2*k_beta**2*mu**2*rho_air*th0 - 23040*R**5*a*c*acc_g**2*k_beta**2*mu**2*rho_air*th1 +
               69120*R**4*a*c*e*acc_g**2*k_beta**2*mu**2*rho_air*th0 + 69120*R**4*a*c*e*acc_g**2*k_beta**2*mu**2*rho_air*th1 - 34560*R**3*a*c*e**2*acc_g**2*k_beta**2*mu**2*rho_air*th0 - 69120*R**3*a*c*e**2*acc_g**2*k_beta**2*mu**2*rho_air*th1 + 23040*R**2*a*c*e**3*acc_g**2*k_beta**2*mu**2*rho_air*th1 -
               46080*Mb**2*R**5*a*alpha_s*c*e**2*mu*rho_air + 69120*Mb**2*R**4*a*alpha_s*c*e**3*mu*rho_air - 23040*Mb**2*R**2*a*alpha_s*c*e**5*mu*rho_air - 92160*Mb*R**5*a*alpha_s*c*e*acc_g*k_beta*mu*rho_air + 138240*Mb*R**4*a*alpha_s*c*e**2*acc_g*k_beta*mu*rho_air - 46080*Mb*R**2*a*alpha_s*c*e**4*acc_g*k_beta*mu*rho_air -
               46080*R**5*a*alpha_s*c*acc_g**2*k_beta**2*mu*rho_air + 69120*R**4*a*alpha_s*c*e*acc_g**2*k_beta**2*mu*rho_air - 23040*R**2*a*alpha_s*c*e**3*acc_g**2*k_beta**2*mu*rho_air - 34560*Mb**2*R**5*a*c*e**2*rho_air*th0 - 27648*Mb**2*R**5*a*c*e**2*rho_air*th1 + 46080*Mb**2*R**4*a*c*e**3*rho_air*th0 +
               69120*Mb**2*R**4*a*c*e**3*rho_air*th1 - 46080*Mb**2*R**3*a*c*e**4*rho_air*th1 - 11520*Mb**2*R*a*c*e**6*rho_air*th0 + 4608*Mb**2*a*c*e**7*rho_air*th1 - 69120*Mb*R**5*a*c*e*acc_g*k_beta*rho_air*th0 - 55296*Mb*R**5*a*c*e*acc_g*k_beta*rho_air*th1 + 92160*Mb*R**4*a*c*e**2*acc_g*k_beta*rho_air*th0 +
               138240*Mb*R**4*a*c*e**2*acc_g*k_beta*rho_air*th1 - 92160*Mb*R**3*a*c*e**3*acc_g*k_beta*rho_air*th1 - 23040*Mb*R*a*c*e**5*acc_g*k_beta*rho_air*th0 + 9216*Mb*a*c*e**6*acc_g*k_beta*rho_air*th1 - 34560*R**5*a*c*acc_g**2*k_beta**2*rho_air*th0 - 27648*R**5*a*c*acc_g**2*k_beta**2*rho_air*th1 +
               46080*R**4*a*c*e*acc_g**2*k_beta**2*rho_air*th0 + 69120*R**4*a*c*e*acc_g**2*k_beta**2*rho_air*th1 - 46080*R**3*a*c*e**2*acc_g**2*k_beta**2*rho_air*th1 - 11520*R*a*c*e**4*acc_g**2*k_beta**2*rho_air*th0 + 4608*a*c*e**5*acc_g**2*k_beta**2*rho_air*th1)*Omega**2 +
               (-23040*Mb**2*R**5*a*c*e**2*mu*pdot*rho_air + 34560*Mb**2*R**4*a*c*e**3*mu*pdot*rho_air - 11520*Mb**2*R**2*a*c*e**5*mu*pdot*rho_air - 46080*Mb*R**5*a*c*e*acc_g*k_beta*mu*pdot*rho_air + 69120*Mb*R**4*a*c*e**2*acc_g*k_beta*mu*pdot*rho_air -
                23040*Mb*R**2*a*c*e**4*acc_g*k_beta*mu*pdot*rho_air - 23040*R**5*a*c*acc_g**2*k_beta**2*mu*pdot*rho_air + 34560*R**4*a*c*e*acc_g**2*k_beta**2*mu*pdot*rho_air - 11520*R**2*a*c*e**3*acc_g**2*k_beta**2*mu*pdot*rho_air + 46080*Mb**2*R**4*a*c*e**2*rho_air*v1 -
                69120*Mb**2*R**3*a*c*e**3*rho_air*v1 + 23040*Mb**2*R*a*c*e**5*rho_air*v1 + 92160*Mb*R**4*a*c*e*acc_g*k_beta*rho_air*v1 - 138240*Mb*R**3*a*c*e**2*acc_g*k_beta*rho_air*v1 + 46080*Mb*R*a*c*e**4*acc_g*k_beta*rho_air*v1 + 46080*R**4*a*c*acc_g**2*k_beta**2*rho_air*v1 -
                69120*R**3*a*c*e*acc_g**2*k_beta**2*rho_air*v1 + 23040*R*a*c*e**3*acc_g**2*k_beta**2*rho_air*v1)*Omega + 276480*Mb**3*R*e**2 + 552960*Mb**2*R*e*acc_g*k_beta +
          276480*Mb*R*acc_g**2*k_beta**2)*acc_g )/(120*R*((-9*Ib *R**8*a**2*c**2*acc_g**3*mu**4*rho_air**2 + 36*Ib *R**7*a**2*c**2*e*acc_g**3*mu**4*rho_air**2 - 54*Ib *R**6*a**2*c**2*e**2*acc_g**3*mu**4*rho_air**2 + 36*Ib *R**5*a**2*c**2*e**3*acc_g**3*mu**4*rho_air**2 - 9*Ib *R**4*a**2*c**2*e**4*acc_g**3*mu**4*rho_air**2 -
                        9*Mb*R**8*a**2*c**2*e*acc_g**2*mu**4*rho_air**2 + 36*Mb*R**7*a**2*c**2*e**2*acc_g**2*mu**4*rho_air**2 - 54*Mb*R**6*a**2*c**2*e**3*acc_g**2*mu**4*rho_air**2 + 36*Mb*R**5*a**2*c**2*e**4*acc_g**2*mu**4*rho_air**2 - 9*Mb*R**4*a**2*c**2*e**5*acc_g**2*mu**4*rho_air**2 +
                        36*Ib *R**8*a**2*c**2*acc_g**3*rho_air**2 - 192*Ib *R**7*a**2*c**2*e*acc_g**3*rho_air**2 + 400*Ib *R**6*a**2*c**2*e**2*acc_g**3*rho_air**2 - 384*Ib *R**5*a**2*c**2*e**3*acc_g**3*rho_air**2 + 120*Ib *R**4*a**2*c**2*e**4*acc_g**3*rho_air**2 + 64*Ib *R**3*a**2*c**2*e**5*acc_g**3*rho_air**2 -
                        48*Ib *R**2*a**2*c**2*e**6*acc_g**3*rho_air**2 + 4*Ib *a**2*c**2*e**8*acc_g**3*rho_air**2 + 36*Mb*R**8*a**2*c**2*e*acc_g**2*rho_air**2 - 192*Mb*R**7*a**2*c**2*e**2*acc_g**2*rho_air**2 + 400*Mb*R**6*a**2*c**2*e**3*acc_g**2*rho_air**2 - 384*Mb*R**5*a**2*c**2*e**4*acc_g**2*rho_air**2 +
                        120*Mb*R**4*a**2*c**2*e**5*acc_g**2*rho_air**2 + 64*Mb*R**3*a**2*c**2*e**6*acc_g**2*rho_air**2 - 48*Mb*R**2*a**2*c**2*e**7*acc_g**2*rho_air**2 + 4*Mb*a**2*c**2*e**9*acc_g**2*rho_air**2)*Omega**6 +
                        (-9*R**8*a**2*c**2*acc_g**3*k_beta*mu**4*rho_air**2 + 36*R**7*a**2*c**2*e*acc_g**3*k_beta*mu**4*rho_air**2 - 54*R**6*a**2*c**2*e**2*acc_g**3*k_beta*mu**4*rho_air**2 + 36*R**5*a**2*c**2*e**3*acc_g**3*k_beta*mu**4*rho_air**2 - 9*R**4*a**2*c**2*e**4*acc_g**3*k_beta*mu**4*rho_air**2 +
                         48*Mb*R**7*a**2*c**2*e**2*acc_g**2*mu**2*rho_air**2 - 168*Mb*R**6*a**2*c**2*e**3*acc_g**2*mu**2*rho_air**2 + 192*Mb*R**5*a**2*c**2*e**4*acc_g**2*mu**2*rho_air**2 - 48*Mb*R**4*a**2*c**2*e**5*acc_g**2*mu**2*rho_air**2 - 48*Mb*R**3*a**2*c**2*e**6*acc_g**2*mu**2*rho_air**2 + 24*Mb*R**2*a**2*c**2*e**7*acc_g**2*mu**2*rho_air**2 +
                         48*R**7*a**2*c**2*e*acc_g**3*k_beta*mu**2*rho_air**2 - 168*R**6*a**2*c**2*e**2*acc_g**3*k_beta*mu**2*rho_air**2 + 192*R**5*a**2*c**2*e**3*acc_g**3*k_beta*mu**2*rho_air**2 - 48*R**4*a**2*c**2*e**4*acc_g**3*k_beta*mu**2*rho_air**2 - 48*R**3*a**2*c**2*e**5*acc_g**3*k_beta*mu**2*rho_air**2 + 24*R**2*a**2*c**2*e**6*acc_g**3*k_beta*mu**2*rho_air**2 +
                         36*R**8*a**2*c**2*acc_g**3*k_beta*rho_air**2 - 192*R**7*a**2*c**2*e*acc_g**3*k_beta*rho_air**2 + 400*R**6*a**2*c**2*e**2*acc_g**3*k_beta*rho_air**2 - 384*R**5*a**2*c**2*e**3*acc_g**3*k_beta*rho_air**2 + 120*R**4*a**2*c**2*e**4*acc_g**3*k_beta*rho_air**2 + 64*R**3*a**2*c**2*e**5*acc_g**3*k_beta*rho_air**2 - 48*R**2*a**2*c**2*e**6*acc_g**3*k_beta*rho_air**2 +
                         4*a**2*c**2*e**8*acc_g**3*k_beta*rho_air**2)*Omega**4 + (2304*Ib *Mb**2*e**2*acc_g + 4608*Ib *Mb*e*acc_g**2*k_beta + 2304*Ib *acc_g**3*k_beta**2 + 2304*Mb**3*e**3 + 4608*Mb**2*e**2*acc_g*k_beta + 2304*Mb*e*acc_g**2*k_beta**2)*Omega**2 + 2304*Mb**2*e**2*acc_g*k_beta + 4608*Mb*e*acc_g**2*k_beta**2 + 2304*acc_g**3*k_beta**3))
      self.state['a0'] = a0

    def _BEM_calc_a1s(self):
      R = blade_params['radius']
      a = blade_params['a']
      c = blade_params['mean_chord']
      e = blade_params['e']
      th0 = blade_params['theta_0']
      th1 = blade_params['theta_1']
      k_beta = blade_params['k_beta']
      Mb = blade_params['Mb']
      Ib = blade_params['Ib']
      K = blade_params['K']
      Omega = self.omega
      v1 = self.vi
      mu = self.mu
      alpha_s = self.alpha_s
      p = self.att[0]
      q = self.att[1]
        
      #WARNING!! DO NOT TOUCH LOOK AT OR ATTEMPT TO UNDERSTAND THIS:
      a1s = -2*((-90*Ib*R**9*a*alpha_s*c*acc_g**2*mu**4*rho_air + 360*Ib*R**8*a*alpha_s*c*e*acc_g**2*mu**4*rho_air - 540*Ib*R**7*a*alpha_s*c*e**2*acc_g**2*mu**4*rho_air +  360*Ib*R**6*a*alpha_s*c*e**3*acc_g**2*mu**4*rho_air - 90*Ib*R**5*a*alpha_s*c*e**4*acc_g**2*mu**4*rho_air - 90*Mb*R**9*a*alpha_s*c*e*acc_g*mu**4*rho_air +  360*Mb*R**8*a*alpha_s*c*e**2*acc_g*mu**4*rho_air - 540*Mb*R**7*a*alpha_s*c*e**3*acc_g*mu**4*rho_air +
             360*Mb*R**6*a*alpha_s*c*e**4*acc_g*mu**4*rho_air - 90*Mb*R**5*a*alpha_s*c*e**5*acc_g*mu**4*rho_air - 120*Ib*R**9*a*c*acc_g**2*mu**3*rho_air*th0 - 90*Ib*R**9*a*c*acc_g**2*mu**3*rho_air*th1 + 420*Ib*R**8*a*c*e*acc_g**2*mu**3*rho_air*th0 + 420*Ib*R**8*a*c*e*acc_g**2*mu**3*rho_air*th1 - 480*Ib*R**7*a*c*e**2*acc_g**2*mu**3*rho_air*th0 - 750*Ib*R**7*a*c*e**2*acc_g**2*mu**3*rho_air*th1 + 120*Ib*R**6*a*c*e**3*acc_g**2*mu**3*rho_air*th0 +
             600*Ib*R**6*a*c*e**3*acc_g**2*mu**3*rho_air*th1 + 120*Ib*R**5*a*c*e**4*acc_g**2*mu**3*rho_air*th0 - 150*Ib*R**5*a*c*e**4*acc_g**2*mu**3*rho_air*th1 - 60*Ib*R**4*a*c*e**5*acc_g**2*mu**3*rho_air*th0 - 60*Ib*R**4*a*c*e**5*acc_g**2*mu**3*rho_air*th1 + 30*Ib*R**3*a*c*e**6*acc_g**2*mu**3*rho_air*th1 - 120*Mb*R**9*a*c*e*acc_g*mu**3*rho_air*th0 - 90*Mb*R**9*a*c*e*acc_g*mu**3*rho_air*th1 + 420*Mb*R**8*a*c*e**2*acc_g*mu**3*rho_air*th0 +
             420*Mb*R**8*a*c*e**2*acc_g*mu**3*rho_air*th1 - 480*Mb*R**7*a*c*e**3*acc_g*mu**3*rho_air*th0 - 750*Mb*R**7*a*c*e**3*acc_g*mu**3*rho_air*th1 + 120*Mb*R**6*a*c*e**4*acc_g*mu**3*rho_air*th0 + 600*Mb*R**6*a*c*e**4*acc_g*mu**3*rho_air*th1 + 120*Mb*R**5*a*c*e**5*acc_g*mu**3*rho_air*th0 - 150*Mb*R**5*a*c*e**5*acc_g*mu**3*rho_air*th1 - 60*Mb*R**4*a*c*e**6*acc_g*mu**3*rho_air*th0 - 60*Mb*R**4*a*c*e**6*acc_g*mu**3*rho_air*th1 +
             30*Mb*R**3*a*c*e**7*acc_g*mu**3*rho_air*th1 - 180*Ib*R**9*a*alpha_s*c*acc_g**2*mu**2*rho_air + 840*Ib*R**8*a*alpha_s*c*e*acc_g**2*mu**2*rho_air - 1500*Ib*R**7*a*alpha_s*c*e**2*acc_g**2*mu**2*rho_air + 1200*Ib*R**6*a*alpha_s*c*e**3*acc_g**2*mu**2*rho_air - 300*Ib*R**5*a*alpha_s*c*e**4*acc_g**2*mu**2*rho_air - 120*Ib*R**4*a*alpha_s*c*e**5*acc_g**2*mu**2*rho_air + 60*Ib*R**3*a*alpha_s*c*e**6*acc_g**2*mu**2*rho_air -
             180*Mb*R**9*a*alpha_s*c*e*acc_g*mu**2*rho_air + 840*Mb*R**8*a*alpha_s*c*e**2*acc_g*mu**2*rho_air - 1500*Mb*R**7*a*alpha_s*c*e**3*acc_g*mu**2*rho_air + 1200*Mb*R**6*a*alpha_s*c*e**4*acc_g*mu**2*rho_air - 300*Mb*R**5*a*alpha_s*c*e**5*acc_g*mu**2*rho_air - 120*Mb*R**4*a*alpha_s*c*e**6*acc_g*mu**2*rho_air + 60*Mb*R**3*a*alpha_s*c*e**7*acc_g*mu**2*rho_air - 240*Ib*R**9*a*c*acc_g**2*mu*rho_air*th0 -
             180*Ib*R**9*a*c*acc_g**2*mu*rho_air*th1 + 1000*Ib*R**8*a*c*e*acc_g**2*mu*rho_air*th0 + 960*Ib*R**8*a*c*e*acc_g**2*mu*rho_air*th1 - 1440*Ib*R**7*a*c*e**2*acc_g**2*mu*rho_air*th0 - 2000*Ib*R**7*a*c*e**2*acc_g**2*mu*rho_air*th1 + 600*Ib*R**6*a*c*e**3*acc_g**2*mu*rho_air*th0 + 1920*Ib*R**6*a*c*e**3*acc_g**2*mu*rho_air*th1 + 400*Ib*R**5*a*c*e**4*acc_g**2*mu*rho_air*th0 - 600*Ib*R**5*a*c*e**4*acc_g**2*mu*rho_air*th1 -
             360*Ib*R**4*a*c*e**5*acc_g**2*mu*rho_air*th0 - 320*Ib*R**4*a*c*e**5*acc_g**2*mu*rho_air*th1 + 240*Ib*R**3*a*c*e**6*acc_g**2*mu*rho_air*th1 + 40*Ib*R**2*a*c*e**7*acc_g**2*mu*rho_air*th0 - 20*Ib*R*a*c*e**8*acc_g**2*mu*rho_air*th1 - 240*Mb*R**9*a*c*e*acc_g*mu*rho_air*th0 - 180*Mb*R**9*a*c*e*acc_g*mu*rho_air*th1 + 1000*Mb*R**8*a*c*e**2*acc_g*mu*rho_air*th0 + 960*Mb*R**8*a*c*e**2*acc_g*mu*rho_air*th1 -
             1440*Mb*R**7*a*c*e**3*acc_g*mu*rho_air*th0 - 2000*Mb*R**7*a*c*e**3*acc_g*mu*rho_air*th1 + 600*Mb*R**6*a*c*e**4*acc_g*mu*rho_air*th0 + 1920*Mb*R**6*a*c*e**4*acc_g*mu*rho_air*th1 + 400*Mb*R**5*a*c*e**5*acc_g*mu*rho_air*th0 - 600*Mb*R**5*a*c*e**5*acc_g*mu*rho_air*th1 - 360*Mb*R**4*a*c*e**6*acc_g*mu*rho_air*th0 - 320*Mb*R**4*a*c*e**6*acc_g*mu*rho_air*th1 + 240*Mb*R**3*a*c*e**7*acc_g*mu*rho_air*th1 +
             40*Mb*R**2*a*c*e**8*acc_g*mu*rho_air*th0 - 20*Mb*R*a*c*e**9*acc_g*mu*rho_air*th1)*Omega**5 +
             (-45*Ib*R**9*a*c*acc_g**2*mu**2*p*rho_air + 150*Ib*R**8*a*c*e*acc_g**2*mu**2*p*rho_air + 90*Ib*R**8*a*c*acc_g**2*mu**3*rho_air*v1 - 165*Ib*R**7*a*c*e**2*acc_g**2*mu**2*p*rho_air - 360*Ib*R**7*a*c*e*acc_g**2*mu**3*rho_air*v1 + 60*Ib*R**6*a*c*e**3*acc_g**2*mu**2*p*rho_air + 540*Ib*R**6*a*c*e**2*acc_g**2*mu**3*rho_air*v1 - 15*Ib*R**5*a*c*e**4*acc_g**2*mu**2*p*rho_air - 360*Ib*R**5*a*c*e**3*acc_g**2*mu**3*rho_air*v1 +
              30*Ib*R**4*a*c*e**5*acc_g**2*mu**2*p*rho_air + 90*Ib*R**4*a*c*e**4*acc_g**2*mu**3*rho_air*v1 - 15*Ib*R**3*a*c*e**6*acc_g**2*mu**2*p*rho_air - 45*Mb*R**9*a*c*e*acc_g*mu**2*p*rho_air + 150*Mb*R**8*a*c*e**2*acc_g*mu**2*p*rho_air + 90*Mb*R**8*a*c*e*acc_g*mu**3*rho_air*v1 - 165*Mb*R**7*a*c*e**3*acc_g*mu**2*p*rho_air - 360*Mb*R**7*a*c*e**2*acc_g*mu**3*rho_air*v1 + 60*Mb*R**6*a*c*e**4*acc_g*mu**2*p*rho_air +
              540*Mb*R**6*a*c*e**3*acc_g*mu**3*rho_air*v1 - 15*Mb*R**5*a*c*e**5*acc_g*mu**2*p*rho_air - 360*Mb*R**5*a*c*e**4*acc_g*mu**3*rho_air*v1 + 30*Mb*R**4*a*c*e**6*acc_g*mu**2*p*rho_air + 90*Mb*R**4*a*c*e**5*acc_g*mu**3*rho_air*v1 - 15*Mb*R**3*a*c*e**7*acc_g*mu**2*p*rho_air - 90*Ib*R**9*a*c*acc_g**2*p*rho_air + 360*Ib*R**8*a*c*e*acc_g**2*p*rho_air + 180*Ib*R**8*a*c*acc_g**2*mu*rho_air*v1 - 500*Ib*R**7*a*c*e**2*acc_g**2*p*rho_air -
              840*Ib*R**7*a*c*e*acc_g**2*mu*rho_air*v1 + 240*Ib*R**6*a*c*e**3*acc_g**2*p*rho_air + 1500*Ib*R**6*a*c*e**2*acc_g**2*mu*rho_air*v1 - 1200*Ib*R**5*a*c*e**3*acc_g**2*mu*rho_air*v1 + 40*Ib*R**4*a*c*e**5*acc_g**2*p*rho_air + 300*Ib*R**4*a*c*e**4*acc_g**2*mu*rho_air*v1 - 60*Ib*R**3*a*c*e**6*acc_g**2*p*rho_air + 120*Ib*R**3*a*c*e**5*acc_g**2*mu*rho_air*v1 - 60*Ib*R**2*a*c*e**6*acc_g**2*mu*rho_air*v1 +
              10*Ib*R*a*c*e**8*acc_g**2*p*rho_air - 90*Mb*R**9*a*c*e*acc_g*p*rho_air + 360*Mb*R**8*a*c*e**2*acc_g*p*rho_air + 180*Mb*R**8*a*c*e*acc_g*mu*rho_air*v1 - 500*Mb*R**7*a*c*e**3*acc_g*p*rho_air - 840*Mb*R**7*a*c*e**2*acc_g*mu*rho_air*v1 + 240*Mb*R**6*a*c*e**4*acc_g*p*rho_air + 1500*Mb*R**6*a*c*e**3*acc_g*mu*rho_air*v1 - 1200*Mb*R**5*a*c*e**4*acc_g*mu*rho_air*v1 + 40*Mb*R**4*a*c*e**6*acc_g*p*rho_air +
              300*Mb*R**4*a*c*e**5*acc_g*mu*rho_air*v1 - 60*Mb*R**3*a*c*e**7*acc_g*p*rho_air + 120*Mb*R**3*a*c*e**6*acc_g*mu*rho_air*v1 - 60*Mb*R**2*a*c*e**7*acc_g*mu*rho_air*v1 + 10*Mb*R*a*c*e**9*acc_g*p*rho_air)*Omega**4 +
              (-90*R**9*a*alpha_s*c*acc_g**2*k_beta*mu**4*rho_air + 360*R**8*a*alpha_s*c*e*acc_g**2*k_beta*mu**4*rho_air - 540*R**7*a*alpha_s*c*e**2*acc_g**2*k_beta*mu**4*rho_air + 360*R**6*a*alpha_s*c*e**3*acc_g**2*k_beta*mu**4*rho_air - 90*R**5*a*alpha_s*c*e**4*acc_g**2*k_beta*mu**4*rho_air + 120*Mb*R**9*a*c*e*acc_g*mu**3*rho_air*th0 + 80*Mb*R**9*a*c*e*acc_g*mu**3*rho_air*th1 - 420*Mb*R**8*a*c*e**2*acc_g*mu**3*rho_air*th0 -
               360*Mb*R**8*a*c*e**2*acc_g*mu**3*rho_air*th1 + 480*Mb*R**7*a*c*e**3*acc_g*mu**3*rho_air*th0 + 600*Mb*R**7*a*c*e**3*acc_g*mu**3*rho_air*th1 - 120*Mb*R**6*a*c*e**4*acc_g*mu**3*rho_air*th0 - 400*Mb*R**6*a*c*e**4*acc_g*mu**3*rho_air*th1 - 120*Mb*R**5*a*c*e**5*acc_g*mu**3*rho_air*th0 + 60*Mb*R**4*a*c*e**6*acc_g*mu**3*rho_air*th0 + 120*Mb*R**4*a*c*e**6*acc_g*mu**3*rho_air*th1 - 40*Mb*R**3*a*c*e**7*acc_g*mu**3*rho_air*th1 -
               10*R**9*a*c*acc_g**2*k_beta*mu**3*rho_air*th1 + 60*R**8*a*c*e*acc_g**2*k_beta*mu**3*rho_air*th1 - 150*R**7*a*c*e**2*acc_g**2*k_beta*mu**3*rho_air*th1 + 200*R**6*a*c*e**3*acc_g**2*k_beta*mu**3*rho_air*th1 - 150*R**5*a*c*e**4*acc_g**2*k_beta*mu**3*rho_air*th1 + 60*R**4*a*c*e**5*acc_g**2*k_beta*mu**3*rho_air*th1 - 10*R**3*a*c*e**6*acc_g**2*k_beta*mu**3*rho_air*th1 + 160*Mb*R**9*a*alpha_s*c*e*acc_g*mu**2*rho_air -
               480*Mb*R**8*a*alpha_s*c*e**2*acc_g*mu**2*rho_air + 360*Mb*R**7*a*alpha_s*c*e**3*acc_g*mu**2*rho_air + 160*Mb*R**6*a*alpha_s*c*e**4*acc_g*mu**2*rho_air - 240*Mb*R**5*a*alpha_s*c*e**5*acc_g*mu**2*rho_air + 40*Mb*R**3*a*alpha_s*c*e**7*acc_g*mu**2*rho_air - 20*R**9*a*alpha_s*c*acc_g**2*k_beta*mu**2*rho_air + 360*R**8*a*alpha_s*c*e*acc_g**2*k_beta*mu**2*rho_air - 1140*R**7*a*alpha_s*c*e**2*acc_g**2*k_beta*mu**2*rho_air +
               1360*R**6*a*alpha_s*c*e**3*acc_g**2*k_beta*mu**2*rho_air - 540*R**5*a*alpha_s*c*e**4*acc_g**2*k_beta*mu**2*rho_air - 120*R**4*a*alpha_s*c*e**5*acc_g**2*k_beta*mu**2*rho_air + 100*R**3*a*alpha_s*c*e**6*acc_g**2*k_beta*mu**2*rho_air + 120*Mb*R**9*a*c*e*acc_g*mu*rho_air*th0 + 96*Mb*R**9*a*c*e*acc_g*mu*rho_air*th1 - 340*Mb*R**8*a*c*e**2*acc_g*mu*rho_air*th0 - 384*Mb*R**8*a*c*e**2*acc_g*mu*rho_air*th1 +
               240*Mb*R**7*a*c*e**3*acc_g*mu*rho_air*th0 + 520*Mb*R**7*a*c*e**3*acc_g*mu*rho_air*th1 + 60*Mb*R**6*a*c*e**4*acc_g*mu*rho_air*th0 - 192*Mb*R**6*a*c*e**4*acc_g*mu*rho_air*th1 - 40*Mb*R**5*a*c*e**5*acc_g*mu*rho_air*th0 - 120*Mb*R**5*a*c*e**5*acc_g*mu*rho_air*th1 - 60*Mb*R**4*a*c*e**6*acc_g*mu*rho_air*th0 + 64*Mb*R**4*a*c*e**6*acc_g*mu*rho_air*th1 + 24*Mb*R**3*a*c*e**7*acc_g*mu*rho_air*th1 +
               20*Mb*R**2*a*c*e**8*acc_g*mu*rho_air*th0 - 8*Mb*R*a*c*e**9*acc_g*mu*rho_air*th1 - 120*R**9*a*c*acc_g**2*k_beta*mu*rho_air*th0 - 84*R**9*a*c*acc_g**2*k_beta*mu*rho_air*th1 + 660*R**8*a*c*e*acc_g**2*k_beta*mu*rho_air*th0 + 576*R**8*a*c*e*acc_g**2*k_beta*mu*rho_air*th1 - 1200*R**7*a*c*e**2*acc_g**2*k_beta*mu*rho_air*th0 - 1480*R**7*a*c*e**2*acc_g**2*k_beta*mu*rho_air*th1 +
               660*R**6*a*c*e**3*acc_g**2*k_beta*mu*rho_air*th0 + 1728*R**6*a*c*e**3*acc_g**2*k_beta*mu*rho_air*th1 + 360*R**5*a*c*e**4*acc_g**2*k_beta*mu*rho_air*th0 - 720*R**5*a*c*e**4*acc_g**2*k_beta*mu*rho_air*th1 - 420*R**4*a*c*e**5*acc_g**2*k_beta*mu*rho_air*th0 - 256*R**4*a*c*e**5*acc_g**2*k_beta*mu*rho_air*th1 + 264*R**3*a*c*e**6*acc_g**2*k_beta*mu*rho_air*th1 + 60*R**2*a*c*e**7*acc_g**2*k_beta*mu*rho_air*th0 - 28*R*a*c*e**8*acc_g**2*k_beta*mu*rho_air*th1)*Omega**3 +
              (80*Mb*R**9*a*c*e*acc_g*mu**2*p*rho_air - 240*Mb*R**8*a*c*e**2*acc_g*mu**2*p*rho_air + 180*Mb*R**7*a*c*e**3*acc_g*mu**2*p*rho_air + 80*Mb*R**6*a*c*e**4*acc_g*mu**2*p*rho_air - 120*Mb*R**5*a*c*e**5*acc_g*mu**2*p*rho_air + 20*Mb*R**3*a*c*e**7*acc_g*mu**2*p*rho_air + 35*R**9*a*c*acc_g**2*k_beta*mu**2*p*rho_air - 90*R**8*a*c*e*acc_g**2*k_beta*mu**2*p*rho_air + 90*R**8*a*c*acc_g**2*k_beta*mu**3*rho_air*v1 + 15*R**7*a*c*e**2*acc_g**2*k_beta*mu**2*p*rho_air -
               360*R**7*a*c*e*acc_g**2*k_beta*mu**3*rho_air*v1 +  40*R**6*a*c*e**3*acc_g**2*k_beta*mu**2*p*rho_air + 540*R**6*a*c*e**2*acc_g**2*k_beta*mu**3*rho_air*v1 - 135*R**5*a*c*e**4*acc_g**2*k_beta*mu**2*p*rho_air - 360*R**5*a*c*e**3*acc_g**2*k_beta*mu**3*rho_air*v1 + 30*R**4*a*c*e**5*acc_g**2*k_beta*mu**2*p*rho_air + 90*R**4*a*c*e**4*acc_g**2*k_beta*mu**3*rho_air*v1 + 5*R**3*a*c*e**6*acc_g**2*k_beta*mu**2*p*rho_air - 160*Mb*R**8*a*c*e*acc_g*mu*rho_air*v1 +
               480*Mb*R**7*a*c*e**2*acc_g*mu*rho_air*v1 - 360*Mb*R**6*a*c*e**3*acc_g*mu*rho_air*v1 - 160*Mb*R**5*a*c*e**4*acc_g*mu*rho_air*v1 + 240*Mb*R**4*a*c*e**5*acc_g*mu*rho_air*v1 - 40*Mb*R**2*a*c*e**7*acc_g*mu*rho_air*v1 - 90*R**9*a*c*acc_g**2*k_beta*p*rho_air + 360*R**8*a*c*e*acc_g**2*k_beta*p*rho_air + 20*R**8*a*c*acc_g**2*k_beta*mu*rho_air*v1 - 500*R**7*a*c*e**2*acc_g**2*k_beta*p*rho_air - 360*R**7*a*c*e*acc_g**2*k_beta*mu*rho_air*v1 +
               240*R**6*a*c*e**3*acc_g**2*k_beta*p*rho_air + 1140*R**6*a*c*e**2*acc_g**2*k_beta*mu*rho_air*v1 - 1360*R**5*a*c*e**3*acc_g**2*k_beta*mu*rho_air*v1 + 40*R**4*a*c*e**5*acc_g**2*k_beta*p*rho_air + 540*R**4*a*c*e**4*acc_g**2*k_beta*mu*rho_air*v1 - 60*R**3*a*c*e**6*acc_g**2*k_beta*p*rho_air + 120*R**3*a*c*e**5*acc_g**2*k_beta*mu*rho_air*v1 - 100*R**2*a*c*e**6*acc_g**2*k_beta*mu*rho_air*v1 + 10*R*a*c*e**8*acc_g**2*k_beta*p*rho_air + 720*Ib*K*Mb*R**4*e*acc_g*v1 -
               1920*Ib*K*Mb*R**3*e**2*acc_g*v1 + 1440*Ib*K*Mb*R**2*e**3*acc_g*v1 - 240*Ib*K*Mb*e**5*acc_g*v1 + 720*Ib*K*R**4*acc_g**2*k_beta*v1 - 1920*Ib*K*R**3*e*acc_g**2*k_beta*v1 + 1440*Ib*K*R**2*e**2*acc_g**2*k_beta*v1 - 240*Ib*K*e**4*acc_g**2*k_beta*v1 - 720*Ib*Mb*R**5*e*acc_g*q + 960*Ib*Mb*R**4*e**2*acc_g*q - 240*Ib*Mb*R*e**5*acc_g*q - 720*Ib*R**5*acc_g**2*k_beta*q + 960*Ib*R**4*e*acc_g**2*k_beta*q - 240*Ib*R*e**4*acc_g**2*k_beta*q + 720*K*Mb**2*R**4*e**2*v1 -
               1920*K*Mb**2*R**3*e**3*v1 + 1440*K*Mb**2*R**2*e**4*v1 - 240*K*Mb**2*e**6*v1 + 720*K*Mb*R**4*e*acc_g*k_beta*v1 - 1920*K*Mb*R**3*e**2*acc_g*k_beta*v1 + 1440*K*Mb*R**2*e**3*acc_g*k_beta*v1 - 240*K*Mb*e**5*acc_g*k_beta*v1 - 720*Mb**2*R**5*e**2*q + 960*Mb**2*R**4*e**3*q - 240*Mb**2*R*e**6*q - 720*Mb*R**5*e*acc_g*k_beta*q + 960*Mb*R**4*e**2*acc_g*k_beta*q - 240*Mb*R*e**5*acc_g*k_beta*q)*Omega**2 +
               (-960*Mb**2*R**5*e*acc_g*mu + 1440*Mb**2*R**4*e**2*acc_g*mu - 480*Mb**2*R**2*e**4*acc_g*mu - 960*Mb*R**5*acc_g**2*k_beta*mu + 1440*Mb*R**4*e*acc_g**2*k_beta*mu - 480*Mb*R**2*e**3*acc_g**2*k_beta*mu)*Omega + 720*K*Mb*R**4*e*acc_g*k_beta*v1 - 1920*K*Mb*R**3*e**2*acc_g*k_beta*v1 + 1440*K*Mb*R**2*e**3*acc_g*k_beta*v1 - 240*K*Mb*e**5*acc_g*k_beta*v1 + 720*K*R**4*acc_g**2*k_beta**2*v1 - 1920*K*R**3*e*acc_g**2*k_beta**2*v1 + 1440*K*R**2*e**2*acc_g**2*k_beta**2*v1 -
            240*K*e**4*acc_g**2*k_beta**2*v1 - 720*Mb*R**5*e*acc_g*k_beta*q + 960*Mb*R**4*e**2*acc_g*k_beta*q - 240*Mb*R*e**5*acc_g*k_beta*q - 720*R**5*acc_g**2*k_beta**2*q + 960*R**4*e*acc_g**2*k_beta**2*q -
            240*R*e**4*acc_g**2*k_beta**2*q)*Omega*a*c*acc_g*rho_air/(5*R*((-9*Ib*R**8*a**2*c**2*acc_g**3*mu**4*rho_air**2 + 36*Ib*R**7*a**2*c**2*e*acc_g**3*mu**4*rho_air**2 - 54*Ib*R**6*a**2*c**2*e**2*acc_g**3*mu**4*rho_air**2 + 36*Ib*R**5*a**2*c**2*e**3*acc_g**3*mu**4*rho_air**2 - 9*Ib*R**4*a**2*c**2*e**4*acc_g**3*mu**4*rho_air**2 - 9*Mb*R**8*a**2*c**2*e*acc_g**2*mu**4*rho_air**2 + 36*Mb*R**7*a**2*c**2*e**2*acc_g**2*mu**4*rho_air**2 - 54*Mb*R**6*a**2*c**2*e**3*acc_g**2*mu**4*rho_air**2 + 36*Mb*R**5*a**2*c**2*e**4*acc_g**2*mu**4*rho_air**2 -
                    9*Mb*R**4*a**2*c**2*e**5*acc_g**2*mu**4*rho_air**2 + 36*Ib*R**8*a**2*c**2*acc_g**3*rho_air**2 - 192*Ib*R**7*a**2*c**2*e*acc_g**3*rho_air**2 + 400*Ib*R**6*a**2*c**2*e**2*acc_g**3*rho_air**2 - 384*Ib*R**5*a**2*c**2*e**3*acc_g**3*rho_air**2 + 120*Ib*R**4*a**2*c**2*e**4*acc_g**3*rho_air**2 + 64*Ib*R**3*a**2*c**2*e**5*acc_g**3*rho_air**2 - 48*Ib*R**2*a**2*c**2*e**6*acc_g**3*rho_air**2 + 4*Ib*a**2*c**2*e**8*acc_g**3*rho_air**2 + 36*Mb*R**8*a**2*c**2*e*acc_g**2*rho_air**2 -
                    192*Mb*R**7*a**2*c**2*e**2*acc_g**2*rho_air**2 + 400*Mb*R**6*a**2*c**2*e**3*acc_g**2*rho_air**2 - 384*Mb*R**5*a**2*c**2*e**4*acc_g**2*rho_air**2 + 120*Mb*R**4*a**2*c**2*e**5*acc_g**2*rho_air**2 + 64*Mb*R**3*a**2*c**2*e**6*acc_g**2*rho_air**2 - 48*Mb*R**2*a**2*c**2*e**7*acc_g**2*rho_air**2 + 4*Mb*a**2*c**2*e**9*acc_g**2*rho_air**2)*Omega**6 +
                    (-9*R**8*a**2*c**2*acc_g**3*k_beta*mu**4*rho_air**2 + 36*R**7*a**2*c**2*e*acc_g**3*k_beta*mu**4*rho_air**2 - 54*R**6*a**2*c**2*e**2*acc_g**3*k_beta*mu**4*rho_air**2 + 36*R**5*a**2*c**2*e**3*acc_g**3*k_beta*mu**4*rho_air**2 - 9*R**4*a**2*c**2*e**4*acc_g**3*k_beta*mu**4*rho_air**2 + 48*Mb*R**7*a**2*c**2*e**2*acc_g**2*mu**2*rho_air**2 - 168*Mb*R**6*a**2*c**2*e**3*acc_g**2*mu**2*rho_air**2 + 192*Mb*R**5*a**2*c**2*e**4*acc_g**2*mu**2*rho_air**2 -
                    48*Mb*R**4*a**2*c**2*e**5*acc_g**2*mu**2*rho_air**2 - 48*Mb*R**3*a**2*c**2*e**6*acc_g**2*mu**2*rho_air**2 + 24*Mb*R**2*a**2*c**2*e**7*acc_g**2*mu**2*rho_air**2 + 48*R**7*a**2*c**2*e*acc_g**3*k_beta*mu**2*rho_air**2 - 168*R**6*a**2*c**2*e**2*acc_g**3*k_beta*mu**2*rho_air**2 + 192*R**5*a**2*c**2*e**3*acc_g**3*k_beta*mu**2*rho_air**2 - 48*R**4*a**2*c**2*e**4*acc_g**3*k_beta*mu**2*rho_air**2 - 48*R**3*a**2*c**2*e**5*acc_g**3*k_beta*mu**2*rho_air**2 +
                    24*R**2*a**2*c**2*e**6*acc_g**3*k_beta*mu**2*rho_air**2 + 36*R**8*a**2*c**2*acc_g**3*k_beta*rho_air**2 - 192*R**7*a**2*c**2*e*acc_g**3*k_beta*rho_air**2 + 400*R**6*a**2*c**2*e**2*acc_g**3*k_beta*rho_air**2 - 384*R**5*a**2*c**2*e**3*acc_g**3*k_beta*rho_air**2 + 120*R**4*a**2*c**2*e**4*acc_g**3*k_beta*rho_air**2 + 64*R**3*a**2*c**2*e**5*acc_g**3*k_beta*rho_air**2 - 48*R**2*a**2*c**2*e**6*acc_g**3*k_beta*rho_air**2 + 4*a**2*c**2*e**8*acc_g**3*k_beta*rho_air**2)*Omega**4 +
                    (2304*Ib*Mb**2*e**2*acc_g + 4608*Ib*Mb*e*acc_g**2*k_beta + 2304*Ib*acc_g**3*k_beta**2 + 2304*Mb**3*e**3 + 4608*Mb**2*e**2*acc_g*k_beta + 2304*Mb*e*acc_g**2*k_beta**2)*Omega**2 + 2304*Mb**2*e**2*acc_g*k_beta + 4608*Mb*e*acc_g**2*k_beta**2 + 2304*acc_g**3*k_beta**3))


      self.state['a1s'] = a1s
    def _BEM_calc_b1s(self):
        R = blade_params['radius']
        a = blade_params['a']
        c = blade_params['mean_chord']
        e = blade_params['e']
        th0 = blade_params['theta_0']
        th1 = blade_params['theta_1']
        k_beta = blade_params['k_beta']
        Mb = blade_params['Mb']
        Ib = blade_params['Ib']
        K = blade_params['K']
        Omega = self.omega
        v1 = self.vi
        mu = self.mu
        p = self.att[0]
        q = self.att[1]
  #WARNING!! DO NOT TOUCH LOOK AT OR ATTEMPT TO UNDERSTAND THIS:
        b1s =-rho_air*Omega*a*c*acc_g*((90*R**13*a**2*c**2*acc_g**2*mu**5*rho_air**2*th0 + 60*R**13*a**2*c**2*acc_g**2*mu**5*rho_air**2*th1 - 495*R**12*a**2*c**2*e*acc_g**2*mu**5*rho_air**2*th0 - 390*R**12*a**2*c**2*e*acc_g**2*mu**5*rho_air**2*th1 + 1080*R**11*a**2*c**2*e**2*acc_g**2*mu**5*rho_air**2*th0 + 1050*R**11*a**2*c**2*e**2*acc_g**2*mu**5*rho_air**2*th1 - 1125*R**10*a**2*c**2*e**3*acc_g**2*mu**5*rho_air**2*th0 - 1470*R**10*a**2*c**2*e**3*acc_g**2*mu**5*rho_air**2*th1 +
                                  450*R**9*a**2*c**2*e**4*acc_g**2*mu**5*rho_air**2*th0 + 1050*R**9*a**2*c**2*e**4*acc_g**2*mu**5*rho_air**2*th1 + 135*R**8*a**2*c**2*e**5*acc_g**2*mu**5*rho_air**2*th0 - 210*R**8*a**2*c**2*e**5*acc_g**2*mu**5*rho_air**2*th1 - 180*R**7*a**2*c**2*e**6*acc_g**2*mu**5*rho_air**2*th0 - 210*R**7*a**2*c**2*e**6*acc_g**2*mu**5*rho_air**2*th1 + 45*R**6*a**2*c**2*e**7*acc_g**2*mu**5*rho_air**2*th0 + 150*R**6*a**2*c**2*e**7*acc_g**2*mu**5*rho_air**2*th1 -
                                  30*R**5*a**2*c**2*e**8*acc_g**2*mu**5*rho_air**2*th1 + 120*R**13*a**2*alpha_s*c**2*acc_g**2*mu**4*rho_air**2 - 960*R**12*a**2*alpha_s*c**2*e*acc_g**2*mu**4*rho_air**2 + 3090*R**11*a**2*alpha_s*c**2*e**2*acc_g**2*mu**4*rho_air**2 - 5100*R**10*a**2*alpha_s*c**2*e**3*acc_g**2*mu**4*rho_air**2 + 4350*R**9*a**2*alpha_s*c**2*e**4*acc_g**2*mu**4*rho_air**2 - 1320*R**8*a**2*alpha_s*c**2*e**5*acc_g**2*mu**4*rho_air**2 - 690*R**7*a**2*alpha_s*c**2*e**6*acc_g**2*mu**4*rho_air**2 +
                                  660*R**6*a**2*alpha_s*c**2*e**7*acc_g**2*mu**4*rho_air**2 - 150*R**5*a**2*alpha_s*c**2*e**8*acc_g**2*mu**4*rho_air**2 - 90*R**13*a**2*c**2*acc_g**2*mu**3*rho_air**2*th0 - 48*R**13*a**2*c**2*acc_g**2*mu**3*rho_air**2*th1 + 195*R**12*a**2*c**2*e*acc_g**2*mu**3*rho_air**2*th0 + 68*R**12*a**2*c**2*e*acc_g**2*mu**3*rho_air**2*th1 + 420*R**11*a**2*c**2*e**2*acc_g**2*mu**3*rho_air**2*th0 + 678*R**11*a**2*c**2*e**2*acc_g**2*mu**3*rho_air**2*th1 -
                                  1650*R**10*a**2*c**2*e**3*acc_g**2*mu**3*rho_air**2*th0 - 2652*R**10*a**2*c**2*e**3*acc_g**2*mu**3*rho_air**2*th1 + 1500*R**9*a**2*c**2*e**4*acc_g**2*mu**3*rho_air**2*th0 + 3948*R**9*a**2*c**2*e**4*acc_g**2*mu**3*rho_air**2*th1 + 240*R**8*a**2*c**2*e**5*acc_g**2*mu**3*rho_air**2*th0 - 2436*R**8*a**2*c**2*e**5*acc_g**2*mu**3*rho_air**2*th1 - 1020*R**7*a**2*c**2*e**6*acc_g**2*mu**3*rho_air**2*th0 - 168*R**7*a**2*c**2*e**6*acc_g**2*mu**3*rho_air**2*th1 +
                                  330*R**6*a**2*c**2*e**7*acc_g**2*mu**3*rho_air**2*th0 + 972*R**6*a**2*c**2*e**7*acc_g**2*mu**3*rho_air**2*th1 + 150*R**5*a**2*c**2*e**8*acc_g**2*mu**3*rho_air**2*th0 - 348*R**5*a**2*c**2*e**8*acc_g**2*mu**3*rho_air**2*th1 - 75*R**4*a**2*c**2*e**9*acc_g**2*mu**3*rho_air**2*th0 - 48*R**4*a**2*c**2*e**9*acc_g**2*mu**3*rho_air**2*th1 + 34*R**3*a**2*c**2*e**10*acc_g**2*mu**3*rho_air**2*th1 - 240*R**13*a**2*alpha_s*c**2*acc_g**2*mu**2*rho_air**2 +
                                  1360*R**12*a**2*alpha_s*c**2*e*acc_g**2*mu**2*rho_air**2 - 2940*R**11*a**2*alpha_s*c**2*e**2*acc_g**2*mu**2*rho_air**2 + 2640*R**10*a**2*alpha_s*c**2*e**3*acc_g**2*mu**2*rho_air**2 - 1680*R**8*a**2*alpha_s*c**2*e**5*acc_g**2*mu**2*rho_air**2 + 840*R**7*a**2*alpha_s*c**2*e**6*acc_g**2*mu**2*rho_air**2 + 240*R**6*a**2*alpha_s*c**2*e**7*acc_g**2*mu**2*rho_air**2 - 240*R**5*a**2*alpha_s*c**2*e**8*acc_g**2*mu**2*rho_air**2 +
                                  20*R**3*a**2*alpha_s*c**2*e**10*acc_g**2*mu**2*rho_air**2 - 180*R**13*a**2*c**2*acc_g**2*mu*rho_air**2*th0 - 144*R**13*a**2*c**2*acc_g**2*mu*rho_air**2*th1 + 990*R**12*a**2*c**2*e*acc_g**2*mu*rho_air**2*th0 + 960*R**12*a**2*c**2*e*acc_g**2*mu*rho_air**2*th1 - 2080*R**11*a**2*c**2*e**2*acc_g**2*mu*rho_air**2*th0 - 2604*R**11*a**2*c**2*e**2*acc_g**2*mu*rho_air**2*th1 + 1890*R**10*a**2*c**2*e**3*acc_g**2*mu*rho_air**2*th0 + 3520*R**10*a**2*c**2*e**3*acc_g**2*mu*rho_air**2*th1 -
                                  360*R**9*a**2*c**2*e**4*acc_g**2*mu*rho_air**2*th0 - 2100*R**9*a**2*c**2*e**4*acc_g**2*mu*rho_air**2*th1 - 420*R**8*a**2*c**2*e**5*acc_g**2*mu*rho_air**2*th0 - 192*R**8*a**2*c**2*e**5*acc_g**2*mu*rho_air**2*th1 + 840*R**7*a**2*c**2*e**6*acc_g**2*mu*rho_air**2*th1 + 180*R**6*a**2*c**2*e**7*acc_g**2*mu*rho_air**2*th0 - 192*R**6*a**2*c**2*e**7*acc_g**2*mu*rho_air**2*th1 + 60*R**5*a**2*c**2*e**8*acc_g**2*mu*rho_air**2*th0 - 120*R**5*a**2*c**2*e**8*acc_g**2*mu*rho_air**2*th1 -
                                  90*R**4*a**2*c**2*e**9*acc_g**2*mu*rho_air**2*th0 + 36*R**3*a**2*c**2*e**10*acc_g**2*mu*rho_air**2*th1 + 10*R**2*a**2*c**2*e**11*acc_g**2*mu*rho_air**2*th0 - 4*R*a**2*c**2*e**12*acc_g**2*mu*rho_air**2*th1)*Omega**5 +
                                  (60*R**13*a**2*c**2*acc_g**2*mu**4*p*rho_air**2 - 300*R**12*a**2*c**2*e*acc_g**2*mu**4*p*rho_air**2 + 555*R**11*a**2*c**2*e**2*acc_g**2*mu**4*p*rho_air**2 - 390*R**10*a**2*c**2*e**3*acc_g**2*mu**4*p*rho_air**2 - 75*R**9*a**2*c**2*e**4*acc_g**2*mu**4*p*rho_air**2 + 240*R**8*a**2*c**2*e**5*acc_g**2*mu**4*p*rho_air**2 - 75*R**7*a**2*c**2*e**6*acc_g**2*mu**4*p*rho_air**2 - 30*R**6*a**2*c**2*e**7*acc_g**2*mu**4*p*rho_air**2 + 15*R**5*a**2*c**2*e**8*acc_g**2*mu**4*p*rho_air**2 -
                                   120*R**13*a**2*c**2*acc_g**2*mu**2*p*rho_air**2 + 500*R**12*a**2*c**2*e*acc_g**2*mu**2*p*rho_air**2 - 120*R**12*a**2*c**2*acc_g**2*mu**3*rho_air**2*v1 - 600*R**11*a**2*c**2*e**2*acc_g**2*mu**2*p*rho_air**2 + 960*R**11*a**2*c**2*e*acc_g**2*mu**3*rho_air**2*v1 - 240*R**10*a**2*c**2*e**3*acc_g**2*mu**2*p*rho_air**2 - 3090*R**10*a**2*c**2*e**2*acc_g**2*mu**3*rho_air**2*v1 + 1140*R**9*a**2*c**2*e**4*acc_g**2*mu**2*p*rho_air**2 + 5100*R**9*a**2*c**2*e**3*acc_g**2*mu**3*rho_air**2*v1 -
                                   960*R**8*a**2*c**2*e**5*acc_g**2*mu**2*p*rho_air**2 - 4350*R**8*a**2*c**2*e**4*acc_g**2*mu**3*rho_air**2*v1 + 300*R**7*a**2*c**2*e**6*acc_g**2*mu**2*p*rho_air**2 + 1320*R**7*a**2*c**2*e**5*acc_g**2*mu**3*rho_air**2*v1 + 690*R**6*a**2*c**2*e**6*acc_g**2*mu**3*rho_air**2*v1 - 60*R**5*a**2*c**2*e**8*acc_g**2*mu**2*p*rho_air**2 - 660*R**5*a**2*c**2*e**7*acc_g**2*mu**3*rho_air**2*v1 + 60*R**4*a**2*c**2*e**9*acc_g**2*mu**2*p*rho_air**2 + 150*R**4*a**2*c**2*e**8*acc_g**2*mu**3*rho_air**2*v1 -
                                   20*R**3*a**2*c**2*e**10*acc_g**2*mu**2*p*rho_air**2 + 240*R**12*a**2*c**2*acc_g**2*mu*rho_air**2*v1 - 1360*R**11*a**2*c**2*e*acc_g**2*mu*rho_air**2*v1 + 2940*R**10*a**2*c**2*e**2*acc_g**2*mu*rho_air**2*v1 - 2640*R**9*a**2*c**2*e**3*acc_g**2*mu*rho_air**2*v1 + 1680*R**7*a**2*c**2*e**5*acc_g**2*mu*rho_air**2*v1 - 840*R**6*a**2*c**2*e**6*acc_g**2*mu*rho_air**2*v1 - 240*R**5*a**2*c**2*e**7*acc_g**2*mu*rho_air**2*v1 + 240*R**4*a**2*c**2*e**8*acc_g**2*mu*rho_air**2*v1 -
                                   20*R**2*a**2*c**2*e**10*acc_g**2*mu*rho_air**2*v1 + 540*Ib*K*R**8*a*c*acc_g**2*mu**2*rho_air*v1 - 2520*Ib*K*R**7*a*c*e*acc_g**2*mu**2*rho_air*v1 + 4500*Ib*K*R**6*a*c*e**2*acc_g**2*mu**2*rho_air*v1 - 3600*Ib*K*R**5*a*c*e**3*acc_g**2*mu**2*rho_air*v1 + 900*Ib*K*R**4*a*c*e**4*acc_g**2*mu**2*rho_air*v1 + 360*Ib*K*R**3*a*c*e**5*acc_g**2*mu**2*rho_air*v1 - 180*Ib*K*R**2*a*c*e**6*acc_g**2*mu**2*rho_air*v1 - 540*Ib*R**9*a*c*acc_g**2*mu**2*q*rho_air +
                                   1800*Ib*R**8*a*c*e*acc_g**2*mu**2*q*rho_air - 1980*Ib*R**7*a*c*e**2*acc_g**2*mu**2*q*rho_air + 720*Ib*R**6*a*c*e**3*acc_g**2*mu**2*q*rho_air - 180*Ib*R**5*a*c*e**4*acc_g**2*mu**2*q*rho_air + 360*Ib*R**4*a*c*e**5*acc_g**2*mu**2*q*rho_air - 180*Ib*R**3*a*c*e**6*acc_g**2*mu**2*q*rho_air + 540*K*Mb*R**8*a*c*e*acc_g*mu**2*rho_air*v1 - 2520*K*Mb*R**7*a*c*e**2*acc_g*mu**2*rho_air*v1 + 4500*K*Mb*R**6*a*c*e**3*acc_g*mu**2*rho_air*v1 - 3600*K*Mb*R**5*a*c*e**4*acc_g*mu**2*rho_air*v1 +
                                   900*K*Mb*R**4*a*c*e**5*acc_g*mu**2*rho_air*v1 + 360*K*Mb*R**3*a*c*e**6*acc_g*mu**2*rho_air*v1 - 180*K*Mb*R**2*a*c*e**7*acc_g*mu**2*rho_air*v1 - 540*Mb*R**9*a*c*e*acc_g*mu**2*q*rho_air + 1800*Mb*R**8*a*c*e**2*acc_g*mu**2*q*rho_air - 1980*Mb*R**7*a*c*e**3*acc_g*mu**2*q*rho_air + 720*Mb*R**6*a*c*e**4*acc_g*mu**2*q*rho_air - 180*Mb*R**5*a*c*e**5*acc_g*mu**2*q*rho_air + 360*Mb*R**4*a*c*e**6*acc_g*mu**2*q*rho_air - 180*Mb*R**3*a*c*e**7*acc_g*mu**2*q*rho_air -
                                   1080*Ib*K*R**8*a*c*acc_g**2*rho_air*v1 + 5760*Ib*K*R**7*a*c*e*acc_g**2*rho_air*v1 - 12000*Ib*K*R**6*a*c*e**2*acc_g**2*rho_air*v1 + 11520*Ib*K*R**5*a*c*e**3*acc_g**2*rho_air*v1 - 3600*Ib*K*R**4*a*c*e**4*acc_g**2*rho_air*v1 - 1920*Ib*K*R**3*a*c*e**5*acc_g**2*rho_air*v1 + 1440*Ib*K*R**2*a*c*e**6*acc_g**2*rho_air*v1 - 120*Ib*K*a*c*e**8*acc_g**2*rho_air*v1 + 1080*Ib*R**9*a*c*acc_g**2*q*rho_air - 4320*Ib*R**8*a*c*e*acc_g**2*q*rho_air + 6000*Ib*R**7*a*c*e**2*acc_g**2*q*rho_air -
                                   2880*Ib*R**6*a*c*e**3*acc_g**2*q*rho_air - 480*Ib*R**4*a*c*e**5*acc_g**2*q*rho_air + 720*Ib*R**3*a*c*e**6*acc_g**2*q*rho_air - 120*Ib*R*a*c*e**8*acc_g**2*q*rho_air - 1080*K*Mb*R**8*a*c*e*acc_g*rho_air*v1 + 5760*K*Mb*R**7*a*c*e**2*acc_g*rho_air*v1 - 12000*K*Mb*R**6*a*c*e**3*acc_g*rho_air*v1 + 11520*K*Mb*R**5*a*c*e**4*acc_g*rho_air*v1 - 3600*K*Mb*R**4*a*c*e**5*acc_g*rho_air*v1 - 1920*K*Mb*R**3*a*c*e**6*acc_g*rho_air*v1 + 1440*K*Mb*R**2*a*c*e**7*acc_g*rho_air*v1 -
                                   120*K*Mb*a*c*e**9*acc_g*rho_air*v1 + 1080*Mb*R**9*a*c*e*acc_g*q*rho_air - 4320*Mb*R**8*a*c*e**2*acc_g*q*rho_air + 6000*Mb*R**7*a*c*e**3*acc_g*q*rho_air - 2880*Mb*R**6*a*c*e**4*acc_g*q*rho_air - 480*Mb*R**4*a*c*e**6*acc_g*q*rho_air + 720*Mb*R**3*a*c*e**7*acc_g*q*rho_air - 120*Mb*R*a*c*e**9*acc_g*q*rho_air)*Omega**4 +
                                   (-720*Mb*R**9*a*c*acc_g**2*mu**3*rho_air + 2520*Mb*R**8*a*c*e*acc_g**2*mu**3*rho_air - 2880*Mb*R**7*a*c*e**2*acc_g**2*mu**3*rho_air + 720*Mb*R**6*a*c*e**3*acc_g**2*mu**3*rho_air + 720*Mb*R**5*a*c*e**4*acc_g**2*mu**3*rho_air - 360*Mb*R**4*a*c*e**5*acc_g**2*mu**3*rho_air + 1440*Mb*R**9*a*c*acc_g**2*mu*rho_air - 6000*Mb*R**8*a*c*e*acc_g**2*mu*rho_air + 8640*Mb*R**7*a*c*e**2*acc_g**2*mu*rho_air - 3600*Mb*R**6*a*c*e**3*acc_g**2*mu*rho_air -
                                    2400*Mb*R**5*a*c*e**4*acc_g**2*mu*rho_air + 2160*Mb*R**4*a*c*e**5*acc_g**2*mu*rho_air - 240*Mb*R**2*a*c*e**7*acc_g**2*mu*rho_air - 17280*Ib*Mb*R**5*alpha_s*e*acc_g*mu**2 + 34560*Ib*Mb*R**4*alpha_s*e**2*acc_g*mu**2 - 17280*Ib*Mb*R**3*alpha_s*e**3*acc_g*mu**2 - 17280*Ib*R**5*alpha_s*acc_g**2*k_beta*mu**2 + 34560*Ib*R**4*alpha_s*e*acc_g**2*k_beta*mu**2 - 17280*Ib*R**3*alpha_s*e**2*acc_g**2*k_beta*mu**2 - 17280*Mb**2*R**5*alpha_s*e**2*mu**2 + 34560*Mb**2*R**4*alpha_s*e**3*mu**2 -
                                    17280*Mb**2*R**3*alpha_s*e**4*mu**2 - 17280*Mb*R**5*alpha_s*e*acc_g*k_beta*mu**2 + 34560*Mb*R**4*alpha_s*e**2*acc_g*k_beta*mu**2 - 17280*Mb*R**3*alpha_s*e**3*acc_g*k_beta*mu**2 - 23040*Ib*Mb*R**5*e*acc_g*mu*th0 - 17280*Ib*Mb*R**5*e*acc_g*mu*th1 + 34560*Ib*Mb*R**4*e**2*acc_g*mu*th0 + 46080*Ib*Mb*R**4*e**2*acc_g*mu*th1 - 34560*Ib*Mb*R**3*e**3*acc_g*mu*th1 - 11520*Ib*Mb*R**2*e**4*acc_g*mu*th0 + 5760*Ib*Mb*R*e**5*acc_g*mu*th1 - 23040*Ib*R**5*acc_g**2*k_beta*mu*th0 -
                                    17280*Ib*R**5*acc_g**2*k_beta*mu*th1 + 34560*Ib*R**4*e*acc_g**2*k_beta*mu*th0 + 46080*Ib*R**4*e*acc_g**2*k_beta*mu*th1 - 34560*Ib*R**3*e**2*acc_g**2*k_beta*mu*th1 - 11520*Ib*R**2*e**3*acc_g**2*k_beta*mu*th0 + 5760*Ib*R*e**4*acc_g**2*k_beta*mu*th1 - 23040*Mb**2*R**5*e**2*mu*th0 - 17280*Mb**2*R**5*e**2*mu*th1 + 34560*Mb**2*R**4*e**3*mu*th0 + 46080*Mb**2*R**4*e**3*mu*th1 - 34560*Mb**2*R**3*e**4*mu*th1 - 11520*Mb**2*R**2*e**5*mu*th0 + 5760*Mb**2*R*e**6*mu*th1 -
                                    23040*Mb*R**5*e*acc_g*k_beta*mu*th0 - 17280*Mb*R**5*e*acc_g*k_beta*mu*th1 + 34560*Mb*R**4*e**2*acc_g*k_beta*mu*th0 + 46080*Mb*R**4*e**2*acc_g*k_beta*mu*th1 - 34560*Mb*R**3*e**3*acc_g*k_beta*mu*th1 - 11520*Mb*R**2*e**4*acc_g*k_beta*mu*th0 + 5760*Mb*R*e**5*acc_g*k_beta*mu*th1)*Omega**3 +
                                    (540*K*R**8*a*c*acc_g**2*k_beta*mu**2*rho_air*v1 - 2520*K*R**7*a*c*e*acc_g**2*k_beta*mu**2*rho_air*v1 + 4500*K*R**6*a*c*e**2*acc_g**2*k_beta*mu**2*rho_air*v1 - 3600*K*R**5*a*c*e**3*acc_g**2*k_beta*mu**2*rho_air*v1 + 900*K*R**4*a*c*e**4*acc_g**2*k_beta*mu**2*rho_air*v1 + 360*K*R**3*a*c*e**5*acc_g**2*k_beta*mu**2*rho_air*v1 - 180*K*R**2*a*c*e**6*acc_g**2*k_beta*mu**2*rho_air*v1 - 540*R**9*a*c*acc_g**2*k_beta*mu**2*q*rho_air + 1800*R**8*a*c*e*acc_g**2*k_beta*mu**2*q*rho_air -
                                     1980*R**7*a*c*e**2*acc_g**2*k_beta*mu**2*q*rho_air + 720*R**6*a*c*e**3*acc_g**2*k_beta*mu**2*q*rho_air - 180*R**5*a*c*e**4*acc_g**2*k_beta*mu**2*q*rho_air + 360*R**4*a*c*e**5*acc_g**2*k_beta*mu**2*q*rho_air - 180*R**3*a*c*e**6*acc_g**2*k_beta*mu**2*q*rho_air - 1080*K*R**8*a*c*acc_g**2*k_beta*rho_air*v1 + 5760*K*R**7*a*c*e*acc_g**2*k_beta*rho_air*v1 - 12000*K*R**6*a*c*e**2*acc_g**2*k_beta*rho_air*v1 + 11520*K*R**5*a*c*e**3*acc_g**2*k_beta*rho_air*v1 - 3600*K*R**4*a*c*e**4*acc_g**2*k_beta*rho_air*v1 -
                                     1920*K*R**3*a*c*e**5*acc_g**2*k_beta*rho_air*v1 + 1440*K*R**2*a*c*e**6*acc_g**2*k_beta*rho_air*v1 - 120*K*a*c*e**8*acc_g**2*k_beta*rho_air*v1 + 1080*R**9*a*c*acc_g**2*k_beta*q*rho_air - 4320*R**8*a*c*e*acc_g**2*k_beta*q*rho_air + 6000*R**7*a*c*e**2*acc_g**2*k_beta*q*rho_air - 2880*R**6*a*c*e**3*acc_g**2*k_beta*q*rho_air - 480*R**4*a*c*e**5*acc_g**2*k_beta*q*rho_air + 720*R**3*a*c*e**6*acc_g**2*k_beta*q*rho_air - 120*R*a*c*e**8*acc_g**2*k_beta*q*rho_air - 8640*Ib*Mb*R**5*e*acc_g*p +
                                     11520*Ib*Mb*R**4*e**2*acc_g*p + 17280*Ib*Mb*R**4*e*acc_g*mu*v1 - 34560*Ib*Mb*R**3*e**2*acc_g*mu*v1 + 17280*Ib*Mb*R**2*e**3*acc_g*mu*v1 - 2880*Ib*Mb*R*e**5*acc_g*p - 8640*Ib*R**5*acc_g**2*k_beta*p + 11520*Ib*R**4*e*acc_g**2*k_beta*p + 17280*Ib*R**4*acc_g**2*k_beta*mu*v1 - 34560*Ib*R**3*e*acc_g**2*k_beta*mu*v1 + 17280*Ib*R**2*e**2*acc_g**2*k_beta*mu*v1 - 2880*Ib*R*e**4*acc_g**2*k_beta*p - 8640*Mb**2*R**5*e**2*p + 11520*Mb**2*R**4*e**3*p + 17280*Mb**2*R**4*e**2*mu*v1 -
                                     34560*Mb**2*R**3*e**3*mu*v1 + 17280*Mb**2*R**2*e**4*mu*v1 - 2880*Mb**2*R*e**6*p - 8640*Mb*R**5*e*acc_g*k_beta*p + 11520*Mb*R**4*e**2*acc_g*k_beta*p + 17280*Mb*R**4*e*acc_g*k_beta*mu*v1 - 34560*Mb*R**3*e**2*acc_g*k_beta*mu*v1 + 17280*Mb*R**2*e**3*acc_g*k_beta*mu*v1 - 2880*Mb*R*e**5*acc_g*k_beta*p)*Omega**2 +
                                    (-17280*Mb*R**5*alpha_s*e*acc_g*k_beta*mu**2 + 34560*Mb*R**4*alpha_s*e**2*acc_g*k_beta*mu**2 - 17280*Mb*R**3*alpha_s*e**3*acc_g*k_beta*mu**2 - 17280*R**5*alpha_s*acc_g**2*k_beta**2*mu**2 + 34560*R**4*alpha_s*e*acc_g**2*k_beta**2*mu**2 - 17280*R**3*alpha_s*e**2*acc_g**2*k_beta**2*mu**2 - 23040*Mb*R**5*e*acc_g*k_beta*mu*th0 - 17280*Mb*R**5*e*acc_g*k_beta*mu*th1 + 34560*Mb*R**4*e**2*acc_g*k_beta*mu*th0 + 46080*Mb*R**4*e**2*acc_g*k_beta*mu*th1 - 34560*Mb*R**3*e**3*acc_g*k_beta*mu*th1 -
                                     11520*Mb*R**2*e**4*acc_g*k_beta*mu*th0 + 5760*Mb*R*e**5*acc_g*k_beta*mu*th1 - 23040*R**5*acc_g**2*k_beta**2*mu*th0 - 17280*R**5*acc_g**2*k_beta**2*mu*th1 + 34560*R**4*e*acc_g**2*k_beta**2*mu*th0 + 46080*R**4*e*acc_g**2*k_beta**2*mu*th1 - 34560*R**3*e**2*acc_g**2*k_beta**2*mu*th1 - 11520*R**2*e**3*acc_g**2*k_beta**2*mu*th0 + 5760*R*e**4*acc_g**2*k_beta**2*mu*th1)*Omega - 8640*Mb*R**5*e*acc_g*k_beta*p + 11520*Mb*R**4*e**2*acc_g*k_beta*p + 17280*Mb*R**4*e*acc_g*k_beta*mu*v1 -
                                 34560*Mb*R**3*e**2*acc_g*k_beta*mu*v1 + 17280*Mb*R**2*e**3*acc_g*k_beta*mu*v1 - 2880*Mb*R*e**5*acc_g*k_beta*p - 8640*R**5*acc_g**2*k_beta**2*p + 11520*R**4*e*acc_g**2*k_beta**2*p + 17280*R**4*acc_g**2*k_beta**2*mu*v1 - 34560*R**3*e*acc_g**2*k_beta**2*mu*v1 + 17280*R**2*e**2*acc_g**2*k_beta**2*mu*v1 -
                                 2880*R*e**4*acc_g**2*k_beta**2*p)/(30*R*((-9*Ib*R**8*a**2*c**2*acc_g**3*mu**4*rho_air**2 + 36*Ib*R**7*a**2*c**2*e*acc_g**3*mu**4*rho_air**2 - 54*Ib*R**6*a**2*c**2*e**2*acc_g**3*mu**4*rho_air**2 + 36*Ib*R**5*a**2*c**2*e**3*acc_g**3*mu**4*rho_air**2 - 9*Ib*R**4*a**2*c**2*e**4*acc_g**3*mu**4*rho_air**2 - 9*Mb*R**8*a**2*c**2*e*acc_g**2*mu**4*rho_air**2 + 36*Mb*R**7*a**2*c**2*e**2*acc_g**2*mu**4*rho_air**2 - 54*Mb*R**6*a**2*c**2*e**3*acc_g**2*mu**4*rho_air**2 + 36*Mb*R**5*a**2*c**2*e**4*acc_g**2*mu**4*rho_air**2 -
                                          9*Mb*R**4*a**2*c**2*e**5*acc_g**2*mu**4*rho_air**2 + 36*Ib*R**8*a**2*c**2*acc_g**3*rho_air**2 - 192*Ib*R**7*a**2*c**2*e*acc_g**3*rho_air**2 + 400*Ib*R**6*a**2*c**2*e**2*acc_g**3*rho_air**2 - 384*Ib*R**5*a**2*c**2*e**3*acc_g**3*rho_air**2 + 120*Ib*R**4*a**2*c**2*e**4*acc_g**3*rho_air**2 + 64*Ib*R**3*a**2*c**2*e**5*acc_g**3*rho_air**2 - 48*Ib*R**2*a**2*c**2*e**6*acc_g**3*rho_air**2 + 4*Ib*a**2*c**2*e**8*acc_g**3*rho_air**2 + 36*Mb*R**8*a**2*c**2*e*acc_g**2*rho_air**2 -
                                          192*Mb*R**7*a**2*c**2*e**2*acc_g**2*rho_air**2 + 400*Mb*R**6*a**2*c**2*e**3*acc_g**2*rho_air**2 - 384*Mb*R**5*a**2*c**2*e**4*acc_g**2*rho_air**2 + 120*Mb*R**4*a**2*c**2*e**5*acc_g**2*rho_air**2 + 64*Mb*R**3*a**2*c**2*e**6*acc_g**2*rho_air**2 - 48*Mb*R**2*a**2*c**2*e**7*acc_g**2*rho_air**2 + 4*Mb*a**2*c**2*e**9*acc_g**2*rho_air**2)*Omega**6 +
                                           (-9*R**8*a**2*c**2*acc_g**3*k_beta*mu**4*rho_air**2 + 36*R**7*a**2*c**2*e*acc_g**3*k_beta*mu**4*rho_air**2 - 54*R**6*a**2*c**2*e**2*acc_g**3*k_beta*mu**4*rho_air**2 + 36*R**5*a**2*c**2*e**3*acc_g**3*k_beta*mu**4*rho_air**2 - 9*R**4*a**2*c**2*e**4*acc_g**3*k_beta*mu**4*rho_air**2 + 48*Mb*R**7*a**2*c**2*e**2*acc_g**2*mu**2*rho_air**2 - 168*Mb*R**6*a**2*c**2*e**3*acc_g**2*mu**2*rho_air**2 + 192*Mb*R**5*a**2*c**2*e**4*acc_g**2*mu**2*rho_air**2 -
                                            48*Mb*R**4*a**2*c**2*e**5*acc_g**2*mu**2*rho_air**2 - 48*Mb*R**3*a**2*c**2*e**6*acc_g**2*mu**2*rho_air**2 + 24*Mb*R**2*a**2*c**2*e**7*acc_g**2*mu**2*rho_air**2 + 48*R**7*a**2*c**2*e*acc_g**3*k_beta*mu**2*rho_air**2 - 168*R**6*a**2*c**2*e**2*acc_g**3*k_beta*mu**2*rho_air**2 + 192*R**5*a**2*c**2*e**3*acc_g**3*k_beta*mu**2*rho_air**2 - 48*R**4*a**2*c**2*e**4*acc_g**3*k_beta*mu**2*rho_air**2 - 48*R**3*a**2*c**2*e**5*acc_g**3*k_beta*mu**2*rho_air**2 +
                                            24*R**2*a**2*c**2*e**6*acc_g**3*k_beta*mu**2*rho_air**2 + 36*R**8*a**2*c**2*acc_g**3*k_beta*rho_air**2 - 192*R**7*a**2*c**2*e*acc_g**3*k_beta*rho_air**2 + 400*R**6*a**2*c**2*e**2*acc_g**3*k_beta*rho_air**2 - 384*R**5*a**2*c**2*e**3*acc_g**3*k_beta*rho_air**2 + 120*R**4*a**2*c**2*e**4*acc_g**3*k_beta*rho_air**2 + 64*R**3*a**2*c**2*e**5*acc_g**3*k_beta*rho_air**2 - 48*R**2*a**2*c**2*e**6*acc_g**3*k_beta*rho_air**2 + 4*a**2*c**2*e**8*acc_g**3*k_beta*rho_air**2)*Omega**4 +
                                             (2304*Ib*Mb**2*e**2*acc_g + 4608*Ib*Mb*e*acc_g**2*k_beta + 2304*Ib*acc_g**3*k_beta**2 + 2304*Mb**3*e**3 + 4608*Mb**2*e**2*acc_g*k_beta + 2304*Mb*e*acc_g**2*k_beta**2)*Omega**2 + 2304*Mb**2*e**2*acc_g*k_beta + 4608*Mb*e*acc_g**2*k_beta**2 + 2304*acc_g**3*k_beta**3))
          self.state['b1s'] =  b1s

      

    def _BEM_T_calc(self, _vi):
        #integrator params. fix later
        prop_state = {
            "a0":   a0,
            "a1s":  a1s,
            "b1s":  b1s,
            "Omega": Omega,
            "mu":   mu,
            "vi":   vi,
            "Vver": Vz,
        }
        T_integral = do_integral(diffThrust, prop_state)
        T_c = 3 * rho_air / (4*np.pi) * T_integral
        return T_c
    def _BEM_H_calc(self):
        #integrator params. fix later
        prop_state = {
            "a0":   a0,
            "a1s":  a1s,
            "b1s":  b1s,
            "Omega": Omega,
            "mu":   mu,
            "vi":   vi,
            "Vver": Vz,
        }
        H_integral = do_integral(diffHforce, prop_state)
        H_c = 3 * rho_air / (4*np.pi) * H_integral
        return H_c
    def _BEM_Q_calc(self):
        #integrator params. fix later
        prop_state = {
            "a0":   a0,
            "a1s":  a1s,
            "b1s":  b1s,
            "Omega": Omega,
            "mu":   mu,
            "vi":   vi,
            "Vver": Vz,
        }
        Q_integral = do_integral(diffTorque, prop_state)
        Q_c = 3 * rho_air / (4*np.pi) * Q_integral
        return Q_c
    
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
      A = blade_params['Area']
      K = blade_params['K']
      R = blade_params['radius']
      Vhor = self.mu*(self.omega*R)
      Vver = self.vel[2] - _vi
      T_BEM = BEM_T_calc(_vi)
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
      return v1
    def _update(self):
        self.vi = _BEM_calc_vi()
      
        self.state['a0'] = _BEM_calc_a0()
        self.state['a1s']= _BEM_calc_a1s()
        self.state['b1s']= _BEM_calc_b1s()
        self.state['T']  = _BEM_T_calc(self.vi)
        self.state['H']  = _BEM_H_calc()
        self.state['Q']  = _BEM_Q_calc()
  
        sgn2 = self.dir
        sgn1 = -self.dir
        self.Force = [
          -(self.state['H'] + np.sin(self.state['a1s'])*self.state['T']),
          sgn1*np.sin(self.state['b1s'])*self.state['T'],
          -self.state['T']*np.cos(self.state['a0']]
        self.Torque = [ 
          sgn1*blade_params['k_beta']*self.state['b1s'],
          blade_params['k_beta']*self.state['a1s'],
          sgn2*self.state['Q']]

    def SetVelocity(self, motOmega, Vel_CM, att_body):
        R = blade_params['radius']
        self.omega = motOmega
        self.vel = FLU_to_FRD(Vel_CM + np.cross(att_body,self.r))
        self.att = FLU_to_FRD(att_body)

        Vhor = np.sqrt(self.vel[0]**2 + self.vel[1]**2)
        self.mu = Vhor/(self.omega*R)
        self.alpha_s = np.atan2(self.vel[2], Vhor)
        _update()
 


