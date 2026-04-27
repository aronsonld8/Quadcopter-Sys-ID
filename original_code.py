from matplotlib.artist import get
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from scipy.spatial.transform import Rotation as R
from scipy import signal
from scipy.integrate import dblquad
from scipy.optimize import least_squares
from scipy.optimize import fsolve
from scipy.interpolate import griddata
import sys

calc_th = True
files = {"vertical_oscillation_vel1_seg_1.csv"}

#constants
'''
#propeller data: https://database.tytorobotics.com/propellers/zdq/hq-prop-5043
                 https://database.tytorobotics.com/tests/wek/dys-samguk-shu-hq-5043
'''
f_s = 400
etol = 1.0e-3
vi_init = 7.62
start_time = 0
stop_time = 5000
acc_g = 9.81
rho_air = 1.204
mass = 0.752
Ixx = 0.00262
Iyy = 0.00214
Izz = 0.0043
armLength = 0.13
dx = 0.078
dy = 0.1
dz = 0.027
TEff = 0.971087084744227
Body_area_x = 0.06 * 0.09      
Body_area_y = 0.1 * 0.09       
Body_area_z = 0.1 * 0.06   
Body_Cd_z = 1.17
Body_Cd_xy = 1.17    

integral_params = {
    "n_r" :   21,
    "n_th":   20
}

blade_params ={
  "radius":                       6.477e-2,    # m
  "pitch":                        21.77,       # deg
  "twist":                        -11.0,       # deg/m
  "chord_inner":                  1.7e-2,      # m
  "chord_outer":                  0.7e-2,      # m
  "mean_chord":                   1.3E-2,      #m
  "solidity":                     0.215,       # 1
  "distance_cog":                 2.7e-3,      # m
  "num_blades":                   3,          # 1
  "mass_blade":                   1.22e-3,     # kg
  "rel_hinge_offset":             0.1,         # 1
  "drag_coefficient":             13.54894,    # 1
  "lift_coefficient":             15.24214,    # 1
  "hinge_spring_constant":        5.89,        # Nm/rad
  "K":                            0,           # thrust distortion factor
  "ef":                           0.1,         # Put a number that gives senseful results, somewhere between 0.05 and 0.2
  "a":                            8.962,       # lift curve (for angles)
  "d":                            0.0833,      # Constant drag (for angles)
  "mb":                           1.22e-3,     # mass of one blade (measured)
  "cg":                           2.7E-3,      # distance between blade CoG and root of the blade
  "k_beta":                       7.571        # in [N*m/rad, based on coning
  }

def set_params():
  #for blade
  blade_params['Area']= np.pi*blade_params['radius']**2
  blade_params['e'] = blade_params['ef']*blade_params['radius']
  blade_params['Ib'] = blade_params['mb'] * blade_params['cg']**2   # blade moment of inertia around hinge (estimate)
  blade_params['Ib'] = blade_params['Ib'] + 1/12*blade_params['mb']*((blade_params['radius']/2)**2) # add moment of lenghty rod
  blade_params['Mb'] = blade_params['mb'] * acc_g * blade_params['cg']            # moment of blade around hinge (estimate)
  blade_params['theta_0'] = blade_params['pitch']*np.pi/180 # pitch angle theta0 of the propeller [deg] (datasheet)
  blade_params['theta_1'] = blade_params['twist']*np.pi/180 # blade twist (measured)

  #for integrator
  points_r, weights_r = np.polynomial.legendre.leggauss(integral_params['n_r'])
  #points_th, weights_th = np.polynomial.legendre.leggauss(integral_params['n_th'])
  integral_params['th_points'] = np.linspace(0,2*np.pi,integral_params['n_th'], endpoint=False)
  integral_params['th_weights'] = (2*np.pi)/(integral_params['n_th'])*np.ones(integral_params['n_th'])
  # Map r from [-1, 1] to [0, R]
  r_div2 = 0.5*blade_params['radius']
  integral_params['r_points'] = r_div2 * points_r + r_div2
  integral_params['r_weights'] = r_div2 * weights_r
  #print(f"th weights = {integral_params['th_weights']}")


def print_stats(d, name):
  if len(d) == 0:
    return
  print(f"{name} range: [{min(d):.2f},{max(d):.2f}], avg = {np.average(d):.2f}, sdev = {np.std(d):.2f}")


def print_stats3(vector3, name):
  dims = {'x','y','z'}
  for d in dims:
    print_stats(vector3[d], f"{name}_{d}")

def plot_path(path) :
  print_stats3(path,"r")
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  ax.plot(path['x'], path['y'], path['z'], label='3D Path')
  ax.set_xlabel('X Axis')
  ax.set_ylabel('Y Axis')
  ax.set_zlabel('Z Axis')
  plt.legend()
  plt.show()
def plot_xy(x,y,name_x,name_y):
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_xlabel(name_x)
  ax.set_ylabel(name_y)
  ax.set_title(f"{name_y} vs {name_x}")
  plt.grid(True)
  ax.plot(x, y)
  plt.show()

def plot_vector(vector, name):
  print_stats(vector, name)
  t = np.linspace(0, len(vector)/f_s, len(vector))
  plot_xy(t, vector, name, "Time (sec)")

def plot_scatter(x,y,name_x,name_y):
  plt.scatter(x, y, color='blue', alpha=0.5, s=25)
  plt.title(f"{name_y} vs {name_x}")
  plt.xlabel(name_x)
  plt.ylabel(name_y)
  plt.show()

def plot_vector3(vector3, name):
  print_stats3(vector3, name)

  t = np.linspace(0, len(vector3['x'])/f_s, len(vector3['x']))

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_xlabel('Time (sec)')
  ax.set_ylabel(name)
  ax.set_title(name)
  plt.grid(True)
  dims = {'x','y','z'}
  for d in dims:
    ax.plot(t, vector3[d], label=f"{name}_{d}")
    ax.legend()
    plt.legend()
    plt.tight_layout()
  plt.show()

def plot_bode(F_w, w, name, f_max = f_s/2):
  n = int(np.floor(len(w)/2))
  w = w[:n]
  F_w = F_w[:n]
  mag = 20 * np.log10(np.abs(F_w))
  phase = np.angle(F_w, deg=True)
  if(f_max != f_s/2):
    w_max = f_max*2*np.pi
    i_end = np.where(w > w_max)[0][0]
    w = w[:i_end]
    mag = mag[:i_end]
    phase = phase[:i_end]

  # Plot Magnitude
  plt.subplot(2, 1, 1)
  plt.semilogx(w, mag)    # Bode magnitude plot
  plt.title(name)
  plt.ylabel('Magnitude (dB)')
  # Plot Phase
  plt.subplot(2, 1, 2)
  plt.semilogx(w, phase)  # Bode phase plot
  plt.ylabel('Phase (deg)')
  plt.xlabel('Frequency (rad/sec)')
  plt.show()
def plot_spectrum(vector, name):
  F_vector = np.fft.fft(vector,len(vector))
  w = 2*np.pi*np.fft.fftfreq(len(F_vector),1/f_s)
  plot_bode(F_vector, w, name)
def plot_heatmap(x, y, F_xy, name_x, name_y, name_F):
  # Define the resolution of your heatmap
  res = 200
  # Create a regular grid
  xi = np.linspace(x.min(), x.max(), res)
  yi = np.linspace(y.min(), y.max(), res)
  X, Y = np.meshgrid(xi, yi)
  # Interpolate the 1D data onto the 2D grid
  Z = griddata((x, y), F_xy, (X, Y), method='linear')
  # Plot the heatmap
  plt.pcolormesh(X, Y, Z, cmap='viridis', shading='auto')
  plt.colorbar(label=name_F)
  plt.scatter(x, y, c='white', s=1, alpha=0.5) # Optional: overlay original data points
  plt.xlabel(name_x)
  plt.ylabel(name_y)
  plt.show()

def x_model(Vx, Theta, Ax):
  X = np.column_stack([Vx,Theta])
  model = LinearRegression().fit(X, Ax)
  ax_hat = model.predict(X)
  plot_vector(Ax, "Ax")
  plot_vector(Vx, "Vx")
  plot_vector(Ax-ax_hat, "residual")
  print(f"Linear Regression Coefficients: {model.coef_}")
  print(f"Linear Regression Intercept: {model.intercept_}")
def y_model(Vy, Phi, Ay):
  X = np.column_stack([Vy,Phi])
  model = LinearRegression().fit(X, Ay)
  ay_hat = model.predict(X)
  plot_vector(Ay, "Ay")
  plot_vector(Vy, "Vy")
  plot_vector(Ay-ay_hat, "residual")
  print(f"Linear Regression Coefficients: {model.coef_}")
  print(f"Linear Regression Intercept: {model.intercept_}")
def filter_data(raw_data, W_c):
  sos = signal.butter(4, W_c, 'low', fs=f_s, output='sos')
  filtered_signal = signal.sosfiltfilt(sos, raw_data)
  return filtered_signal
def verify_derivative(U, Udot):
  U_w = np.fft.fft(U)
  #do derivative
  w = 2*np.pi*np.fft.fftfreq(len(U_w), 1/f_s)
  Udot_w = 1j*w*U_w
  #inverse and get MSE
  Udot_hat = np.real(np.fft.ifft(U_w))
  MSE = np.mean((Udot-Udot_hat)**2)
  #plot_vector(Udot, "Udot")
  #plot_vector(Udot_hat, "Udot_hat")
  return MSE


def Az_model(params, u, v, w):
  vi_hov, Cd0, Cd1 = params
  Vh_2 = (u**2 + v**2)
  vi_dyn = (vi_hov**2)/np.sqrt(Vh_2 + (w - vi_hov**2)**2)
  T = 2*(rho_air/mass)*vi_dyn*np.sqrt(Vh_2 + (w - vi_dyn**2)**2)
  #T = 2*(rho_air)*blade_area*vi_dyn*np.sqrt(Vh_2 + (w - vi_dyn**2)**2)
  D = 0.5*(rho_air/mass)*(Cd0 + (Cd1*np.abs(w)))*w
  #D = Cd*w*np.abs(w)
  return (4*T + D)

def Az_cost_function(params, u, v, w, Az_meas):
  Az_hat = Az_model(params, u, v, w)
  return Az_meas - Az_hat


#integral stuff
def diffElement(r, Psi, prop_state, type='T'):
      #print(f"diffThrust(r={r},th={Psi*(180/np.pi)})")
      #get prop params
      a0 = prop_state['a0']
      a1s = prop_state['a1s']
      b1s = prop_state['b1s']
      Omega = prop_state['Omega']
      mu = prop_state['mu']
      Vz = prop_state['Vz']

      #get blade params
      R = blade_params['radius']
      K = blade_params['K']
      ci = blade_params['chord_inner']
      co = blade_params['chord_outer']
      cd = blade_params['drag_coefficient']
      cl = blade_params['lift_coefficient']
      theta_0 = blade_params['theta_0']
      theta_1 = blade_params['theta_1']

      r_ratio = r/R
      sPsi = np.sin(Psi)
      cPsi = np.cos(Psi)

      #induced velocity
      if prop_state['PP'] is True:
        vi = (R*Omega)*(prop_state['lamda_0'] + (r_ratio*prop_state['lamda_c']*cPsi) + (r_ratio*prop_state['lamda_s']*sPsi))
      else:
        vi = prop_state['vi']

      beta = a0 -  a1s*cPsi - b1s*sPsi
      U_T = Omega * (r + R*mu*sPsi)
      corr = K*r_ratio*cPsi
      U_P = Vz - vi*(1+corr) - r*Omega*(a1s*sPsi + b1s*cPsi) - Vz*beta*cPsi
      phi = np.atan2(U_P,U_T)
      #f = 1.5 * (1 - r_ratio) / (r_ratio*np.abs(np.sin(phi)) + 1.0e-8)
      #print(f"f = {f}, e^-f = {np.exp(-f)}")
      F = 1.0# 2/np.pi * np.arccos(np.exp(-f))
      U2 = U_T**2 + U_P**2
      #print(f"Up = {U_P}, Ut = {U_T}, phi = {phi/np.pi}pi, U2 = {U2}")
      alpha = theta_0 + (theta_1 * r_ratio) + phi;
      _cl = cl*np.sin(alpha)*np.cos(alpha);
      _cd = cd * np.sin(alpha)**2;
      _c = ci + (r_ratio * (co - ci))
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
def diffThrust(r,Psi, prop_state):
  return diffElement(r, Psi, prop_state, 'T')
def diffHforce(r,Psi, prop_state):
  return diffElement(r, Psi, prop_state, 'H')
def diffQtorque(r,Psi, prop_state):
  return diffElement(r, Psi, prop_state, 'Q')
def diffPitchMoment(r, Psi, prop_state):
  return diffThrust(r, Psi, prop_state) * r *np.cos(Psi)
def diffRollMoment(r, Psi, prop_state):
  return -diffThrust(r, Psi, prop_state) * r* np.sin(Psi)
def plot_quadrature(prop_state):
    x = np.empty(0)
    y = np.empty(0)
    F_xy =  np.empty(0)
    for i in range(integral_params['n_r']):
        for j in range(integral_params['n_th']):
          ri = integral_params['r_points'][i]
          thj = integral_params['th_points'][j]
          x = np.append(x, ri*np.cos(thj))
          y = np.append(y, ri*np.sin(thj))
          F_xy = np.append(F_xy, diffElement(ri, thj, prop_state,'T'))
    plot_heatmap(x, y, F_xy, "x(m)", "y(m)", "Thrust(N)")
    sys.exit(0)

def do_integral(func, args):
  ans_calc = 0
  for i in range(integral_params['n_r']):
    for j in range(integral_params['n_th']):
        r = integral_params['r_points'][i]
        th = integral_params['th_points'][j]
        ans_calc += integral_params['r_weights'][i] * integral_params['th_weights'][j] * func(r,th,args)

  return ans_calc


'''
BEM calculation functions.
reverse engineered from UZH matlab code
'''
def BEM_calc_a0(Omega, v1, mu, alpha_s, K, pdot, qdot):
  R = blade_params['radius']
  a = blade_params['a']
  c = blade_params['mean_chord']
  e = blade_params['e']
  th0 = blade_params['theta_0']
  th1 = blade_params['theta_1']
  k_beta = blade_params['k_beta']
  Mb = blade_params['Mb']
  Ib = blade_params['Ib']
  #WARNING!! DO NOT TOUCH LOOK OR ATTEMPT TO UNDERSTAND THIS:
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
  return a0

def BEM_calc_a1s(Omega, v1, mu, alpha_s, K, p, q):
  R = blade_params['radius']
  a = blade_params['a']
  c = blade_params['mean_chord']
  e = blade_params['e']
  th0 = blade_params['theta_0']
  th1 = blade_params['theta_1']
  k_beta = blade_params['k_beta']
  Mb = blade_params['Mb']
  Ib = blade_params['Ib']
  #WARNING!! DO NOT TOUCH LOOK OR ATTEMPT TO UNDERSTAND THIS:
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


  return a1s
def BEM_calc_b1s(Omega, v1, mu, alpha_s, K, p, q):
  R = blade_params['radius']
  a = blade_params['a']
  c = blade_params['mean_chord']
  e = blade_params['e']
  th0 = blade_params['theta_0']
  th1 = blade_params['theta_1']
  k_beta = blade_params['k_beta']
  Mb = blade_params['Mb']
  Ib = blade_params['Ib']
  #WARNING!! DO NOT TOUCH LOOK OR ATTEMPT TO UNDERSTAND THIS:
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
  return b1s


def BEM_T_calc(Vx, Vy, Vz, K, mu, Omega, vi, a0, a1s, b1s):
    #prop state
    prop_state = {
        "a0":   a0,
        "a1s":  a1s,
        "b1s":  b1s,
        "Omega": Omega,
        "mu":   mu,
        "vi":   vi,
        "Vz":   Vz,
        "PP" : False
    }
    #plot_quadrature(prop_state)

    #do integral
    #f_int_T = lambda Psi, r: _diffThrust(r,Psi)
    #T_integral, error = dblquad(f_int_T, 0, blade_params['radius'], 0, 2*np.pi)
    T_integral = do_integral(diffThrust, prop_state)
    #print(f"diff integrals = {np.abs(T_integral-T_integral2)}")
    T_c = 3 * rho_air / (4*np.pi) * T_integral
    return T_c
  
def BEM_H_calc(Vx, Vy, Vz, K, mu, Omega, vi, a0, a1s, b1s):
    #prop state
    prop_state = {
        "a0":   a0,
        "a1s":  a1s,
        "b1s":  b1s,
        "Omega": Omega,
        "mu":   mu,
        "vi":   vi,
        "Vz":   Vz,
        "PP" : False
    }
    #plot_quadrature(prop_state)

    #do integral
    H_integral = do_integral(diffHforce, prop_state)
    H_c = 3 * rho_air / (4*np.pi) * H_integral
    return H_c
  
def BEM_Q_calc(Vx, Vy, Vz, K, mu, Omega, vi, a0, a1s, b1s):
    #prop state
    prop_state = {
        "a0":   a0,
        "a1s":  a1s,
        "b1s":  b1s,
        "Omega": Omega,
        "mu":   mu,
        "vi":   vi,
        "Vz":   Vz,
        "PP" : False
    }
    #plot_quadrature(prop_state)

    #do integral
    Q_integral = do_integral(diffQtorque, prop_state)
    Q_c = 3 * rho_air / (4*np.pi) * Q_integral
    return Q_c

def BEM_calc_vh(T):
  v_hover = np.sqrt(T/(2.0 * blade_params['Area'] * rho_air))
  return v_hover
def BEM_gao_correction(Vz, Vh):
    k0=1;
    k1=-1.125
    k2=-1.372
    k3=-1.718
    k4=-0.655
    r = -Vz/Vh
    v2 = Vh*(k0 + k1*(r) + k2*(r)**2 + k3*(r)**3 + k4*(r)**4)
    #print(f"vortex ring, v1 = {v1:.2f}, Vz = {vz:.2f}, vh= {vh:.2f}, v2 = {v2:.2f} ")
    return v2
def BEM_fun(vi,Vx,Vy,Vz,Omega):
      A = blade_params['Area']
      K = blade_params['K']
      R = blade_params['radius']
      Vhor = np.sqrt(Vx**2 + Vy**2)
      Vver = Vz + vi
      mu = np.sqrt(Vhor/(Omega*R))
      T_BEM = BEM_T_calc(Vx, Vy, Vz ,K,mu,Omega,vi,0,0,0)
      T_M =  2 * rho_air * A * vi * np.sqrt(Vhor**2 + Vver**2)
      res = T_BEM - T_M
      return res
def BEM_calc_vi(Vx, Vy, Vz, Omega):
    vel = np.sqrt(np.dot(np.array([Vx, Vy, Vz]), np.array([Vx, Vy, Vz])))
    R = blade_params['radius'];
    mu = np.sqrt(Vx**2 + Vy**2)/(Omega*R)

    #print(f"start , vel = <{Vx:.6f},{Vy:.6f},{Vz:.6f}>, w= {Omega}")
    v1 = fsolve(BEM_fun,vi_init, args = (Vx,Vy,Vz,Omega),xtol=etol, maxfev=500)[0]
    v_ratio = Vz/v1
    if v_ratio >= 0 and v_ratio <= 2:
        vz = -Vz
        Vz = 0;
        vel = np.sqrt(np.dot(np.array([Vx, Vy]), np.array([Vx, Vy])))
        vh=fsolve(BEM_fun,vi_init,args = (Vx,Vy,Vz,Omega),xtol=etol, maxfev=500)[0]
        v2 = BEM_gao_correction(Vz,vh)
        #print(f"vortex ring, v1 = {v1:.2f}, Vz = {vz:.2f}, vh= {vh:.2f}, v2 = {v2:.2f} ")
        v1 = max(v1,v2)

    return v1

def get_coning(v_induced, Vx, Vy, Vz, Omega, rollrate, pitchrate, yawrate):
  n = len(Vx)
  R = blade_params['radius']
  K = blade_params['K']
  coning = np.zeros(n)

  for i in range(n):
    Vhor = np.sqrt(Vx[i]**2 + Vy[i]**2)
    mu = Vhor/(Omega[i]*R)
    alpha_s = np.atan2(-Vz[i], Vhor)
    coning[i] = BEM_calc_a0(Omega[i], v_induced[i], mu, alpha_s, K, rollrate[i], pitchrate[i])
  return coning
def get_flapping(v_induced, Vx, Vy, Vz, Omega, rollrate, pitchrate, yawrate):
  n = len(Vx)
  R = blade_params['radius']
  K = blade_params['K']
  flapping = {
      "lat"  : np.zeros(n),
      "long" : np.zeros(n)}
  for i in range(n):
    Vhor = np.sqrt(Vx[i]**2 + Vy[i]**2)
    mu = Vhor/(Omega[i]*R)
    alpha_s = np.atan2(-Vz[i], Vhor)
    flapping['lat'][i]  = BEM_calc_a1s(Omega[i], v_induced[i], mu, alpha_s, K, rollrate[i], pitchrate[i])
    flapping['long'][i] = BEM_calc_b1s(Omega[i], v_induced[i], mu, alpha_s, K, rollrate[i], pitchrate[i])
  return flapping

def BEM_get_prop_state(v_induced, Vel, Omega, att, dir):
  R = blade_params['radius']
  K = blade_params['K']

  Vhor = np.sqrt(Vel[0]**2 + Vel[1]**2)
  mu = Vhor/(Omega*R)
  alpha_s = np.atan2(Vel[2], Vhor)
  a0 = BEM_calc_a0(Omega, v_induced, mu, alpha_s, K, att[0], att[1])
  a1s = BEM_calc_a1s(Omega, v_induced, mu, alpha_s, K,att[0], att[1])
  b1s = BEM_calc_b1s(Omega, v_induced, mu, alpha_s, K,att[0], att[1])
  T = BEM_T_calc(Vel[0], Vel[1], Vel[2], K, mu, Omega, v_induced, a0, a1s, b1s)
  H = BEM_H_calc(Vel[0], Vel[1], Vel[2], K, mu, Omega, v_induced, a0, a1s, b1s)
  Q = BEM_Q_calc(Vel[0], Vel[1], Vel[2], K, mu, Omega, v_induced, a0, a1s, b1s)
  sgn2 = dir
  sgn1 = -dir

  res = [
      -(H + np.sin(a1s)*T),
      -(sgn1*np.sin(b1s)*T),
      T*np.cos(a0),
      sgn1*blade_params['k_beta']*b1s,
      -blade_params['k_beta']*a1s,
      -sgn2*Q
  ]
  return res

#pitt peters stuff
#https://www.researchgate.net/publication/238355026_Theoretical_prediction_of_dynamic-in_ow_derivatives
#https://ntrs.nasa.gov/api/citations/19900006622/downloads/19900006622.pdf

PP_Mass_Mtx = np.array([
      [8/(3*np.pi),0,0],
       [0,-16/(45*np.pi),0],
        [0,0,-16/(45*np.pi)]])

def PP_get_Vm(lamda0, Vel, Omega):
  R = blade_params['radius']
  RW  = R*Omega
  Vhor = np.sqrt(Vel[0]**2 + Vel[1]**2)
  mu = Vhor/RW
  Vax = Vel[2]/RW
  lamda = lamda0+Vax
  VT = np.sqrt(mu**2 + lamda**2)
  res = (mu**2 + lamda*(lamda + lamda0))/VT
  return res, VT
def PP_get_WS_angle(lamda0, Vel, Omega):
  R = blade_params['radius']
  RW  = R*Omega
  Vhor = np.sqrt(Vel[0]**2 + Vel[1]**2)
  mu = Vhor/RW
  Vax = Vel[2]/RW
  res = np.atan2(mu,Vax + lamda0)
  return res
def PP_calc_load(lamda, Vel, att, Omega):
  R = blade_params['radius']
  K = blade_params['K']
  A = blade_params['Area']
  RW  = R*Omega
  temp = 3 * rho_air / (4*np.pi)
  temp2 = rho_air*A*RW**2

  Vhor = np.sqrt(Vel[0]**2 + Vel[1]**2)
  mu = Vhor/RW
  alpha_s = np.atan2(Vel[2], Vhor)
  vi = lamda[0]*RW

  prop_state = {
        "a0":   BEM_calc_a0(Omega, vi, mu, alpha_s, K, att[0], att[1]),
        "a1s":  BEM_calc_a1s(Omega, vi, mu, alpha_s, K, att[0], att[1]),
        "b1s":  BEM_calc_b1s(Omega, vi, mu, alpha_s, K, att[0], att[1]),
        "Omega": Omega,
        "mu":   mu,
        "Vz":   Vel[2],
        "PP" : True,
        "lamda_0" : lamda[0],
        "lamda_s" : lamda[1],
        "lamda_c" : lamda[2]
    }
  Thrust = temp * do_integral(diffThrust, prop_state)
  CoefThrust  = Thrust/temp2
  CoefRoll  = temp * do_integral(diffRollMoment, prop_state)/(temp2*R)
  CoefPitch  = temp * do_integral(diffPitchMoment, prop_state)/(temp2*R)
  return CoefThrust, CoefRoll, CoefPitch, Thrust*np.cos(prop_state['a0'])
def PP_get_gain_mtx(lamda0, Vel, Omega):
  Vm, Vt = PP_get_Vm(lamda0,Vel,Omega)
  Chi = PP_get_WS_angle(lamda0, Vel, Omega)

  temp = (15*np.pi)/64
  tan = np.tan(0.5*Chi)
  cos = np.cos(Chi)

  LMtx = np.array([
      [0.5/Vt,0,(temp/Vm)*tan],
       [0,-(4/Vm)/(1+cos),0],
        [(temp/Vt)*tan,0,-(4/Vm)*cos/(1+cos)]])
  return LMtx

def PP_step(lamda, Vel, att, Omega):
  dt = 1/f_s
  RW = blade_params['radius']*Omega
  
  VRS = False
  lamda0_true = lamda[0]
  if Vel[2] >= 0.0:
    vh=fsolve(BEM_fun,vi_init,args = (Vel[0],Vel[1],0,Omega), xtol=etol, maxfev=500)[0]
    vratio = Vel[2]/vh
    if vratio <= 2.0:
      v2 = BEM_gao_correction(-Vel[2],vh)
      print(f"VRS zone, vh = {vh}, Vz = {Vel[2]}, Vcorr = {v2}")
      lamda0_true = v2/RW
      VRS = True

  CoefThrust, CoefRoll, CoefPitch, Thrust = PP_calc_load([lamda0_true,lamda[1],lamda[2]], Vel, att, Omega)
  Cvec = np.array([CoefThrust,CoefRoll,CoefPitch])

  #print(f"V0 = {lamda0_true*RW} Vm = {Vm*RW}, Chi = {Chi*180/np.pi} deg")

  LMtx = PP_get_gain_mtx(lamda0_true, Vel, Omega)
  TauMtx = LMtx@PP_Mass_Mtx
  #print(f"(1/v0){TauMtx*lamda[0]}")
  invTauMtx = np.linalg.inv(TauMtx)
  #print(f"eigs = {np.linalg.eigvals(-invTauMtx)}")
  lhs = np.identity(3) + dt*invTauMtx
  rhs = lamda + dt*(invTauMtx@LMtx)@Cvec
  lamda_new = np.linalg.inv(lhs)@rhs
  if VRS:
    lamda_new[0] = lamda0_true
  #print(f"lamda_new = {lamda_new}, lamda_old = {lamda}, Omega = {Omega}")
  #print("")
  return lamda_new,Thrust

def main():
  set_params()

  time_s = np.empty(0)
  pos = {
      "x" : np.empty(0),
      "y" : np.empty(0),
      "z" : np.empty(0)
  }
  pos_dot = {
      "x" : np.empty(0),
      "y" : np.empty(0),
      "z" : np.empty(0)
  }
  pos_dotdot = {
      "x" : np.empty(0),
      "y" : np.empty(0),
      "z" : np.empty(0)
  }
  omega = {
      "x" : np.empty(0),
      "y" : np.empty(0),
      "z" : np.empty(0)
  }
  '''
          ^+x
    mot4  |  mot2
   (CW)   |  (CCW)
   +y <---|------
          |
   mot3   |  mot1
   (CCW)  |  (CW)
          |

    '''
  mot= {
      "r1" : np.array([-dx, -dy, dz]),
      "r2" : np.array([ dx, -dy, dz]),
      "r3" : np.array([-dx,  dy, dz]),
      "r4" : np.array([ dx,  dy, dz]),
      "dir1" : -1,
      "dir2" : 1,
      "dir3" : 1,
      "dir4" : -1,
      "w1" : np.empty(0),
      "w2" : np.empty(0),
      "w3" : np.empty(0),
      "w4" : np.empty(0),
      "T1" : np.empty((0,6)),
      "T2" : np.empty((0,6)),
      "T3" : np.empty((0,6)),
      "T4" : np.empty((0,6))
  }
  euler_angles = {
      "roll" : np.empty(0),
      "pitch" : np.empty(0),
      "yaw" : np.empty(0)
  }
  grav = {
      "x" : np.empty(0),
      "y" : np.empty(0),
      "z" : np.empty(0)
  }

  #open files
  for f in files:

    df = pd.read_csv(f)
    df = df[(df['t'] >= start_time) & (df['t'] <= stop_time)]
    if len(df['pos x']) == 0:
       continue

    #time
    time_s = np.append(time_s,df['t'])
    print(f"{len(time_s)} samples")

    #position
    pos['x'] = np.append(pos['x'],df['pos x'])
    pos['y'] = np.append(pos['y'],df['pos y'])
    pos['z'] = np.append(pos['z'],df['pos z'])
    plot_path(pos)

    #omega
    omega['x'] = np.append(omega['x'],df['ang vel x'])
    omega['y'] = np.append(omega['y'],df['ang vel y'])
    omega['z'] = np.append(omega['z'],df['ang vel z'])

    #acceleration
    pos_dotdot['x'] = np.append(pos_dotdot['x'],df['acc x'])
    pos_dotdot['y'] = np.append(pos_dotdot['y'],df['acc y'])
    pos_dotdot['z'] = np.append(pos_dotdot['z'],df['acc z'])

    #motor
    mot['w1'] = np.append(mot['w1'],df['mot 1'])
    mot['w2'] = np.append(mot['w2'],df['mot 2'])
    mot['w3'] = np.append(mot['w3'],df['mot 3'])
    mot['w4'] = np.append(mot['w4'],df['mot 4'])


    #quaternions
    quaternions_xyzw = df[['quat x','quat y','quat z','quat w']].to_numpy()
    rot = R.from_quat(quaternions_xyzw)
    euler_angles_radians = rot.as_euler('xyz', degrees=False)
    euler_angles['roll'] = np.append(euler_angles['roll'],euler_angles_radians[:,0])
    euler_angles['pitch'] = np.append(euler_angles['pitch'],euler_angles_radians[:,1])
    euler_angles['yaw'] = np.append(euler_angles['yaw'],euler_angles_radians[:,2])

    #velocity
    pos_dot['x'] = np.append(pos_dot['x'], df['vel x'])
    pos_dot['y'] = np.append(pos_dot['y'], df['vel y'])
    pos_dot['z'] = np.append(pos_dot['z'], df['vel z'])
    plot_xy(time_s, pos_dot['z'],"time", "Vz")

    #grav
    grav_b = rot.inv().apply(np.array([0,0,acc_g]))
    grav['x'] = np.append(grav['x'], grav_b[:,0])
    grav['y'] = np.append(grav['y'], grav_b[:,1])
    grav['z'] = np.append(grav['z'], grav_b[:,2])

    #Aero Acc
    omega_cross_v = {
      "x" : omega['y']*pos_dot['z'] - omega['z']*pos_dot['y'],
      "y" : omega['z']*pos_dot['x'] - omega['x']*pos_dot['z'],
      "z" : omega['x']*pos_dot['y'] - omega['y']*pos_dot['x']
    }
    Ax_aero = pos_dotdot['x'] - grav['x'] - omega_cross_v['x']
    Ay_aero = pos_dotdot['y'] - grav['y'] - omega_cross_v['y']
    Az_aero = pos_dotdot['z'] - grav['z'] - omega_cross_v['z']

    #thrust
    if calc_th is True:
      vrel_1 = {
          "x" : pos_dot['x'] + omega['y']*mot['r1'][2] - omega['z']*mot['r1'][1],
          "y" : pos_dot['y'] + omega['z']*mot['r1'][0] - omega['x']*mot['r1'][2],
          "z" : pos_dot['z'] + omega['x']*mot['r1'][1] - omega['y']*mot['r1'][0]
      }
      vrel_2 = {
          "x" : pos_dot['x'] + omega['y']*mot['r2'][2] - omega['z']*mot['r2'][1],
          "y" : pos_dot['y'] + omega['z']*mot['r2'][0] - omega['x']*mot['r2'][2],
          "z" : pos_dot['z'] + omega['x']*mot['r2'][1] - omega['y']*mot['r2'][0]
      }
      vrel_3 = {
          "x" : pos_dot['x'] + omega['y']*mot['r3'][2] - omega['z']*mot['r3'][1],
          "y" : pos_dot['y'] + omega['z']*mot['r3'][0] - omega['x']*mot['r3'][2],
          "z" : pos_dot['z'] + omega['x']*mot['r3'][1] - omega['y']*mot['r3'][0]
      }
      vrel_4 = {
          "x" : pos_dot['x'] + omega['y']*mot['r4'][2] - omega['z']*mot['r4'][1],
          "y" : pos_dot['y'] + omega['z']*mot['r4'][0] - omega['x']*mot['r4'][2],
          "z" : pos_dot['z'] + omega['x']*mot['r4'][1] - omega['y']*mot['r4'][0]
      }


      T1 = np.empty((len(time_s),6))
      T2 = np.empty((len(time_s),6))
      T3 = np.empty((len(time_s),6))
      T4 = np.empty((len(time_s),6))
      for i in range(len(time_s)):
        att = np.array([omega['x'][i], omega['y'][i], omega['z'][i]])
        Vel1 = np.array([vrel_1['x'][i], -vrel_1['y'][i], -vrel_1['z'][i]])
        Vel2 = np.array([vrel_2['x'][i], -vrel_2['y'][i], -vrel_2['z'][i]])
        Vel3 = np.array([vrel_3['x'][i], -vrel_3['y'][i], -vrel_3['z'][i]])
        Vel4 = np.array([vrel_4['x'][i], -vrel_4['y'][i], -vrel_4['z'][i]])
        
        T1[i] = BEM_get_prop_state(BEM_calc_vi(Vel1[0], Vel1[1], Vel1[2], mot['w1'][i]), Vel1, mot['w1'][i], att, mot['dir1'])
        T2[i] = BEM_get_prop_state(BEM_calc_vi(Vel2[0], Vel2[1], Vel2[2], mot['w2'][i]), Vel2, mot['w2'][i], att, mot['dir2'])
        T3[i] = BEM_get_prop_state(BEM_calc_vi(Vel3[0], Vel3[1], Vel3[2], mot['w3'][i]), Vel3, mot['w3'][i], att, mot['dir3'])
        T4[i] = BEM_get_prop_state(BEM_calc_vi(Vel4[0], Vel4[1], Vel4[2], mot['w4'][i]), Vel4, mot['w4'][i], att, mot['dir4'])
        
        print(f"{time_s[i]:.4f}sec DONE")


      #append df
      df['Fp_x_1'] = T1[:,0]
      df['Fp_x_2'] = T2[:,0]
      df['Fp_x_3'] = T3[:,0]
      df['Fp_x_4'] = T4[:,0]
      df['Fp_y_1'] = T1[:,1]
      df['Fp_y_2'] = T2[:,1]
      df['Fp_y_3'] = T3[:,1]
      df['Fp_y_4'] = T4[:,1]
      df['Fp_z_1'] = T1[:,2]
      df['Fp_z_2'] = T2[:,2]
      df['Fp_z_3'] = T3[:,2]
      df['Fp_z_4'] = T4[:,2]
      df['Tp_x_1'] = T1[:,3]
      df['Tp_x_2'] = T2[:,3]
      df['Tp_x_3'] = T3[:,3]
      df['Tp_x_4'] = T4[:,3]
      df['Tp_y_1'] = T1[:,4]
      df['Tp_y_2'] = T2[:,4]
      df['Tp_y_3'] = T3[:,4]
      df['Tp_y_4'] = T4[:,4]
      df['Tp_z_1'] = T1[:,5]
      df['Tp_z_2'] = T2[:,5]
      df['Tp_z_3'] = T3[:,5]
      df['Tp_z_4'] = T4[:,5]
      #save file
      df.to_csv(f"proc_{f}", index=False)

    mot['T1'] = np.append(mot['T1'], np.array(df['Fp_x_1'], df['Fp_y_1'], df['Fp_z_1'], df['Tp_x_1'], df['Tp_y_1'], df['Tp_z_1']))
    mot['T2'] = np.append(mot['T2'], np.array(df['Fp_x_2'], df['Fp_y_2'], df['Fp_z_2'], df['Tp_x_2'], df['Tp_y_2'], df['Tp_z_2']))
    mot['T3'] = np.append(mot['T3'], np.array(df['Fp_x_3'], df['Fp_y_3'], df['Fp_z_3'], df['Tp_x_3'], df['Tp_y_3'], df['Tp_z_3']))
    mot['T4'] = np.append(mot['T4'], np.array(df['Fp_x_4'], df['Fp_y_4'], df['Fp_z_4'], df['Tp_x_4'], df['Tp_y_4'], df['Tp_z_4']))


  #Fz = T - W - D
  plot_xy(time_s, mass*Az_aero,"time", "Fz")
  Thrust_sim = TEff*(mot['T1'][:,2] + mot['T2'][:,2] + mot['T3'][:,2] + mot['T4'][:,2])
  Dragz =  0.5*rho_air*Body_area_z*Body_Cd_z*pos_dot['z']*np.abs(pos_dot['z'])
  Fz_hat = Thrust_sim - mass*grav['z'] - Dragz

  plot_xy(time_s, Thrust_sim, "time", "Thrust_sim")
  plot_xy(time_s, (mass*Az_aero) - Fz_hat,"time", "residual")

  '''
  plot_xy(time_s, Thrust - mass*grav['z'],"time", "Thrust - W")
  plot_xy(time_s, Dragz,"time", "Dragz")
  vrel_1 = {
          "x" : pos_dot['x'] + omega['y']*mot['r1'][2] - omega['z']*mot['r1'][1],
          "y" : pos_dot['y'] + omega['z']*mot['r1'][0] - omega['x']*mot['r1'][2],
          "z" : pos_dot['z'] + omega['x']*mot['r1'][1] - omega['y']*mot['r1'][0]
  }
  vrel_2 = {
          "x" : pos_dot['x'] + omega['y']*mot['r2'][2] - omega['z']*mot['r2'][1],
          "y" : pos_dot['y'] + omega['z']*mot['r2'][0] - omega['x']*mot['r2'][2],
          "z" : pos_dot['z'] + omega['x']*mot['r2'][1] - omega['y']*mot['r2'][0]
  }
  vrel_3 = {
          "x" : pos_dot['x'] + omega['y']*mot['r3'][2] - omega['z']*mot['r3'][1],
          "y" : pos_dot['y'] + omega['z']*mot['r3'][0] - omega['x']*mot['r3'][2],
          "z" : pos_dot['z'] + omega['x']*mot['r3'][1] - omega['y']*mot['r3'][0]
  }
  vrel_4 = {
          "x" : pos_dot['x'] + omega['y']*mot['r4'][2] - omega['z']*mot['r4'][1],
          "y" : pos_dot['y'] + omega['z']*mot['r4'][0] - omega['x']*mot['r4'][2],
          "z" : pos_dot['z'] + omega['x']*mot['r4'][1] - omega['y']*mot['r4'][0]
  }
  mdf = pd.DataFrame()
  CT1 = TEff*mot['T1']/(rho_air*blade_params['Area']*(blade_params['radius']*mot['w1'])**2)
  CT2 = TEff*mot['T2']/(rho_air*blade_params['Area']*(blade_params['radius']*mot['w2'])**2)
  CT3 = TEff*mot['T3']/(rho_air*blade_params['Area']*(blade_params['radius']*mot['w3'])**2)
  CT4 = TEff*mot['T4']/(rho_air*blade_params['Area']*(blade_params['radius']*mot['w4'])**2)
  mdf['CT'] = np.concatenate((CT1,CT2,CT3,CT4))
  mu1 = np.sqrt(vrel_1['x']**2 + vrel_1['y']**2)/(blade_params['radius']*mot['w1'])
  mu2 = np.sqrt(vrel_2['x']**2 + vrel_2['y']**2)/(blade_params['radius']*mot['w2'])
  mu3 = np.sqrt(vrel_3['x']**2 + vrel_3['y']**2)/(blade_params['radius']*mot['w3'])
  mu4 = np.sqrt(vrel_4['x']**2 + vrel_4['y']**2)/(blade_params['radius']*mot['w4'])
  mdf['mu'] = np.concatenate((mu1,mu2,mu3,mu4))
  alpha1 = np.arctan2(-vrel_1['z'],np.sqrt(vrel_1['x']**2 + vrel_1['y']**2))*180/np.pi
  alpha2 = np.arctan2(-vrel_2['z'],np.sqrt(vrel_2['x']**2 + vrel_2['y']**2))*180/np.pi
  alpha3 = np.arctan2(-vrel_3['z'],np.sqrt(vrel_3['x']**2 + vrel_3['y']**2))*180/np.pi
  alpha4 = np.arctan2(-vrel_4['z'],np.sqrt(vrel_4['x']**2 + vrel_4['y']**2))*180/np.pi
  mdf['alpha'] = np.concatenate((alpha1,alpha2,alpha3,alpha4))
  vz1 = vrel_1['z']/(blade_params['radius']*mot['w1'])
  vz2 = vrel_2['z']/(blade_params['radius']*mot['w2'])
  vz3 = vrel_3['z']/(blade_params['radius']*mot['w3'])
  vz4 = vrel_4['z']/(blade_params['radius']*mot['w4'])
  mdf['vz'] = np.concatenate((vz1,vz2,vz3,vz4))
  mdf['w'] = np.concatenate((mot['w1'],mot['w2'],mot['w3'],mot['w4']))

  plot_heatmap(mdf['alpha'], mdf['mu'], mdf['CT'], "a_tpp", "mu", "CT")
  plot_scatter(mdf['vz'], mdf['CT'], "vz", "CT")
  '''









'''
W_c = 10
Vz_f = filter_data(pos_dot['z'],W_c)
Vy_f = filter_data(pos_dot['y'],W_c)
Vx_f = filter_data(pos_dot['x'],W_c)
Thrust_f = filter_data(Thrust,W_c)
Az_aero_f = filter_data(Az_aero,W_c)

initial_guess = [vi_hov0, 0.15, -0.015]
res = least_squares(
    Az_cost_function,
    initial_guess,
    args=(Vx_f, Vy_f, Vz_f, Az_aero_f),
    loss='linear',
    method='lm'
)
vi_hov, Cd0, Cd1 = res.x
print(f"vi_hov: {vi_hov}, Cd: {Cd0} + {Cd1}Vz")

Kw = 0.5*(Cd0 + Cd1*np.abs(Vz_f))*(rho_air/mass)
plot_xy(Vz_f, Kw,"Vz", "Kw")
'''

if __name__ == "__main__":
  main()


