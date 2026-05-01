from matplotlib.artist import get
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
from scipy import signal
from scipy.optimize import fsolve
from scipy.interpolate import griddata
import sys

#main constants
f_s = 400
etol = 1.0e-3
vi_init = 7.62
acc_g = 9.81
rho_air = 1.204
dx = 0.078  #m
dy = 0.1    #m  
dz = 0.027  #m

#plotting funcs
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

def plot_xy3(t,x,y,z,name):
  plt.figure(figsize=(10, 5))
  plt.plot(t, x, label=f"{name}_x", color='blue')
  plt.plot(t, y, label=f"{name}_y", color='red')
  plt.plot(t, z, label=f"{name}_z", color='green')
  plt.title(f"{name} vs time(sec)")
  plt.grid(True)
  plt.legend()
  plt.show()

def plot_scatter(x,y,name_x,name_y):
  plt.scatter(x, y, color='blue', alpha=0.5, s=25)
  plt.title(f"{name_y} vs {name_x}")
  plt.grid(True)
  plt.xlabel(name_x)
  plt.ylabel(name_y)
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
def plot_PSD(f_t, name):
  Pxx,freq = plt.psd(f_t, NFFT=256, Fs=f_s)
  plt.title(f"Power Spectral Density {name}")
  plt.show()
  cutoff_index = np.where(Pxx < 0.5 * np.max(Pxx))[0][0]
  fc = freq[cutoff_index]
  print(f"f_c = {fc}")

#signal processing stuff
def filter_data(f_t):
  freq,Pxx = signal.welch(f_t,
                          fs=f_s,
                          window='hann',   # plt.psd default window
                          nperseg=256,     # Match NFFT
                          noverlap=0,      # plt.psd default is 0; welch default is nperseg//2
                          detrend=False,   # plt.psd default is no detrend
                          scaling='density')
  max_idx = np.argmax(Pxx)
  temp = np.where(Pxx < 0.5 * np.max(Pxx))[0]
  cutoff_index=0
  for idx in temp:
    if idx > max_idx:
      cutoff_index = idx
      break
  f_c = freq[cutoff_index]
  arr = np.column_stack((freq, Pxx))
  print(f"f_c = {f_c}")
  sos = signal.butter(4, f_c, 'low', fs=f_s, output='sos')
  filtered_signal = signal.sosfiltfilt(sos, f_t)
  return filtered_signal
def get_derivative(U):
  U_w = np.fft.fft(U)
  #do derivative
  w = 2*np.pi*np.fft.fftfreq(len(U_w), 1/f_s)
  Udot_w = 1j*w*U_w
  Udot = np.real(np.fft.ifft(U_w))
  return Udot
