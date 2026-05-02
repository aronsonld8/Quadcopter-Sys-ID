import common

class state:
  def __init__(self):
    #declarations
    self.pos_x = np.empty()
    self.pos_y = np.empty()
    self.pos_z = np.empty()
    self.u = np.empty()
    self.v = np.empty()
    self.w = np.empty()
    self.roll = np.empty()
    self.pitch = np.empty()
    self.yaw = np.empty()
    self.p = np.empty()
    self.q = np.empty()
    self.r = np.empty()
    self.mot_omega  = np.empty((0,4))
    self.mot_domega = np.empty((0,4))
    self.Force = np.empty((0,3))
    self.Torque = np.empty((0,3))
    self.grav = np.empty((0,3))
    self.F_res = np.empty((0,3))
    self.Tau_res = np.empty((0,3))

  def parse_datalog(self, file):
    df = pd.read_csv(file)
    #position
    self.pos_x = np.append(self.pos_x,df['pos x'])
    self.pos_y = np.append(self.pos_y,df['pos y'])
    self.pos_z = np.append(self.pos_z,df['pos z'])

    #velocity
    self.u = np.append(self.u, df['vel x'])
    self.v = np.append(self.v, df['vel y'])
    self.w = np.append(self.w, df['vel z'])
    
    #omega
    self.p = np.append(self.p,df['ang vel x'])
    self.q = np.append(self.q,df['ang vel y'])
    self.r = np.append(self.r,df['ang vel z'])

    #Force
    self.Force = np.append(self.Force, df[['acc x','acc y','acc z'].to_numpy(),axis =0)
    self.Force *= mass
    self.F_res = np.zeros((len(self.Force),3))
    
    #Torque
    self.Torque = np.append(self.Force, df[['ang acc x','ang acc y','ang acc z'].to_numpy(),axis =0)
    self.Torque[:,0] *= Ixx
    self.Torque[:,1] *= Iyy
    self.Torque[:,2] *= Izz
    self.Tau_res = np.zeros((len(self.Force),3))

    #motor
    self.mot_omega = np.append(self.Force, df[['mot 1','mot 2','mot 3','mot 4'].to_numpy(),axis =0)
    self.mot_domega = np.append(self.Force, df[['dmot 1','dmot 2','dmot 3','dmot 4'].to_numpy(),axis =0)

    #angles
    quaternions_xyzw = df[['quat x','quat y','quat z','quat w']].to_numpy()
    rot = R.from_quat(quaternions_xyzw)
    euler_angles_radians = rot.as_euler('xyz', degrees=False)
    self.roll = np.append(self.roll ,euler_angles_radians[:,0])
    self.pitch = np.append(self.pitch,euler_angles_radians[:,1])
    self.yaw = np.append(self.yaw,euler_angles_radians[:,2])
    
    #grav
    grav_b = rot.inv().apply(np.array([0,0,-acc_g]))
    self.grav = np.append(self.grav, grav_b)

  def write_output(self, file):
    out_df = pd.DataFrame({
        "x" : self.pos_x,
        "y" : self.pos_y,
        "z" : self.pos_y,
        "roll" : self.roll,
        "pitch" : self.pitch,
        "yaw" : self.yaw,
        "u"   : self.u,
        "v"   : self.v,
        "w"   : self.w,
        "p"   : self.p,
        "q"   : self.q,
        "r"   : self.r,
        "mot1" : self.mot_omega[:,0],
        "mot2" : self.mot_omega[:,1],
        "mot3" : self.mot_omega[:,2],
        "mot4" : self.mot_omega[:,3],
        "dmot1" : self.mot_domega[:,0],
        "dmot2" : self.mot_domega[:,1],
        "dmot3" : self.mot_domega[:,2],
        "dmot4" : self.mot_domega[:,3],
        "Fx" : self.Force[:,0],
        "Fy" : self.Force[:,1],
        "Fz" : self.Force[:,2],
        "F_resx" : self.F_res[:,0],
        "F_resy" : self.F_res[:,1],
        "F_resz" : self.F_res[:,2],
        "Taux" : self.Torque[:,0],
        "Tauy" : self.Torque[:,1],
        "Tauz" : self.Torque[:,2],
        "Tau_resx" : self.Tau_res[:,0],
        "Tau_resy" : self.Tau_res[:,1],
        "Tau_resz" : self.Tau_res[:,2],
        "grav_x" : self.grav[:,0],
        "grav_y" : self.grav[:,1],
        "grav_z" : self.grav[:,2]})
      out_df.to_csv(outfile, index=False)
    
