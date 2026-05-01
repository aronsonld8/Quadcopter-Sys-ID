import common

#blade
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


#motor
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
mot_params= {
      "r1" : np.array([-dx, -dy, dz]),
      "r2" : np.array([ dx, -dy, dz]),
      "r3" : np.array([-dx,  dy, dz]),
      "r4" : np.array([ dx,  dy, dz]),
      "dir1" : -1,
      "dir2" : 1,
      "dir3" : 1,
      "dir4" : -1
}

#body
body_params= {
  "mass" : 0.752,  #kg
  "Ixx"  : 0.00262,  #kgm2
  "Iyy"  : 0.00214,  #kgm2
  "Izz"  : 0.0043,  #kgm2
  "Ax" : 0.06 * 0.09,  #m2
  "Ay" :  0.1 * 0.09, #m2
  "Az" : 0.1 * 0.06, #m2
  "Cx" : 1.17,
  "Cy" : 1.17,
  "Cz" : 1.17
}
