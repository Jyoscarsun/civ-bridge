# Import necessary modules to be implemented in the program
import numpy as np
import matplotlib.pyplot as plt
import os
import copy
import math


# MUST CHANGE THE FOLLOWING VARIABLES
# Overall bridge dimension
L = 1270 # length of the bridge
h = 76.27  # height of the bridge


# DO NOT CHANGE THE FOLLOWING VARIABLES (matboard and glue parameters)
t = 1.27        # thickness of matboard
sigma_t = 30    # matboard tensile strength
sigma_c = 6     # matboard compressive strength
tau_m = 4       # matboard shear strength
E_m = 4000      # matboard Young's modulus
nu_m = 0.2      # matboard Poisson's ratio 
tau_g = 2       # Ideal value of glue shear strength given cured property


def A_tot(shapes_and_dim):
    # Calculate the total cross-sectional area
    res = 0
    for number in shapes_and_dim.keys():
        res += shapes_and_dim[number][0] * shapes_and_dim[number][1]
    return res

def A_cent(shapes_and_dim):
    # Multiple each area in the cross section with the distances between their centroids and the bottom
    res = 0
    for number in shapes_and_dim.keys():
        res += shapes_and_dim[number][0] * shapes_and_dim[number][1] * shapes_and_dim[number][2]
    return res

def y_bar(shapes_and_dim):
    # Calculate height of the centroidal axis
    return A_cent(shapes_and_dim)/A_tot(shapes_and_dim)

def I_calc(shapes_and_dim):
    # Calculate the second moment of area of the centroidal axis
    res = 0
    cent_axis = y_bar(shapes_and_dim)
    for number in shapes_and_dim.keys():
        # Calculate second moments of area using bh^3/12
        res += (shapes_and_dim[number][0]**3) * (shapes_and_dim[number][1]) / 12
        # Use parallel axis theorem Ad^2 and account for the distance 
        # between centroids of rectangles and centroid of overall cross-sectional area
        res += shapes_and_dim[number][0] * shapes_and_dim[number][1] * (shapes_and_dim[number][2]-cent_axis)**2
    return res

def Q_calc(shapes_and_dim, depth_of_int):
    # Calculate value of Q at the depth of interest of a cross sectional area
    res = 0
    # Take out areas that fall beneath the depth of interest (areas of interest)
    areas_of_int = {}
    for number in shapes_and_dim.keys():
        if shapes_and_dim[number][3]+shapes_and_dim[number][0] <= depth_of_int:
            # Depth of interest lies above the rectangle
            areas_of_int[number] = shapes_and_dim[number]
        elif shapes_and_dim[number][3] < depth_of_int:
            # Depth of interest cuts through some rectangles
            # Modified entries before adding to into the new areas_of_int dictionary
            en1 = depth_of_int-shapes_and_dim[number][3] 
            en2 = shapes_and_dim[number][1]
            en3 = shapes_and_dim[number][3] + en1/2
            en4 = shapes_and_dim[number][3]
            areas_of_int[number] = [en1, en2, en3, en4]
    y_cent = y_bar(shapes_and_dim)  # centroid of the original cross-sectional area
    for number in areas_of_int.keys():
        res += areas_of_int[number][0] * areas_of_int[number][1] * abs(y_cent - areas_of_int[number][2]) # Calculate Q = Ad
    return res

def Q_calc_glue(shapes_and_dim, depth_of_int):
    # Calculate value of Q at glue location 
    # (isolate areas from top down instead of bottom up)
    res = 0
    # Take out areas that lie above the glue location (areas of interest)
    areas_of_int = {}
    for number in shapes_and_dim.keys():
        if shapes_and_dim[number][3] >= depth_of_int:
            # Glue location lies beneath the rectangle
            areas_of_int[number] = shapes_and_dim[number]
        elif shapes_and_dim[number][3]+shapes_and_dim[number][0] > depth_of_int:
            # Glue height cuts through some rectangles
            # Modified entries before adding to into the new areas_of_int dictionary
            en1= shapes_and_dim[number][3]+shapes_and_dim[number][0]-depth_of_int
            en2=shapes_and_dim[number][1]
            en3=depth_of_int+en1/2
            en4=depth_of_int
            areas_of_int[number] = [en1, en2, en3, en4]
    y_cent = y_bar(shapes_and_dim) # centroid of the original cross-sectional area
    for number in areas_of_int.keys():
        res += areas_of_int[number][0] * areas_of_int[number][1] * abs(areas_of_int[number][2]-y_cent) # Calculate Q = Ad
    return res

def h_calc(x_profile, h_profile):
    '''
    Return the height portoflio along the bridge length if the bridge design has variable heights
    x_profile[] represents where variable heights occur 
    h_profile[] represents the heights at the x locations
    For example, x = [0, 200, 400] and h = [160, 200, 200] represent how h switches to 160 at x = 200
    '''
    # [0, 600, 600, 700, 700, 1270], [120, 120, 160, 160, 120, 120]
    h=np.zeros(x_profile[len(x_profile)-1]+1)
    for i in range(0, len(x_profile)):
        h[x_profile[i]] = h_profile[i]
        if i != len(x_profile) - 1 and x_profile[i+1] - x_profile[i] != 0:
            dh = (h_profile[i+1] - h_profile[i]) / (x_profile[i+1] - x_profile[i])
            for j in range(x_profile[i] + 1, x_profile[i + 1]):
                h[j] = h[j-1]+dh
    return h
                
def I_profile(h_calc, shapes_and_dim):
    '''
    Use the height portfolio to generate the second moment of inertia portfolio along x-direction
    '''
    I_along_x=np.zeros(len(h_calc))
    for i in range(0, len(h_calc)):
        modified_dim = copy.deepcopy(shapes_and_dim)
        for number in modified_dim.keys():
            if modified_dim[number][0] != 1.27:
                modified_dim[number][0] += (h_calc[i]-h_calc[0])
            if modified_dim[number][3] > 1.27:
                modified_dim[number][3] += (h_calc[i]-h)
                modified_dim[number][2] = modified_dim[number][3] + modified_dim[number][0]/2
        I_along_x[i] = I_calc(modified_dim)
    return I_along_x

def Q_profile(h_calc, shapes_and_dim):
    '''
    Version of Q calculation that is applicable to designs 3 and 5 
    or when bridge changes height midway, the centroidal axis shifts up/ down
    Use the height portfolio to generate Q portfolio along x-direction
    '''
    Q_along_x=np.zeros(len(h_calc))
    for i in range(0, len(h_calc)):
        modified_dim = copy.deepcopy(shapes_and_dim)
        for number in modified_dim.keys():
            if modified_dim[number][0] != 1.27:
                modified_dim[number][0] += (h_calc[i]-h_calc[0])
            if modified_dim[number][3] > 1.27:
                modified_dim[number][3] += (h_calc[i]-h)
                modified_dim[number][2] = modified_dim[number][3] + modified_dim[number][0]/2
        Q_along_x[i] = Q_calc(modified_dim, y_bar(modified_dim))
    return Q_along_x

def Q_profile_glue(h_calc, shapes_and_dim, depth_of_int):
    '''
    Version of Q calculation that is applicable to designs 3 and 5 
    or when bridge changes height midway, the glue location is shifted up by how much bridge height changes
    Use the height portfolio to generate Q portfolio along x-direction
    '''
    Q_along_x=np.zeros(len(h_calc))
    for i in range(0, len(h_calc)):
        modified_dim = copy.deepcopy(shapes_and_dim)
        for number in modified_dim.keys():
            if modified_dim[number][0] != 1.27:
                modified_dim[number][0] += (h_calc[i]-h)
            if modified_dim[number][3] > 1.27:
                modified_dim[number][3] += (h_calc[i]-h)
                modified_dim[number][2] = modified_dim[number][3] + modified_dim[number][0]/2
        Q_along_x[i] = Q_calc_glue(modified_dim, depth_of_int+(h_calc[i]-h_calc[0]))
    return Q_along_x

def y_bar_profile(h_calc, shapes_and_dim):
    '''This function calculates centroidal axis height profile along the length of the bridge.
    When the overall bridge length increases, y_bar should be shifted as well.
    Function first modifies the dictionary containing the dimensions, 
    and calculates the centoidal axis height of the new modified dictionary'''
    y_bar_along_x=np.zeros(len(h_calc))
    for i in range(0, len(h_calc)):
        modified_dim = copy.deepcopy(shapes_and_dim)
        for number in modified_dim.keys():
            if modified_dim[number][0] != 1.27:
                modified_dim[number][0] += (h_calc[i]-h_calc[0])
            if modified_dim[number][3] > 1.27:
                modified_dim[number][3] += (h_calc[i]-h)
                modified_dim[number][2] = modified_dim[number][3] + modified_dim[number][0]/2
        y_bar_along_x[i] = y_bar(modified_dim)
    return y_bar_along_x

####################################################################################################################
'''
The general steps included in the functions for calculating FOS of different failure modes
1. Obtain the relevant parameters for calculation (I, Q, b), which could change if the bridge is of a variable height
2. Calculate the stress (sigma) or shear stress (tau) profile along the bridge at every location given V_max and M_max
For matboard shear failure, shear stress is only calculated at the centroidal axis because it is where it's maximized
3. Take the given strengths of matboard and glue (max stress/ shear stress before failure) or calculate
the max strengths of matboard and glue before failure specifically for the buckling cases
4. Calculate FOS based on the ratio of failure stress and actual stress -> return value
5. Calculate the max moment and shear force portfolio that a given bridge design is capable to withstand -> return value
'''

def F_tension(M_max, I, y_bar):
    # Calculate factor of safety associated with matboard tensile failure
    # and maximum moment the bridge can resist given the tensile failure 
    if not isinstance(I, np.ndarray): 
       # If I is not passed into the function an array (I is constant throughout bridge)
       # Create an array that matches dimension of M_max array
        I_arr = np.linspace(I, I, 1201) 
        y_bar = np.linspace(y_bar, y_bar, 1201)
    else:
        I_arr = I
    # Lowest FOS along the bridge length
    sigma_bot = M_max*y_bar/I_arr 
    print(f"The maximum flexural tensile stress is {max(abs(sigma_bot))} MPa")           
    FOS = abs(sigma_t/sigma_bot)
    min_FOS = min(FOS)
    # Returns array of maximum moment that can be withheld by the beam
    M_fail = -1* sigma_t * I_arr / y_bar
    # OR
    # M_fail = FOS * M_max
    # Return allowable stress along length (array), and min FOS (a number)
    return M_fail, min_FOS

def F_compression(M_max, shapes_and_dim, I):
    # Calculate factor of safety associated with matboard compressive failure
    # and maximum moment the bridge can resist given the compressive failure 
    if not isinstance(I, np.ndarray): 
       # If I is not passed into the function an array (I is constant throughout bridge)
       # Create an array that matches dimension of M_max array
       I_arr = np.linspace(I, I, 1201)
    else:
        I_arr = I
    global h
    # Lowest FOS along the bridge length
    sigma_top = M_max*(h-y_bar(shapes_and_dim))/I_arr
    print(f"The maximum flexural compressive stress is {max(abs(sigma_top))} MPa")           
    FOS = abs(sigma_c/sigma_top)
    min_FOS = min(FOS)
    # Returns array of maximum moment that can be withheld by the beam
    M_fail = -1* sigma_c * I_arr / (h-y_bar(shapes_and_dim))
    # OR
    # M_fail = FOS * M_max
    # Return allowable stress along length (array), and min FOS (a number)
    return M_fail, min_FOS

def F_shear_cent(V_max, I, Q, b):
    # Calculate factor of safety associated with matboard shear failure
    # and maximum moment the bridge can resist given the shear failure 
    # Calculate the shear at centroid of cross-sectional area because it will be max shear
    if not isinstance(I, np.ndarray): 
       # If I is not passed into the function an array (I is constant throughout bridge)
       # Create an array that matches dimension of M_max array
        I_arr = np.linspace(I, I, 1201)
        Q_arr = np.linspace(Q, Q, 1201)
    else:
        I_arr = I
        Q_arr = Q
    b_arr = np.linspace(b, b, 1201)
    tau = abs(V_max*Q_arr/I_arr/b_arr)
    print(f"The maximum matboard shear stress is {max(tau)} MPa")           
    FOS = tau_m/tau
    min_FOS = min(FOS)
    V_fail = tau_m*I_arr*b_arr/Q_arr
    # OR
    # V_fail = FOS * V_max
    return V_fail, min_FOS

def F_shear_glue(V_max, I, Q, b):
    # Calculate factor of safety associated with contact cement shear failure
    # and maximum shear force the bridge can resist given the shear failure 
    if not isinstance(I, np.ndarray): 
       # If I is not passed into the function an array (I is constant throughout bridge)
       # Create an array that matches dimension of V_max array
        I_arr = np.linspace(I, I, 1201)
        Q_arr = np.linspace(Q, Q, 1201)
    else:
        I_arr = I
        Q_arr = Q
    b_arr = np.linspace(b, b, 1201)
    tau = abs(V_max*Q_arr/I_arr/b_arr)
    print(f"The maximum glue shear stress is {max(tau)} MPa")           
    FOS = tau_g/tau
    min_FOS = min(FOS)
    V_fail = tau_g*I_arr*b_arr/Q_arr
    # OR
    # V_fail = FOS * V_max
    return V_fail, min_FOS

####################################################################################################################

def F_flange_buckle(M_max, shapes_and_dim, I, b_flange, t):
    # Calculate factor of safety associated with horizontal flange buckling
    # and maximum moment the bridge can resist given the buckling failure 
    global h
    if not isinstance(I, np.ndarray): 
        I_arr = np.linspace(I, I, 1201) 
    else:
        I_arr = I
    sigma_fail = 4*math.pi**2*E_m/12/(1-nu_m**2)*(t/b_flange)**2
    sigma = abs(M_max*(h-y_bar(shapes_and_dim))/I_arr)
    print(f"The maximum stress to cause flexural buckling is {max(sigma)} MPa")    
    print(f"The flexural buckling failure stress is {sigma_fail} MPa")                  
    FOS = abs(sigma_fail/sigma)
    min_FOS = min(FOS)
    M_fail = -1*sigma_fail * I_arr / (h-y_bar(shapes_and_dim))

    return M_fail, min_FOS

def F_tip_buckle(M_max, shapes_and_dim, I, b_tip):
    # Calculate factor of safety associated with tip buckling (three constrained, one free side)
    # and maximum moment the bridge can resist given the buckling failure 
    global h, t
    if not isinstance(I, np.ndarray): 
       I_arr = np.linspace(I, I, 1201) 
    else:
        I_arr = I
    sigma_fail = 0.425*math.pi**2*E_m/12/(1-nu_m**2)*(t/b_tip)**2
    sigma = abs(M_max*(h-y_bar(shapes_and_dim))/I_arr)
    print(f"The maximum stress to cause tip buckling is {max(sigma)} MPa") 
    print(f"The tip buckling failure stress is {sigma_fail} MPa")                     
    FOS = abs(sigma_fail/sigma)
    min_FOS = min(FOS)
    M_fail = -1*sigma_fail * I_arr / (h-y_bar(shapes_and_dim))

    return M_fail, min_FOS

def F_web_flex_shear_buckle(M_max, shapes_and_dim, I):
    # Calculate factor of safety associated with bridge web buckling (vertical)
    # and maximum moment the bridge can resist given the buckling failure 
    global h, t
    if not isinstance(I, np.ndarray): 
       I_arr = np.linspace(I, I, 1201)
    else:
        I_arr = I 
    sigma_fail = 6*math.pi**2*E_m/12/(1-nu_m**2)*(t/(h-y_bar(shapes_and_dim)))**2
    sigma = abs(M_max*(h-y_bar(shapes_and_dim))/I_arr)
    print(f"The maximum stress to cause web flexural buckling is {max(sigma)} MPa")   
    print(f"The web flexural buckling failure stress is {sigma_fail} MPa")                             
    FOS = abs(sigma_fail/sigma)
    min_FOS = min(FOS)
    M_fail = -1*sigma_fail * I_arr / (h-y_bar(shapes_and_dim))

    return M_fail, min_FOS

def F_web_shear_buckle(V_max, diaph_loc, I, Q, b, h):
    # Calculate factor of safety associated with shear buckling associated with diaphragms
    # and maximum shear force the bridge can resist given the buckling failure 
    global t
    diaph_dist = [0] * (len(diaph_loc)-1)
    for i in range(0, len(diaph_dist)):
        diaph_dist[i] = diaph_loc[i+1] - diaph_loc[i]
    
    b_arr = np.linspace(b, b, 1201)
    if not isinstance(I, np.ndarray): 
       I_arr = np.linspace(I, I, 1201) 
       Q_arr = np.linspace(Q, Q, 1201)
       h_arr = np.linspace(h, h, 1201)
    else:
        I_arr = I 
        Q_arr = Q
        h_arr = h

    a = np.zeros(len(V_max))
    for i in range(0, len(diaph_loc)-1):
        a[diaph_loc[i]: diaph_loc[i+1]] = diaph_dist[i]

    tau_fail = 5*(math.pi)**2*E_m/12/(1-nu_m**2)*((t/h_arr)**2 + (t/a)**2)
    # tau_fail = 5*(math.pi)**2*E_m/12/(1-nu_m**2)*((t/75)**2 + (t/400)**2)
    # tau_fail = np.linspace(tau_fail, tau_fail, 1201)
    tau = abs(V_max*Q_arr/I_arr/b_arr)
    print(f"The maximum stress to cause web shear buckling is {max(tau)} MPa") 
    print(f"The web shear buckling failure stress is {min(tau_fail)} MPa")                     

    FOS = abs(tau_fail/tau)
    min_FOS = min(FOS)
    V_fail = tau_fail*I_arr*b_arr/Q_arr

    V_fail[0 : diaph_loc[0]] = None
    V_fail[diaph_loc[len(diaph_loc)-1] : len(tau_fail)-1] = None

    return V_fail, min_FOS
