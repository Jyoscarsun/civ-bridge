import numpy as np
import matplotlib.pyplot as plt 
import os
import bridge_calc

L = 1200
n = L # Graphing 0-1200 m along length of the bridge
P = 400 # Basic load case
x = np.arange(0, n+1, 1) # x-coordinates 

x_train = np.array([52, 228, 392, 568, 732, 908])
# difference between wheels of the train (contact points of the train and matboard top flange)
x_train_diff = np.array([0, 176, 340, 516, 680, 856]) 

# There are 2057 possible train locations (i.e front wheel at index 0 to last wheel at index 1201)
V_all = np.zeros([2057, 1201])
M_all = np.zeros([2057, 1201])
# Maximum shear force and moment profile
V_max = np.zeros(1201)
M_max = np.zeros(1201)


def fail_load(FOS):
    # Given all FOS's against all modes of failure, calculate the P (dominated by the smallest FOS)
    global P
    return P * min(FOS)

def FOS_2_load(FOS):
    # Practically, my team expects the bridge to fail at half of the theoretical calculated load,
    # Which requires a factor of safety of 2
    global P
    FOS_adjusted = []
    for e in FOS:
        FOS_adjusted.append(e / 2)
    return P * min(FOS_adjusted)

for i in range(0,2057):
    # Depending on the load case, uncomment the actual load case
    P_train = np.ones(6) * (P/6) # Load case 1
    # P_train = np.array([1.35/2/3.35, 1.35/2/3.35, 1/2/3.35, 1/2/3.35, 1/2/3.35, 1/2/3.35]) * P # locomotive first car
    # P_train = np.array([1/2/3.35, 1/2/3.35, 1/2/3.35, 1/2/3.35, 1.35/2/3.35, 1.35/2/3.35]) * P # locomotive last car


    # Matching distance between each part of the train and starting location specified in i
    x_train_loc = np.linspace(i, i, 6) - x_train_diff
    valid_load = np.logical_and(x_train_loc > 0, x_train_loc < n) # load that is actually present on the bridge
    x_train_loc = x_train_loc[valid_load]
    P_train = P_train[valid_load]
    moments = x_train_loc * P_train
    

    # Reaction force at supports
    By = sum(moments)/L
    Ay = sum(P_train) - By

    # Calculate shear force (non-cumulative), 
    # i.e only locations of train wheels and supports will have non zero values
    V_train = np.zeros(P_train.size+2) # shear force in the beam
    Vx_train = np.zeros(P_train.size+2) # location of shear force in the beam
    V_train[0] = Ay
    V_train[1:1+P_train.size] = np.flip(-1*P_train)
    V_train[P_train.size+1] = By
    Vx_train[0] = 0
    Vx_train[1:1+P_train.size] = np.flip(x_train_loc)
    Vx_train[P_train.size+1] = L
    
    # Cumulative shear force (equivalent to the up and down method)
    V_train_cul = np.zeros(V_train.size)
    for k in range(V_train_cul.size-1, -1, -1):
        V_train_cul[k] = sum(V_train[0:k+1])


    V_sum = 0
    V = np.zeros(L+1)
    Vx = np.arange(0, 1201, 1)
    index = 0
    while(index < Vx_train.size):
        V[int(Vx_train[index])] = V_train_cul[index]
        if index+1 < Vx_train.size:
            V[int(Vx_train[index]+1):int(Vx_train[index+1])] = V_train_cul[index]
        index += 1
    
    V_all[i,:]=V
    
    # Moment calculations
    M = np.zeros(L+1)
    for k in range(0, V.size-1, 1):
        M[k] = -1*sum(V[0:k+1])
    M_all[i,:] = M

# Calculate the max possible moment and shear force at every point location
for i in range(1201):
    M_max[i] = min(M_all[:,i])
    pos_max = 0
    neg_min = 0
    for k in range(0, 1201):
        if V_all[k][i] > pos_max:
            pos_max = V_all[k][i]
        elif V_all[k][i] < neg_min:
            neg_min = V_all[k][i]
    if abs(pos_max) >= abs(neg_min):
        V_max[i] = pos_max
    else:
        V_max[i] = neg_min


if __name__ == "__main__":
    # Cross-sectional designs are divided into rectangles and passed into a dictionary
    # The parameters in lists in the dictionary are 
    # [height, width, distance from centroid of rectangle to bottom, distance from bottom of rectangle to bottom]

    # Design 0: provided in the assignment
    design0 = {
        1: [1.27, 80, 0.635, 0],
        2: [72.46, 1.27, 37.5, 1.27],
        3: [72.46, 1.27, 37.5, 1.27],
        4: [1.27, 6.27, 74.365, 73.73],
        5: [1.27, 6.27, 74.365, 73.73],
        6: [1.27, 100, 75.635, 75]
    }
    y_bar0 = bridge_calc.y_bar(design0)
    print(f"The height of centroidal axis of design 0 cross section is {y_bar0}")
    I0 = bridge_calc.I_calc(design0)
    print(f"The second moment of area of design 0 cross section is {I0}")
    Q0c = bridge_calc.Q_calc(design0, y_bar0)
    print(f"The Q value at centroid of cross section equals to {Q0c}")
    Q0g = bridge_calc.Q_calc_glue(design0, 75)
    print(f"The Q value at glue height of cross section equals to {Q0g}")
    diaph0 = [435, 835]

    # # Design 1: I-shaped cross section
    # design1 = {
    #     1: [1.27, 100, 0.635, 0],
    #     2: [1.27, 42.54, 1.905, 1.27],
    #     3: [154.92, 1.27, 80, 2.54], 
    #     4: [154.92, 1.27, 80, 2.54],  
    #     5: [1.27, 30, 158.095, 157.46], 
    #     6: [1.27, 30, 158.095, 157.46],  
    #     7: [1.27, 100, 159.365, 158.73] 
    # }
    # I1 = bridge_calc.I_calc(design1)
    # y_bar1 = bridge_calc.y_bar(design1)
    # Q1c = bridge_calc.Q_calc(design1, y_bar1)
    # Q1g = bridge_calc.Q_calc_glue(design1, 158.73)

    # # Design 2: Cross-section with multiple middle U-shaped webs
    # design2 = {
    #     1: [1.27, 100, 0.635, 0],
    #     2: [1.27, 30, 1.905, 1.27],
    #     3: [1.27, 30, 1.905, 1.27],
    #     4: [1.27, 30, 1.905, 1.27],
    #     5: [34.92, 1.27, 20, 2.54],
    #     6: [34.92, 1.27, 20, 2.54],
    #     7: [34.92, 1.27, 20, 2.54],
    #     8: [34.92, 1.27, 20, 2.54],
    #     9: [34.92, 1.27, 20, 2.54],
    #     10: [34.92, 1.27, 20, 2.54],
    #     11: [1.27, 30, 38.095, 37.46],
    #     12: [1.27, 30, 38.095, 37.46],
    #     13: [1.27, 30, 38.095, 37.46],
    #     14: [1.27, 100, 39.365, 38.73]
    # }
    # I2 = bridge_calc.I_calc(design2)
    # y_bar2 = bridge_calc.y_bar(design2)
    # Q2c = bridge_calc.Q_calc(design2, y_bar2)
    # Q2g = bridge_calc.Q_calc_glue(design2, 1.27)

    # # Design 3: Trapezoidal Bridge Shape
    # design3 = {
    #     1: [1.27, 80, 0.635, 0],
    #     2: [116.19, 1.27, 59.365, 1.27],
    #     3: [116.19, 1.27, 59.365, 1.27],
    #     4: [1.27, 30, 118.095, 117.46],
    #     5: [1.27, 30, 118.095, 117.46],
    #     6: [1.27, 100, 119.365, 118.73]
    # }
    # h3 = bridge_calc.h_calc([0, 420, 850, 1270], [120, 200, 200, 120])
    # I3 = bridge_calc.I_profile(h3, design3)
    # start = int((len(I3)-(L+1))/2)
    # end = int((len(I3)-(L+1))/2+L)
    # I3 = I3[start : end+1]
    # Q3c = bridge_calc.Q_profile(h3, design3)
    # Q3g = bridge_calc.Q_profile_glue(h3, design3, 118.73)
    # Q3c = Q3c[start : end+1]
    # Q3g = Q3g[start : end+1]
    # y_bar3 = bridge_calc.y_bar_profile(h3, design3)
    # y_bar3 = y_bar3[start : end+1]
    # diaph3 = [50, 150, 270, 420, 563, 707, 850, 1000, 1120, 1220]
    # diaph3 = [x - start for x in diaph3]

    # # Design 4: Box Girder With Constant Height
    # design4_org = {
    #     1: [1.25, 85, 0.635, 0],
    #     2: [117.46, 1.25, 60, 1.27],
    #     3: [117.46, 1.25, 60, 1.27],
    #     4: [1.27, 25, 118.095, 117.46],
    #     5: [1.27, 25, 118.095, 117.46],
    #     6: [1.27, 100, 119.35, 118.73]
    # }
    # design4 = {
    #     1: [1.25, 85, 0.635, 0],
    #     2: [117.46, 1.25, 60, 1.27],
    #     3: [117.46, 1.25, 60, 1.27],
    #     4: [1.27, 25, 118.095, 117.46],
    #     5: [1.27, 25, 118.095, 117.46],
    #     6: [1.27, 100, 119.35, 118.73]
    # }
    # I4 = bridge_calc.I_calc(design4)
    # y_bar4 = bridge_calc.y_bar(design4)
    # Q4c = bridge_calc.Q_calc(design4, y_bar4)
    # Q4g = bridge_calc.Q_calc_glue(design4, 118.73)
    # diaph4 = [50, 150, 270, 420, 563, 707, 850, 1000, 1120, 1220]
    # diaph4 = [x - 35 for x in diaph4]

    # Design 5: Box Girder With Step Heights
    # design5 = {
    #     1: [1.25, 85, 0.635, 0],
    #     2: [117.46, 1.25, 60, 1.27],
    #     3: [117.46, 1.25, 60, 1.27],
    #     4: [1.27, 15, 118.095, 117.46],
    #     5: [1.27, 15, 118.095, 117.46],
    #     6: [1.27, 100, 119.35, 118.73]
    # }
    # h5 = bridge_calc.h_calc([0, 420, 420, 850, 850, 1270], [120, 120, 160, 160, 120, 120])
    # I5 = bridge_calc.I_profile(h5, design5)
    # start = int((len(I5)-(L+1))/2)
    # end = int((len(I5)-(L+1))/2+L)
    # I5 = I5[start : end+1]
    # Q5c = bridge_calc.Q_profile(h5, design5)
    # Q5g = bridge_calc.Q_profile_glue(h5, design5, 118.73)
    # Q5c = Q5c[start : end+1]
    # Q5g = Q5g[start : end+1]
    # y_bar5 = bridge_calc.y_bar_profile(h5, design5)
    # y_bar5 = y_bar5[start : end+1]
    # diaph5 = [50, 150, 270, 420, 528, 635, 743, 850, 1000, 1120, 1220]
    # diaph5 = [x - start for x in diaph5]

    # # Design 6: C-Shaped Webs
    # design6 = {
    #     1: [1.27, 100, 0.635, 0],
    #     2: [1.27, 10, 1.905, 1.27], 
    #     3: [1.27, 10, 1.905, 1.27],  
    #     4: [157.46, 1.27, 80, 1.27],  
    #     5: [157.46, 1.27, 80, 1.27],  
    #     6: [1.27, 10, 158.095, 157.46], 
    #     7: [1.27, 10, 158.095, 157.46], 
    #     8: [1.27, 100, 159.365, 158.73] 
    # }
    # I6 = bridge_calc.I_calc(design6)
    # y_bar6 = bridge_calc.y_bar(design6)
    # Q6c = bridge_calc.Q_calc(design6, y_bar6)
    # Q6g = bridge_calc.Q_calc_glue(design6, 158.73)
    # diaph6 = [50, 150, 270, 420, 528, 635, 743, 850, 1000, 1120, 1220]
    # diaph6 = [x - 35 for x in diaph6]

    # Design 7
    # design7 = {
    #     1: [1.27, 76.27, 0.635, 0],
    #     2: [1.27, 10, 1.905, 1.27],
    #     3: [1.27, 10, 1.905, 1.27],
    #     4: [114.92, 1.27, 60, 2.54],
    #     5: [114.92, 1.27, 60, 2.54],
    #     6: [1.27, 38.135, 118.095, 117.46],
    #     7: [1.27, 38.135, 118.095, 117.46],
    #     8: [1.27, 100, 119.365, 118.73]
    # }
    # I7 = bridge_calc.I_calc(design7)
    # y_bar7 = bridge_calc.y_bar(design7)
    # Q7c = bridge_calc.Q_calc(design7, y_bar7)
    # Q7g = bridge_calc.Q_calc(design7, 1.27)
    # diaph7 = [50, 150, 270, 420, 528, 635, 743, 850, 1000, 1120, 1220]
    # diaph7 = [x - 35 for x in diaph7]

    # Design 8
    # design8 = {
    #     1: [1.27, 76.27, 0.635, 0],
    #     2: [114.92, 1.27, 58.73, 1.27],
    #     3: [114.92, 1.27, 58.73, 1.27],
    #     4: [1.27, 10, 116.825, 116.19],
    #     5: [1.27, 10, 116.825, 116.19],
    #     6: [1.27, 76.27, 118.095, 117.46],
    #     7: [1.27, 100, 119.365, 118.73]
    # }
    # I8 = bridge_calc.I_calc(design8)
    # y_bar8 = bridge_calc.y_bar(design8)
    # Q8c = bridge_calc.Q_calc(design8, y_bar8)
    # Q8g = bridge_calc.Q_calc_glue(design8, 117.46)
    # diaph8 = [50, 150, 270, 420, 528, 635, 743, 850, 1000, 1120, 1220]
    # diaph8 = [x - 35 for x in diaph8]

    ####################################################################################################################

    figure, axis = plt.subplots(2, 3)
    # The first row are the shear force diagrams (matboard shear failure, glue shear failure, shear buckling failure)
    for i in range(0, 3):
        axis[0][i].axhline(y=0, color='black', linestyle='-', linewidth=1)
        axis[0][i].plot([0, 0], [0, V_all[i][0]], color='#1f77b4')
        axis[0][i].set_xlim(0, 1200)
        axis[0][i].set_xlabel("x (mm)")
        axis[0][i].set_ylabel("V (N)")
        axis[1][i].set_xlim(0, 1200)
        axis[1][i].axhline(y=0, color='black', linestyle='-', linewidth=1)
        axis[1][i].set_xlabel("x (mm)")
        axis[1][i].set_ylabel("M (Nmm)")

    # Graph shear force and moment envelope
    for i in range(1201):
        axis[0][0].plot(Vx, V_all[i], color='#1f77b4')
        axis[0][1].plot(Vx, V_all[i], color='#1f77b4')
        axis[0][2].plot(Vx, V_all[i], color='#1f77b4')
        axis[1][0].plot(Vx, M_all[i], color='red')
        axis[1][1].plot(Vx, M_all[i], color='red')
        axis[1][2].plot(Vx, M_all[i], color='red')

    # axis[0].set_title("Shear Forces V(x) Along Bridge", fontsize = 10)
    # axis[1].set_title("Bending Moment M(x) Along Bridge", fontsize = 10)

    # Calculate FOS associated with different types of failure modes
    M_t, FOS_t = bridge_calc.F_tension(M_max, I0, y_bar0)
    M_c, FOS_c = bridge_calc.F_compression(M_max, design0, I0)
    V_m, FOS_m = bridge_calc.F_shear_cent(V_max, I0, Q0c, 2.54)
    V_g, FOS_g = bridge_calc.F_shear_glue(V_max, I0, Q0g, 75)
    M_fb, FOS_fb = bridge_calc.F_flange_buckle(M_max, design0, I0, 77.46, 1.27)
    M_tb, FOS_tb = bridge_calc.F_tip_buckle(M_max, design0, I0, 10)
    M_wfb, FOS_wfb = bridge_calc.F_web_flex_shear_buckle(M_max, design0, I0)
    V_wsb, FOS_wsb = bridge_calc.F_web_shear_buckle(V_max, diaph0, I0, Q0c, 2.54, 75)


    # Print out the calculated FOS
    print(f"The FOS of flexural tensile failure is {FOS_t}.")
    print(f"The FOS of flexural compressive failure is {FOS_c}.")
    print(f"The FOS of shear matboard failure is {FOS_m}.")
    print(f"The FOS of shear glue failure is {FOS_g}.")
    print(f"The FOS of flange buckling is {FOS_fb}.")
    print(f"The FOS of tip buckling is {FOS_tb}.")
    print(f"The FOS of web flexural buckling is {FOS_wfb}.")
    print(f"The FOS of web shear buckling is {FOS_wsb}.")
    FOS = [FOS_t, FOS_c, FOS_m, FOS_g, FOS_fb, FOS_tb, FOS_wfb, FOS_wsb]
    FOS = [FOS_t, FOS_c, FOS_m, FOS_g, FOS_fb, FOS_wfb]
    print(f"The theoretical failure load is {fail_load(FOS)}")
    # print(f"The experimental failure load (2.0 FOS) is {FOS_2_load(FOS)}")

    # Graph out the max shear force moment a FOS can resist on the correct SFD and BMD
    axis[0][0].plot(Vx, V_m, label="Matboard Shear Failure", color = "orange")
    axis[0][0].plot(Vx, -1*V_m, color = "orange")
    
    axis[0][1].plot(Vx, V_g, label="Glue Shear Failure", color = "orange")
    axis[0][1].plot(Vx, -1*V_g, color = "orange")

    axis[0][2].plot(Vx, V_wsb, label="Web Shear Buckling Failure", color = "orange")
    axis[0][2].plot(Vx, -1*V_wsb, color = "orange")

    axis[1][0].plot(Vx, M_t, label="Matboard Tension Failure")
    axis[1][0].plot(Vx, M_c, label="Matboard Compression Failure")
    axis[1][1].plot(Vx, M_fb, label="Matboard Buckling Failure - Mid")
    axis[1][1].plot(Vx, M_tb, label="Top Flange Tip Buckling Failure")
    axis[1][2].plot(Vx, M_wfb, label="Matboard Buckling Failure - Webs")


    axis[0][0].legend(fontsize = 7, loc="lower right")
    axis[0][1].legend(fontsize = 7, loc="lower right")
    axis[0][2].legend(fontsize = 7, loc="lower right")
    axis[1][0].legend(fontsize = 7, loc="lower right")
    axis[1][1].legend(fontsize = 7, loc="lower right")
    axis[1][2].legend(fontsize = 7, loc="upper right")
    figure.tight_layout()
    plt.show()
