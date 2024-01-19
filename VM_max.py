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
    # P_train = np.ones(6) * (P/6) # Load case 1
    P_train = np.array([1.35/2/3.35, 1.35/2/3.35, 1/2/3.35, 1/2/3.35, 1/2/3.35, 1/2/3.35]) * P # locomotive first car
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
max_moment = 10000
max_moment_ind = -1
max_shear = 0
max_shear_ind = -1
for i in range(1201):
    for j in range(2056):
        if M_all[j][i] < max_moment:
            max_moment = M_all[j][i]
            max_moment_ind = j
        if V_all[j][i] > max_shear:
            max_shear = V_all[j][i]
            max_shear_ind = j

print(f"The max moment which equals to {max_moment} occurs at location {max_moment_ind} mm")
print(f"The max shear which equals to {max_shear} occurs at location {max_shear_ind} mm")


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
    figure, axis = plt.subplots(1, 2)
    axis[0].plot([0, 0], [0,  V_all[max_shear_ind][0]], color='#1f77b4')
    axis[0].plot(Vx, V_all[max_shear_ind], label="SFD at location of max shear", color='#1f77b4')
    axis[0].plot(Vx, V_max, label="Max shear at all locations on the bridge", color='red')
    axis[0].plot([0, 0], [0,  V_max[0]], color='red')
    axis[0].set_title("Max Shear Stress Throughout Bridge", fontsize = 10)
    axis[1].plot(Vx, M_all[max_moment_ind], label="BMD at location of max moment" ,color='#1f77b4')
    axis[1].plot(Vx, M_max, label = "Max moment at all locations on the bridge",color='red')
    axis[1].set_title("Max Bending Moment Throughout Bridge", fontsize = 10)

    axis[0].legend(fontsize = 7, loc="lower right")
    axis[1].legend(fontsize = 7, loc="lower right")
    axis[0].set_xlabel("x (mm)")
    axis[0].set_ylabel("V (N)")
    axis[1].set_xlabel("x (mm)")
    axis[1].set_ylabel("M (Nmm)")
    figure.tight_layout()
    plt.show()
