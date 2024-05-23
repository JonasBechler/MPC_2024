import numpy as np
import matplotlib.pyplot as plt
from control.matlab import *

from deconv_n import deconv_n
from poly_div import poly_div

d = 1 # plant delay
A = [1, -1.2254, 0.5711, -0.3507, 0.005]
A = [1, -0.8]
B = [0.0281, 0.1278, 0.0513, 0.0013]
B = [0.4, 0.6]

plant_tf = tf(A, B, dt=True)
delay_poly = np.zeros(d + 1)
delay_poly[0] = 1
z = tf('z')
delay_tf = z**-d
pzmap(plant_tf*delay_tf)




C = [1]
D = [1, -1]


# h_p = prediction horizon
# h_i = initialization horizon
# h_c = control horizon
# l   = input weighting
h_p = np.zeros(4).astype(int)
h_i = np.zeros(4).astype(int)
h_c = np.zeros(4).astype(int)
l = np.zeros(4)

# property 1 (Minimal Variance Control)
# h_p = h_i = plant_delay + 1
# h_c = 1
# l = 0
h_p[0] = d + 1
h_i[0] = d + 1
h_c[0] = 1
l[0] = 0

# property 2 (One step Ahead Predictive Control)
# h_p = h_i = plant_delay + 1
# h_c = 1
# l != 0
h_p[1] = d + 1
h_i[1] = d + 1
h_c[1] = 1
l[1] = 1 # not 0


# property 3 (Pole placement)
# h_c = order(A) + order(D)
# h_i = order(B) + plant_delay + 1
# h_p > h_i + h_c
# l = 0
h_i[2] = len(B) + d
h_c[2] = len(A) + len(D) - 2
h_p[2] = h_i[2] + h_c[2] + 1 
l[2] = 0

# property 4 (Internal Model Control)
# h_p -> inf (100+)
# h_c = 1
# h_i = plant_delay + 1
# l = 0
# D = 1 - q^-1 
h_p[3] = 100
h_i[3] = d + 1
h_c[3] = 1
l[3] = 0



N = 3

# N1 = start prediction
# N2 = end prediction

h_i[:] = d+1
h_p[:] = d+N
h_c[:] = N
l[:] = 0.8



def compute_RST(d, A, B, C, D, h_i, h_c, h_p, l):
    #E_, F_ = deconv_n(C, np.polymul(A, D), d)
    #F_ = F_[d:]
    E_ges = []
    F_ges = []
    G_ges = []
    H_ges = []

    for i in range(0, h_c):
        E, F = poly_div(C, np.polymul(A, D), i + 1)
        E_ges.append(E)
        F_ges.append(F)

        G, H = poly_div(np.polymul(B, E), C, i + d + 1)
        G_ges.append(G)
        H_ges.append(H)

        pass

    print (f"E_ges: {E_ges}")
    print (f"F_ges: {F_ges}")
    print (f"G_ges: {G_ges}")
    print (f"H_ges: {H_ges}")


    psi_size = h_c
    psi = np.zeros((psi_size, psi_size))

    for i in range(0, psi_size): 

        psi[i, 0:i+1] = G_ges[i][i::-1]
    
    print(f"psi: {psi}")

    psi_t = np.transpose(psi)
    gamma = np.matmul(np.linalg.inv(np.matmul(psi_t, psi) + l * np.eye(psi_size)), psi_t)
    gamma = gamma[0, :]

    R = []
    for i in range(0, h_c):
        R = np.polyadd(R, np.polymul(gamma[i], F_ges[i]))

    S = []
    for i in range(0, h_c):
        S = np.polyadd(S, np.polymul(gamma[i], H_ges[i]))
    S = np.append(S, 0)
    S = np.polymul(D, np.polyadd(C, S))
    
    T = []
    for i in range(0, h_c):
        T = np.polyadd(T, gamma[i])
    T = np.polymul(C, T)

    #R = np.polymul(gamma, F)
    #S = np.polymul(D, np.polyadd(C, np.polymul(gamma, H_ges[0])))

    print(f"R: {R}")
    print(f"S: {S}")
    print(f"T: {T}")

    return R, S, T

def calculate_tf(A, B, R, S, T, plant_tf, delay_tf):

    p_c = 1 + tf(R, S, dt=True) * delay_tf * plant_tf

    y_ys = tf(np.polymul(B, T), 1, dt=True) / p_c
    y_vy = tf(np.polymul(A, S), 1, dt=True) / p_c               # Sensitivity function
    y_vu = tf(np.polymul(B, S), 1, dt=True) * delay_tf / p_c    # Sensitivity function x System
    y_n  = tf(np.polymul(B, R), 1, dt=True) * delay_tf / p_c    # Complementary sensitivity function

    u_ys = tf(np.polymul(A, T), 1, dt=True) / p_c
    u_vy = - tf(np.polymul(A, R), 1, dt=True) / p_c             # Sensitivity function x Controller
    u_vu = tf(np.polymul(A, S), 1, dt=True) / p_c               # Sensitivity function
    u_n  = - tf(np.polymul(A, R), 1, dt=True) / p_c             # Sensitivity function x Controller

    return {y_ys, y_vy, y_vu, y_n}, {u_ys, u_vy, u_vu, u_n}

plt.figure("Comparison of properties")
for i in range(1):
    [R, S, T] = compute_RST(d, A, B, C, D, h_i[i], h_c[i], h_p[i], l[i])
    [y_tf, u_tf] = calculate_tf(A, B, R, S, T, plant_tf, delay_tf)

    [y_ys, y_vy, y_vu, y_n] = y_tf
    [u_ys, u_vy, u_vu, u_n] = u_tf



    
    
    plt.subplot(4, 4, i+1)
    [poles, zeros] = pzmap(y_ys, plot=False)
    plt.plot(np.real(zeros), np.imag(zeros), 'bo')
    plt.plot(np.real(poles), np.imag(poles), 'rx')
    plt.grid()
    plt.title(f"Property {i+1}")
    plt.xlabel("Re")
    plt.ylabel("Im")
    plt.plot(np.cos(np.linspace(0, 2*np.pi, 100)), np.sin(np.linspace(0, 2*np.pi, 100)), '-')
    plt.legend(["Poles", "Zeros"])

    plt.subplot(4, 4, i+5)
    [y, t] = step(y_ys)
    plt.plot(t, y)
    plt.grid()
    plt.title("Step response")
    plt.xlabel("Time")
    plt.ylabel("Amplitude")
    

plt.tight_layout()
plt.show()



plt.subplot(2, 2, 2)
pzmap(ComplementarySensitivity_tf)

plt.subplot(2, 2, 3)
pzmap(SensitivityController_tf)

plt.subplot(2, 2, 4)
pzmap(SensitivityPlant_tf)
plt.show()
pass








