import numpy as np
import matplotlib.pyplot as plt
from control.matlab import *

from deconv_n import deconv_n
from poly_div import poly_div

d = 0 # plant delay
B = [0.0281, 0.1278, 0.0513, 0.0013]
B = [0.4, 0.6]
A = [1, -1.2254, 0.5711, -0.3507, 0.005]
A = [1, -0.8]

plt.figure("Model Predictive Control")
plt.subplot(4, 5, 1)
plant_tf = tf(B, np.append(A, np.zeros(d+1)), dt=True)
[poles, zeros] = pzmap(plant_tf, plot=False)
plt.plot(np.real(zeros), np.imag(zeros), 'bo')
plt.plot(np.real(poles), np.imag(poles), 'rx')
plt.grid()
plt.title(f"Plant Poles and Zeros")
plt.plot(np.cos(np.linspace(0, 2*np.pi, 100)), np.sin(np.linspace(0, 2*np.pi, 100)), '-')


plt.subplot(4, 5, 5+1)
[y, t] = step(plant_tf, T=12)
plt.plot(t, y)
plt.grid()
plt.title("Step response")
plt.tight_layout()
#plt.show()





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
h_i[2] = len(B) + d -1
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

h_i[0] = d
h_p[0] = d+N
h_c[0] = N
l[0] = 0.8



def compute_RST(d, A, B, C, D, h_i, h_c, h_p, l):
    #E_, F_ = deconv_n(C, np.polymul(A, D), d)
    #F_ = F_[d:]
    E_ges = []
    F_ges = []
    G_ges = []
    H_ges = []

    for i in range(h_i, h_p):
        E, F = poly_div(C, np.polymul(A, D), i + 1)
        E_ges.append(E)
        F_ges.append(F)

        G, H = poly_div(np.polymul(B, E), C, i + 1)
        G_ges.append(G)
        H_ges.append(H)

        pass

    print (f"E_ges: {E_ges}")
    print (f"F_ges: {F_ges}")
    print (f"G_ges: {G_ges}")
    print (f"H_ges: {H_ges}")


    psi_size = h_c
    psi = np.zeros((psi_size, psi_size))

    # for i in range(0, psi_size): 
    #     psi[i, 0:i+1] = G_ges[i][i::-1]
    

    # psi:
    # [g_0_0,       0,      0,      0]
    # [g_1_1,   g_1_0,      0,      0]
    # [g_2_2,   g_2_1,  g_2_0,      0]
    # [g_3_3,   g_3_2,  g_3_1,  g_3_0]






    for i in range(h_i, h_p):
        print(f"i: {i}")
        print("psi[i-h_i, 0:min(i+1, h_c)]")
        print(f"= psi[{i}-{h_i}, 0:min({i}+1, {h_c})]")
        print("= " + str(psi[i-h_i, 0:min(i+1, h_c)]))

        print("np.flip(G_ges[i-h_i][max(0, i - h_c + 1):i+1])")
        print(f"np.flip(G_ges[{i}-{h_i}][max(0, {i} - {h_c} + 1):{i}+1])")
        print("= " + str(np.flip(G_ges[i-h_i][max(0, i - h_c + 1):i+1])))

        psi[i-h_i, 0:min(i+1, h_c)] = np.flip(G_ges[i-h_i][max(0, i - h_c + 1):i+1])


    print(f"psi:\n{psi}")

    psi_t = np.transpose(psi)
    gamma = np.matmul(np.linalg.inv(np.matmul(psi_t, psi) + l * np.eye(psi_size)), psi_t)
    gamma = gamma[0, :]

    print(f"gamma: {gamma}")

    # R * y + S * u = T * y_s
    # see https://www.mathworks.com/help/sps/ref/rstcontroller.html#d126e252765

    # h_c = h_p - h_i
    # R = sum(j = [0, ..., h_c])(gamma_j * F_j)
    # correct as is
    # plus q^-j?
    R = []
    for i in range(0, h_c):
        R = np.polyadd(R, np.polymul(gamma[i], F_ges[i]))



    # is it q^-1 or q^-j?
    # S = D * { C + sum(j = [0, ..., h_c])(gamma_j * H_j-d) * q^-1}
    #   = D * C + D * sum(j = [0, ..., h_c])(gamma_j * H_j-d) * q^-1
    #   = S1 + S2

    S1 = np.flip(np.polymul(D, C))
    S2 = []
    for i in range(0, h_c):
        S2 = np.polyadd(S2, np.polymul(gamma[i], H_ges[i]))
    S2 = np.polymul(D, S2)
    S2 = np.flip(np.append(0, S2))

    #S = np.polymul(D, np.polyadd(np.append(C, 0), S))
    S = np.flip(np.polyadd(S1, S2))
    
    # T = C * sum(j = [0, ..., h_c])(gamma_j * q^-j)
    # T = C * gamma
    # correct as is
    T = []
    T = np.polymul(C, gamma)

    #R = np.polymul(gamma, F)
    #S = np.polymul(D, np.polyadd(C, np.polymul(gamma, H_ges[0])))

    print(f"R: {R}")
    print(f"S: {S}")
    print(f"T: {T}")

    return R, S, T

def calculate_tf(d, A, B, R, S, T):
    # p_c = A * S + B * R * q^(-d-1)
    p_c = np.polyadd(np.append(np.polymul(A, S), np.zeros(d+1)), np.polymul(B, R))
    #p_c = np.polyadd(np.append(np.polymul(np.polymul(A, S), D), np.zeros(d+1)), np.polymul(B, R))
    #p_c = np.polyadd(np.polymul(A, S), np.append(np.polymul(B, R), np.zeros(d+1)))
    #p_c = np.polyadd(np.append(np.polymul(A, S), np.zeros(d+1)), np.polymul(B, R))
    #p_c = np.polyadd(np.polymul(np.polymul(A, S), D), np.append(np.polymul(B, R), np.zeros(d+1)))
    #p_c = np.polyadd(np.polymul(A, S), np.append(np.polymul(B, R), np.zeros(d+1)))

    p_c1 = np.polymul(A, S)
    p_c2 = np.append(np.zeros(d+1), np.polymul(B, R))
    # make same length by appending zeros
    if len(p_c1) > len(p_c2):
        p_c2 = np.append(p_c2, np.zeros(len(p_c1) - len(p_c2)))
    else:
        p_c1 = np.append(p_c1, np.zeros(len(p_c2) - len(p_c1)))
    
    p_c = np.polyadd(p_c1, p_c2)
    


    
    #y_ys = tf(np.append(np.zeros(d+1), np.polymul(B, T)), p_c, dt=True)
    y_ys = tf(np.polymul(B, T), p_c, dt=True)
    y_ys = tf(np.polymul(B, T), np.append(p_c, np.zeros(d+1)), dt=True)

    #y_ys = tf(np.polymul(B, T), p_c, dt=True)
    y_vy = tf(np.polymul(A, S), p_c, dt=True)                               # Sensitivity function
    y_vu = tf(np.append(np.polymul(B, S), np.zeros(d+1)), p_c, dt=True)     # Sensitivity function x System
    y_n  = tf(np.append(np.polymul(B, R), np.zeros(d+1)), p_c, dt=True)     # Complementary sensitivity function

    u_ys =   tf(np.polymul(A, T), np.append(p_c, np.zeros(d+1)), dt=True)
    u_ys =   tf(np.polymul(A, T), p_c, dt=True)
    #u_ys =   tf(np.append(np.polymul(A, T), np.zeros(d+2)), p_c,  dt=True)
    u_vy = - tf(np.polymul(A, R), p_c, dt=True)                         # Sensitivity function x Controller
    u_vu =   tf(np.polymul(A, S), p_c, dt=True)                           # Sensitivity function
    u_n  = - tf(np.polymul(A, R), p_c, dt=True)                         # Sensitivity function x Controller

    return {y_ys, y_vy, y_vu, y_n}, {u_ys, u_vy, u_vu, u_n}

#plt.figure("Comparison of properties")
for i in range(0, 1):
    [R, S, T] = compute_RST(d, A, B, C, D, h_i[i], h_c[i], h_p[i], l[i])
    [y_tf, u_tf] = calculate_tf(d, A, B, R, S, T)

    [y_ys_tf, y_vy_tf, y_vu_tf, y_n_tf] = y_tf
    [u_ys_tf, u_vy_tf, u_vu_tf, u_n_tf] = u_tf



    
    
    plt.subplot(4,5, i+2)
    [poles, zeros] = pzmap(y_ys_tf, plot=False)
    plt.plot(np.real(zeros), np.imag(zeros), 'bo')
    plt.plot(np.real(poles), np.imag(poles), 'rx')
    plt.grid()
    plt.title(f"Poles and Zeros Y/Y*")
    
    plt.plot(np.cos(np.linspace(0, 2*np.pi, 100)), np.sin(np.linspace(0, 2*np.pi, 100)), '-')

    try:
        plt.subplot(4,5, i+2 +5)
        plt.grid()
        plt.title("Step response Y/Y* and U/U*")
        plt.xlabel("Time")
        plt.ylabel("Amplitude")
        [y_ys_step, t_step_y] = step(y_ys_tf, T=52)
        [u_ys_step, t_step_u] = step(u_ys_tf, T=52)
        plt.plot(t_step_y, y_ys_step, "r", t_step_u, u_ys_step, "b")
        #plt.plot(t, u_ys)
        pass
    except:
        pass



    
plt.show()



pass








