
clear;
clc;
close all;


d = 0; % plant delay
B = [0.0281, 0.1278, 0.0513, 0.0013];
B = [0.4, 0.6];

A = [1, -1.2254, 0.5711, -0.3507, 0.005];
A = [1, -0.8];


figure();
subplot(4, 5, 1);
plant_tf = tf(B, [A, zeros(1, d+1)], -1);
pzmap(plant_tf);
title("Plant Poles and Zeros");

subplot(4, 5, 5+1);
step(plant_tf, 12);
grid on
title("Step response");


C = [1];
D = [1, -1];

h_p = zeros(1, 4);
h_i = zeros(1, 4);
h_c = zeros(1, 4);
l = zeros(1, 4);

% property 1 (Minimal Variance Control)
% h_p = h_i = plant_delay + 1
% h_c = 1
% l = 0
h_p(1) = d + 1;
h_i(1) = d + 1;
h_c(1) = 1;
l(1) = 0;

% property 2 (One step Ahead Predictive Control)
% h_p = h_i = plant_delay + 1
% h_c = 1
% l != 0
h_p(2) = d + 1;
h_i(2) = d + 1;
h_c(2) = 1;
l(2) = 0.9; % not 0

% property 3 (Pole placement)
% h_c = order(A) + order(D)
% h_i = order(B) + plant_delay + 1
% h_p > h_i + h_c
% l = 0
h_c(3) = length(A) + length(D) - 2;
h_i(3) = length(B) + d;
h_p(3) = h_i(3) + h_c(3) + 1;
l(3) = 0.9;

% property 4 (Internal Model Control)
% h_p -> inf (100+)
% h_c = 1
% h_i = plant_delay + 1
% l = 0
% D = 1 - q^-1
h_p(4) = 100;
h_i(4) = d + 1;
h_c(4) = 100;
l(4) = 1.5;

N = 2;

% N1 = start prediction
% N2 = end prediction

% h_i(1) = d+1;
% h_p(1) = d+N+1;
% h_c(1) = N+1;
% l(1) = 0.8;

for i = 1:4

    E_ges = zeros(h_c(i), h_c(i));
    F_size = length(A) + length(D)-2;
    F_ges = zeros(h_c(i), F_size);
    G_ges = zeros(h_c(i), h_c(i));
    H_size = length(C);
    H_ges = zeros(h_c(i), H_size);

    for j = h_i(i):h_p(i)
        [E, F] = PolyDiv(C, conv(A, D), j);
        E_ges(j-h_i(i)+1, 1:length(E)) = E;
        F_ges(j-h_i(i)+1, 1:length(F)) = F;

        [G, H] = PolyDiv(conv(B, E), C, j);
        G_ges(j-h_i(i)+1, 1:length(G)) = G;
        H_ges(j-h_i(i)+1, 1:length(H)) = H;
    end

    disp("E_ges:")
    disp(E_ges)
    disp("F_ges:")
    disp(F_ges)
    disp("G_ges:")
    disp(G_ges)
    disp("H_ges:")
    disp(H_ges)

    % psi:
    % [g_0_0,       0,      0,      0]
    % [g_1_1,   g_1_0,      0,      0]
    % [g_2_2,   g_2_1,  g_2_0,      0]
    % [g_3_3,   g_3_2,  g_3_1,  g_3_0]
    
    psi = zeros(h_c(i), h_c(i));
    
    for j = h_i(i):h_p(i)
        psi(j-h_i(i)+1, min(j, h_c(i)):-1:1) = G_ges(j-h_i(i)+1, max(1, j - h_c(i) + 1):j);
    end
    
    disp("psi:")
    disp(psi)
    
    gamma = (psi' * psi + l(i) * eye(h_c(i))) \ psi';
    gamma = gamma(1, :);
    
    disp("gamma:")
    disp(gamma)
    
    
    % R * y + S * u = T * y_s
    % see https://www.mathworks.com/help/sps/ref/rstcontroller.html#d126e252765
    
    % R = sum(j = [0, ..., h_c])(gamma_j * F_j)
    R = zeros(1, length(F_ges(1, :)));
    for j = 1:h_c(i)
        R = R + gamma(j) * F_ges(j, :);
    end
    
    % S = D * { C + sum(j = [0, ..., h_c])(gamma_j * H_j-d) * q^-1}
    %   = D * C + D * sum(j = [0, ..., h_c])(gamma_j * H_j-d) * q^-1
    
    S1 = conv(D, C);
    S2 = zeros(1, length(H_ges(1,:)));
    for j = 1:h_c(i)
        S2 = S2 + gamma(j) * H_ges(j, :);
    end
    S2 = conv(D, S2);
    S2 = [0, S2];
    if length(S1) < length(S2)
        S1 = [S1, zeros(1, length(S2) - length(S1))];
    else
        S2 = [S2, zeros(1, length(S1) - length(S2))];
    end
    S = S1 + S2;
    
    % T = C * sum(j = [0, ..., h_c])(gamma_j * q^-j)
    T = C * gamma;
    
    disp("R:")
    disp(R)
    disp("S:")
    disp(S)
    disp("T:")
    disp(T)
    
    % p_c = A * S + B * R * q^(-d-1)
    p_c1 = conv(A, S);
    p_c2 = [zeros(1, d+1), conv(B, R)];
    if length(p_c1) < length(p_c2)
        p_c1 = [p_c1, zeros(1, length(p_c2) - length(p_c1))];
    else
        p_c2 = [p_c2, zeros(1, length(p_c1) - length(p_c2))];
    end
    p_c = p_c1 + p_c2;
    
    
    y_ys = tf(conv(B, T), [p_c, zeros(1, d+1)], -1);
    y_vy = tf(conv(A, S), p_c, -1);                               % Sensitivity function
    y_vu = tf([conv(B, S), zeros(1, d+1)], p_c, -1);     % Sensitivity function x System
    %y_n  = tf(S, p_c, -1);     % Complementary sensitivity function
    y_n  = tf([conv(B, R), zeros(1, d+1)], p_c, -1);     % Complementary sensitivity function
    
    u_ys =   tf(conv(A, T), p_c, -1);
    u_vy = - tf(conv(A, R), p_c, -1);                         % Sensitivity function x Controller
    % ? u_vy = - tf(conv(A, R), p_c, -1);                         % Sensitivity function x Controller
    u_vu =   tf(conv(B, S), p_c, -1);                           % Sensitivity function
    % ? u_vu =   tf(conv(A, S), p_c, -1);                           % Sensitivity function
    u_n  = - tf(conv(A, R), p_c, -1);                         % Sensitivity function x Controller
    % ? u_n  = - tf(conv(A, R), p_c, -1);                         % Sensitivity function x Controller
    
    subplot(4,5, i+1);
    pzmap(y_ys);
    title("Poles and Zeros Y/Y*");
    
    subplot(4,5, i+1 +5);
    step(y_ys, 12);
    hold on;
    step(u_ys, 12);
    grid on;
    title("Step response Y/Y* and U/Y*");
    xlabel("Time");
    ylabel("Amplitude");
    ylim([-0.2, 1.2]);
    
    
    
    
    
    
    
    waitforbuttonpress;
    





































end


















































































































































