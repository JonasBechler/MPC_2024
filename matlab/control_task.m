
clear;
clc;
close all;

plant = 1;

if plant == 1
    % Main Plant for Task
    d = 1; % plant delay
    B = [0.0281, 0.1278, 0.0513, 0.0013];
    A = [1, -1.2254, 0.5711, -0.3507, 0.005];
    C = [1 -0.9];
    D = [1, -1];

elseif plant == 2
    % Simple Plant Example from Camacho S56
    d = 0; % plant delay
    B = [0.4, 0.6];
    A = [1, -0.8];
    C = [1];
    D = [1, -1];
end


o_a = length(A) - 1;
o_b = length(B) - 1;


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

o_c = length(C) - 1;
o_d = length(D) - 1;







test_delay = false;
test_hc = true;

for i = 1:4
    
    if test_delay
        d = i-1;
        N = 5;
        hi = d + 1;
        hp = d + N + 1;
        hc = N + 1;
        lambda = 0.8;

    elseif test_hc
        d = 0;
        N = i*i*i
        hi = d + 1;
        hp = d + N + 1;
        hc = N+1;
        lambda = 0.8;


    else
        
        
        if i == 1
            % property 1 (Minimal Variance Control)
            % h_p = h_i = plant_delay + 1
            % h_c = 1
            % l = 0
            hp = d + 1;
            hi = d + 1;
            hc = 1;
            lambda = 0.00;
            
        elseif i == 2
            % property 2 (One step Ahead Predictive Control)
            % h_p = h_i = plant_delay + 1
            % h_c = 1
            % l != 0
            hp = d + 1;
            hi = d + 1;
            hc = 1;
            lambda = 0.8;
            
        elseif i == 3
            % property 3 (Pole placement)
            % h_c = order(A) + order(D)
            % h_i = order(B) + plant_delay + 1
            % h_p > h_i + h_c
            % l = 0
            hc = o_a + o_d;
            hi = o_b + d + 1;
            hp = hi + hc + 1;
            lambda = 0.00;
            
        elseif i == 4
            % property 4 (Internal Model Control)
            % h_p -> inf (100+)
            % h_c = 1
            % h_i = plant_delay + 1
            % l = 0
            % D = 1 - q^-1
            hp = 10;
            hc = d + 1;
            hi = d + 1;
            lambda = 0.00;
        end
        
    end
    
    disp("d: " + d + " hp: " + hp + " hi: " + hi + " hc: " + hc + " lambda: " + lambda);
    
    
    
    E_ges = zeros(hp, hp);
    F_size = o_a + o_d;
    F_ges = zeros(hp, F_size);
    G_ges = zeros(hp, hp);
    G_ges_flipped = zeros(hp, hp);
    H_size = length(C);
    H_ges = zeros(hp, H_size);
    
    for j = 1:hp
        [E, F] = PolyDiv(C, conv(A, D), j);
        E_ges(j, 1:length(E)) = E;
        F_ges(j, 1:length(F)) = F;
        
        [G, H] = PolyDiv(conv(B, E), C, j);
        G_ges(j, 1:length(G)) = G;
        H_ges(j, 1:length(H)) = H;

        G_ges_flipped(j, 1:length(G)) = fliplr(G);
    end
    
    disp("E_ges:")
    disp(E_ges)
    disp("F_ges:")
    disp(F_ges)
    disp("G_ges:")
    disp(G_ges)
    disp("H_ges:")
    disp(H_ges)

    % G_ges:

    % [g_0_0,       0,      0,      0]
    % [g_1_0,   g_1_1,      0,      0]
    % [g_2_0,   g_2_1,  g_2_2,      0]
    % [g_3_0,   g_3_1,  g_3_2,  g_3_3]

    % G_ges_flipped:
    
    % [g_0_0,       0,      0,      0]
    % [g_1_1,   g_1_0,      0,      0]
    % [g_2_2,   g_2_1,  g_2_0,      0]
    % [g_3_3,   g_3_2,  g_3_1,  g_3_0]
    
    % psi:
    % [g_1_1,   g_1_0,      0]
    % [g_2_2,   g_2_1,  g_2_0]
    % [g_3_3,   g_3_2,  g_3_1]

    
    psi = G_ges_flipped(hi:hp, 1:hc);
    
    % ?????????????????
    %for j = hi:hp
    %    buf = G_ges(j, :);
    %    buf = buf.
    %    % psi(j-hp+1, min(j, hp):-1:1) = G_ges(j, max(1, j - hp + 1):j);
    %end
    
    disp("psi:")
    disp(psi)
    
    gamma = (psi' * psi + lambda * eye(hc)) \ psi';
    gamma = gamma(1, :);
    
    disp("gamma:")
    disp(gamma)
    
    
    % R * y + S * u = T * y_s
    % see https://www.mathworks.com/help/sps/ref/rstcontroller.html#d126e252765
    
    % R = sum(j = [0, ..., h_c])(gamma_j * F_j)
    R = zeros(1, length(F_ges(1, :)));
    for j = hp:hp
        R = R + gamma(j - hp + 1) * F_ges(j, :);
    end
    disp("R:")
    disp(R)
    
    % S = D * { C + sum(j = [0, ..., h_c])(gamma_j * H_j-d) * q^-1}
    
    S1 = C;
    S2 = zeros(1, length(H_ges(1,:)));
    for j = hp:hp
        S2 = S2 + gamma(j - hp + 1) * H_ges(j - d, :);
    end
    S2 = [0, S2];
    
    if length(S1) < length(S2)
        S1 = [S1, zeros(1, length(S2) - length(S1))];
    else
        S2 = [S2, zeros(1, length(S1) - length(S2))];
    end
    S = S1 + S2;
    S = conv(D, S);
    
    disp("S:")
    disp(S)
    
    % T = C * sum(j = [0, ..., h_c])(gamma_j * q^-j)
    T = conv(C, gamma);
    disp("T:")
    disp(T)
    
    
    
    % p_c = A * S + B * R * q^(-d-1)
    p_c1 = conv(A, S);
    p_c2 = [zeros(1, d+1), conv(B, R)];
    disp("p_c1:")
    disp(p_c1)
    disp("p_c2:")
    disp(p_c2)
    if length(p_c1) < length(p_c2)
        p_c1 = [p_c1, zeros(1, length(p_c2) - length(p_c1))];
    else
        p_c2 = [p_c2, zeros(1, length(p_c1) - length(p_c2))];
    end
    p_c = p_c1 + p_c2;
    %p_c = p_c(1:find(p_c, 1, 'last'));

    disp("p_c:")
    disp(p_c)
    
    % y_ys = B * T * q^(-d-1) / p_c
    y_ys_n = [zeros(1, d+1), conv(B, T)];
    y_ys_d = p_c;
    if length(y_ys_n) < length(y_ys_d)
        y_ys_n = [y_ys_n, zeros(1, length(y_ys_d) - length(y_ys_n))];
    else
        y_ys_d = [y_ys_d, zeros(1, length(y_ys_n) - length(y_ys_d))];
    end
    disp("y_ys:")
    disp(y_ys_n)
    disp(y_ys_d)
    y_ys = tf(y_ys_n, y_ys_d, -1);

    % u_ys = A * T / p_c
    u_ys_n = conv(A, T);
    u_ys_d = p_c;
    if length(u_ys_n) < length(u_ys_d)
        u_ys_n = [u_ys_n, zeros(1, length(u_ys_d) - length(u_ys_n))];
    else
        u_ys_d = [u_ys_d, zeros(1, length(u_ys_n) - length(u_ys_d))];
    end
    disp("u_ys:")
    disp(u_ys_n)
    disp(u_ys_d)
    u_ys = tf(u_ys_n, u_ys_d, -1);




    AT = conv(A, T);
    AS = conv(A, S);
    BS = conv(B, S);
    AR = conv(A, R);
    BR = conv(B, R);

    y_vy = tf(conv(A, S), p_c, -1);                               % Sensitivity function
    y_vu = tf(conv(B, S), [p_c, zeros(1, d+1)], -1);     % Sensitivity function x System
    y_n  = tf(conv(B, R), [p_c, zeros(1, d+1)], -1);     % Complementary sensitivity function
    
    %u_ys =   tf(conv(A, T), p_c, -1);
    u_vy = - tf(conv(A, R), p_c, -1);                         % Sensitivity function x Controller
    u_vu =   tf(conv(B, R), [p_c, zeros(1, d+1)], -1);                           % Sensitivity function
    u_n  = - tf(conv(A, R), p_c, -1);                         % Sensitivity function x Controller
    
    subplot(4,5, i+1);
    pzmap(y_ys);
    title("Poles and Zeros Y/Y*");
    
    
    subplot(4,5, i+1 +5);
    step(y_ys);
    hold on;
    step(u_ys);
    grid on;
    title("Step response Y/Y* and U/Y*");
    xlabel("Time");
    ylabel("Amplitude");
    
    
    
    
end
