%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Predictive Control
% HTWG Konstanz
% Task 2024
% Jonas Bechler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 


clear;
clc;
close all;

plant = 2;

if plant == 1
    % Main Plant for Task
    d = 1; % plant delay
    B = [0.0281, 0.1278, 0.0513, 0.0013];
    A = [1, -1.2254, 0.5711, -0.3507, 0.05];
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
plant_tf = tf(B, A, -1, 'Variable', 'z^-1', 'InputDelay', d);
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
test_hc = false;

for i = 1:4
    
    if test_delay
        d = i-1;
        N = 500;
        hi = d + 1;
        hp = d + N + 1;
        hc = N + 1;
        lambda = 0.8;

    elseif test_hc
        d = 0;
        N = i
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
            hi = d+1;
            hp = d+1+2;
            hc = 1+2;
            lambda = 0.80;
            
        elseif i == 2
            % property 2 (One step Ahead Predictive Control)
            % h_p = h_i = plant_delay + 1
            % h_c = 1
            % l != 0
            hi = d + 1;
            hp = d + 1;
            hc = 1;
            lambda = 0.8;
            
        elseif i == 3
            % property 3 (Pole placement)
            % h_c = order(A) + order(D)
            % h_i = order(B) + plant_delay + 1
            % h_p > h_i + h_c
            % l = 0
            hi = o_b + d + 1;
            hc = o_a + o_d;
            hp = hi + hc+1
            lambda = 0.00;
            
        elseif i == 4
            % property 4 (Internal Model Control)
            % h_p -> inf (100+)
            % h_c = 1
            % h_i = plant_delay + 1
            % l = 0
            % D = 1 - q^-1
            hi = d + 1;
            hp = d + 1 + 50;
            hc = 50;
            lambda = 0.80;
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


    print_array("E_ges", E_ges)
    print_array("F_ges", F_ges)
    print_array("G_ges", G_ges)
    print_array("H_ges", H_ges)


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
    
    print_matrix("psi", psi)
    
    gamma = (psi' * psi + lambda * eye(hc)) \ psi';
    gamma = gamma(1, :);
    
    print_array("gamma", gamma)
    
    
    % R * y + S * u = T * y_s
    % see https://www.mathworks.com/help/sps/ref/rstcontroller.html#d126e252765
    
    % R = sum(j = [0, ..., h_c])(gamma_j * F_j)
    R = zeros(1, length(F_ges(1, :)));
    for j = hp:hp
        R = R + gamma(j - hp + 1) * F_ges(j, :);
    end
    print_array("R", R)

    
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
    
    print_array("S", S)
    
    % T = C * sum(j = [0, ..., h_c])(gamma_j * q^-j)
    T = conv(C, gamma);

    print_array("T", T)

    
    
    
    % p_c = A * S + B * R * q^(-d-1)
    p_c1 = conv(A, S);
    p_c2 = [zeros(1, d+1), conv(B, R)];
    print_array("p_c1", p_c1)
    print_array("p_c2", p_c2)

    if length(p_c1) < length(p_c2)
        p_c1 = [p_c1, zeros(1, length(p_c2) - length(p_c1))];
    else
        p_c2 = [p_c2, zeros(1, length(p_c1) - length(p_c2))];
    end
    p_c = p_c1 + p_c2;

    print_array("p_c", p_c)

    
    % y_ys = B * T * q^(-d-1) / p_c
    y_ys =  tf(conv(B, T), p_c, -1, 'InputDelay', d, 'Variable', 'z^-1');

    % y_vy = A * S / p_c
    y_vy =  tf(conv(A, S), p_c, -1, 'Variable', 'z^-1');     
    
    % y_vu = B * S * q^(-d-1) / p_c
    % Sensitivity function x System?
    y_vu =  tf(conv(B, S), p_c, -1, 'InputDelay', d, 'Variable', 'z^-1');     
    
    % y_n = B * R * q^(-d-1) / p_c
    % Complementary sensitivity function?
    y_n  =  tf(conv(B, R), p_c, -1, 'InputDelay', d, 'Variable', 'z^-1');     


    % u_ys = A * T / p_c
    u_ys =  tf(conv(A, T), p_c, -1, 'Variable', 'z^-1');

    % u_vy = - A * R / p_c
    % Sensitivity function x Controller?
    u_vy = -tf(conv(A, R), p_c, -1, 'Variable', 'z^-1');                      

    % u_vu = B * R / p_c
    % u_vu = B * R * q^(-d-1) / p_c
    % Sensitivity function?
    u_vu =  tf(conv(B, R), p_c, -1, 'InputDelay', d, 'Variable', 'z^-1');           
    
    % u_n = - A * R / p_c
    % Sensitivity function x Controller?
    u_n  = -tf(conv(A, R), p_c, -1, 'Variable', 'z^-1');
    
    p_f = deconv(p_c, C);
    if i == 1
        % property 1 test
        % p_f == 1/b0 * B
        disp("p_f == 1/b0 * B")
        disp(p_f)
        disp(1 / B(1) * B)
        
        
    elseif i == 2
        % property 2 test
        % no test for this property
        
    elseif i == 3
        % property 3 test
        % p_f == 1
        % p_c == C
        disp("p_f == 1")
        disp(p_f)
        disp(1)
        disp("p_c == C")
        disp(p_c)
        disp(C)
        
    elseif i == 4
        % property 4 test
        % p_f == A
        disp("p_f == A")
        disp(p_f)
        disp(A)

    end
    
    
    subplot(4,5, i+1);
    pzmap(y_ys);
    title("Poles and Zeros Y/Y*");
    
    
    subplot(4,5, i+1 +5);
    step(y_ys);
    hold on;
    step(u_ys);
    grid on;
    title("Step Response Y/Y* and U/U*");

end



function print_array(name, value)
    if length(value) < 10
        disp(name + ":");
        disp(value);
        return;
    end
end

function print_matrix(name, value)
    if length(value(1, :)) < 10
        disp(name + ":");
        disp(value);
        return;
    end
end
    
    

