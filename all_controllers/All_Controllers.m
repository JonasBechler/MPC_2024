% All Controllers
close all
clear



A = [1 -1.0646 0.3679];
B = [0.7585 -0.4552];
plant_delay = 5;

disp("plant poles:")
disp(roots(A))
disp("plant zeros:")
disp(roots(B))


C = [1 -1];

%step response
D = [1 -1]; 
% ramp response
%D = [1 -2 1];

disp("disturbance poles:")
disp(roots(D))
disp("disturbance zeros:")
disp(roots(C))

[E, F] = deconv_n(C, conv(A, D), plant_delay);

disp("E:")
disp(E)
disp("F:")
disp(F)


figure(1)

% ----------------------------
% Minimal Variance Control
subplot(3, 1, 1)



S = conv(B, conv(E, D));
R = F;
T = C;

sim("All_Controllers_Sim.slx");

t = tout;
y = yout(:,1);
u = yout(:,2);
y_star = yout(:,3);
noise = yout(:,4);


plot(t, y, "-")
hold on
plot(t, u, "-")
plot(t, y_star, "*")
plot(t, noise, "-")
hold off
title('Minimal Variance Control')
legend('y', 'u', 'y*', 'noise')
xlabel('Time (s)')
ylabel('Value')
ylim([-5, 15])
grid on

% ----------------------------
% One-Step Ahead Predictive Control
% K: A*D + B(1)/input_weight * B <= must be stable
subplot(3, 1, 2)

input_weight = 0.1;

K1 = conv(A, D);
K2 = B(1)/input_weight * B;

if length(K1) < length(K2)
    K1 = [zeros(1, length(K2) - length(K1)), K1];
else
    K2 = [zeros(1, length(K1) - length(K2)), K2];
end

K = K1 + K2;
disp("Roots of A*D + B(1)/input_weight * B:")
disp(" - must be stable")
roots(K)



S1 = B(1)/input_weight * conv(B, conv(E, D));
S2 = conv(C,D);

if length(S1) < length(S2)
    S1 = [zeros(1, length(S2) - length(S1)), S1];
elseif length(S2) < length(S1)
    S2 = [zeros(1, length(S1) - length(S2)), S2];
end

S = S1 + S2;
R = B(1)/input_weight * F;
T = B(1)/input_weight * C;

sim("All_Controllers_Sim.slx");

t = tout;
y = yout(:,1);
u = yout(:,2);
y_star = yout(:,3);
noise = yout(:,4);

plot(t, y, "-")
hold on
plot(t, u, "-")
plot(t, y_star, "*")
plot(t, noise, "-")
hold off
title('One-Step Ahead Predictive Control')
legend('y', 'u', 'y*', 'noise')
xlabel('Time (s)')
ylabel('Value')
ylim([-5, 15])
grid on

% ----------------------------
% One-Step Ahead Predictive Control
% with input frequency weighting
subplot(3, 1, 3)

H = [1, -2, 1]
W = conv(D, [1, -1]);

S1 = conv(C, W);
S2 = conv(H, conv(E, conv(B, H)));

if length(S1) < length(S2)
    S1 = [zeros(1, length(S2) - length(S1)), S1];
elseif length(S2) < length(S1)
    S2 = [zeros(1, length(S1) - length(S2)), S2];
end

S = S1 + S2;
R = conv(H, F);
T = conv(H, C);

sim("All_Controllers_Sim.slx");

t = tout;
y = yout(:,1);
u = yout(:,2);
y_star = yout(:,3);
noise = yout(:,4);

plot(t, y, "-")
hold on
plot(t, u, "-")
plot(t, y_star, "*")
plot(t, noise, "-")
hold off
title('One-Step Ahead Predictive Control with input frequency weighting')
legend('y', 'u', 'y*', 'noise')
xlabel('Time (s)')
ylabel('Value')
ylim([-5, 15])
grid on




























function [C, D] = deconv_n(A, B, n)
    % A / B = C + q^-n * D / B

    len_A = length(A);
    len_B = length(B);
    max_len = max(max(len_A, len_B), n);

    if len_A < max_len
        A = [A zeros(1, max_len - len_A)];
    end

    [x, r] = deconv(A(1:len_B), B);

    A = [r, A(len_B+1:end)];
    
    if n == 1
        C = x;
        D = [0, A(2:end)];
    else
        [C_next, D_next] = deconv_n(A(2:end), B, n - 1);
        C = [x, C_next];
        D = [0, D_next];
    end
    
end
