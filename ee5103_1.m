%a0260045r
clc
clear all

N = 10000;   
n = 1;    
sw = 0.1;    
sv = 1;      
A = 1;       
C = 1;
R1 = sw^2;   
R2 = sv^2;   
P_n(:, 1) = 1e6;   
w = random('normal', 0, sw, n, N);  
v = random('normal', 0, sv, n, N);  

for k = 1:N
    Kf(:, k) = P_n(:, k) * C' * (C * P_n(:, k) * C' + R2)^(-1);
    K(:, k) = (A * P_n(:, k) * C') * (C * P_n(:, k) * C' + R2)^(-1);
    P(:, k) = P_n(:, k) - (P_n(:, k) * C') * (C * P_n(:, k) * C' + R2)^(-1) * C * P_n(:, k);
    P_n(:, k + 1) = A * P_n(:, k) * A' + R1 - K(:, k) * (C * P_n(:, k) * C' + R2) * K(:, k)';
end

for run = 1:n
    x_m(:, 1) = 0;   
    x(:, 1) = 5;    
    for k = 1:N
        x(:, k + 1) = A * x(:, k) + w(run, k);
        y(k) = C * x(:, k) + v(run, k);
        xh(:, k) = x_m(:, k) + Kf(:, k) * (y(k) - C * x_m(:, k));
        x_m(:, k + 1) = A * x_m(:, k) + K(:, k) * (y(k) - C * x_m(:, k));
    end
end

figure;  
subplot(3, 1, 1);
plot(0:N-1, x(1:N), 'b-','DisplayName', '$x(k)$'), hold on;
plot(0:N-1, xh(1:N), 'r:', 'DisplayName', '$\hat{x}(k|k)$'), hold on;
xlabel('k')
legend('Interpreter', 'latex')
title('Graph 1')

subplot(3, 1, 2); 
plot(0:N-1, P, 'bo','DisplayName', '$P(k|k)$'), hold on;
xlabel('k')
legend('Interpreter', 'latex')
title('Graph 2')

subplot(3, 1, 3); 
plot(0:N-1, Kf,'bo', 'DisplayName', '$K_{f}(k)$'), hold on;
xlabel('k')
legend('Interpreter', 'latex')
title('Graph 3')

Bias = sum(x(1:N) - xh) / (N)
Var = sum((x(1:N) - xh).^2) / (N)



