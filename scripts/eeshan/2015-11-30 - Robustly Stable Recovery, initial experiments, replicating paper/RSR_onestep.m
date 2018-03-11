
clear; clc; 
rng(5);

%% PARAMETERS
n = 200;
m = 80;
k = 10;
r = 0.1;
tolerance = 1e-6; % difference in consecutive estimates of x
epsilon = 1e-2; % for optimization threshold
epsilon_noise = epsilon;

%% GENERATE DATA

A = normrnd(2, 5, [m, n]);
B = normrnd(2, 5, [m, n]);

beta_actual = r*(2*rand(n,1) - 1); 
delta_actual = diag(beta_actual);
phi_actual = (A + (B * delta_actual));

x_actual = zeros(n,1);
x_actual(randsample(n, k)) = normrnd(0, 1, [k, 1]);

y_actual = phi_actual * x_actual;

%noise = zeros(n, 1);
noise = normrnd(0, 1, [m,1]);
noise = noise * epsilon_noise / norm(noise, 2);

y_measured = y_actual + noise;


%% ESTIMATE

beta_estimated = zeros(n,1);
x_plus = zeros(n,1);
x_minus = zeros(n,1);
p = zeros(n,1);
M = [A -A B];
R = r * ones(n,1);

cvx_begin
    variables x_plus(n) x_minus(n) p(n)
    minimize(sum(x_plus + x_minus))

    subject to
        norm(y_measured - M*[x_plus; x_minus; p], 2) <= epsilon

        x_plus >= 0
        x_minus >= 0
        -r*(x_plus + x_minus) <= p <= r*(x_plus + x_minus) 
cvx_end

assert(all(x_plus .* x_minus < 1e-8), 'CHECK FAIL: x+ and x- overlap')
x_estimated = x_plus - x_minus;
beta_estimated(x_estimated > 1e-8) = p(x_estimated > 1e-8) ./ x_estimated(x_estimated > 1e-8); % Avoid div by 0 locations


%% EVALUATE RESULTS

figure;
plot(1:n, x_actual, '-r', 1:n, x_estimated, '-b');
legend('actual', 'recovered');

norm(x_actual - x_estimated, 1)    % AA-

    

beep; pause(1); beep; pause(1); beep; 
