
clear; clc; close all;
setenv('TEMP','C:\Users\Eeshan Malhotra\Temp');
rng(5);


%% PARAMETERS
n = 200;
k = 10; % k-sparse signal
r = 0.1;
tolerance = 1e-6; % difference in consecutive estimates of x
epsilon = 1e-4; % for optimization threshold
epsilon_noise = epsilon;
maxiter = 200;


%% GENERATE DATA

A = normrnd(5, 100, [n, n]);
B = normrnd(5, 100, [n, n]);

x_actual = zeros(n,1);
x_actual(randsample(n, k)) = normrnd(0, 1, [k, 1]);
beta_actual = r*(2*rand(n,1) - 1); 
delta_actual = diag(beta_actual);
phi_actual = (A + (B * delta_actual));
y_actual = phi_actual * x_actual;

%noise = zeros(n, 1);
noise = normrnd(0,1,[n,1]);
noise = noise * epsilon_noise / norm(noise, 2);

y_measured = y_actual + noise;


%% ESTIMATE
tic

beta_estimated = zeros(n,1);
x_estimated = zeros(n,1);

rel_diff_x = Inf;
R = r * ones(n,1);
iter = 0;

list_rel_diff_x = zeros(maxiter, 1);
list_norm_diff_x = zeros(maxiter, 1);


while rel_diff_x > tolerance && iter < maxiter
    iter = iter + 1
    
    %Update phi
	delta_estimated = diag(beta_estimated);
    phi_estimated = (A + (B * delta_estimated));

    %% 1. Find x
    disp('Optimizing x');
    prev_x_estimated = x_estimated;
    
    cvx_begin
        variable x_estimated(n)
        minimize(norm(x_estimated,1))
        subject to
            norm(y_measured - phi_estimated*x_estimated, 2) <= epsilon
    cvx_end
    
    norm_diff_x = norm(x_estimated - prev_x_estimated, 1);
    rel_diff_x = abs(norm(x_estimated,1) - norm(prev_x_estimated, 1)) / norm(prev_x_estimated, 1);
    
    list_rel_diff_x(iter)  = rel_diff_x;
    list_norm_diff_x(iter) = norm_diff_x;
    
    %% 2. Find beta (assuming x fixed)
    disp('Optimizing beta');
    prev_beta_estimated = beta_estimated;
    
    cvx_begin
        variable beta_estimated(n)
        minimize(norm(y_measured - ((A + B*diag(beta_estimated)) * x_estimated), 2))
        subject to
            -R <= beta_estimated <= R
    cvx_end
    
    diff_beta = norm(beta_estimated - prev_beta_estimated, 1);
    

end

delta_estimated = diag(beta_estimated);
phi_estimated = (A + (B * delta_estimated));

toc

iter

%% EVALUATE RESULTS

close all;
Error_x = (x_actual - x_estimated);
norm(Error_x, 1);

figure;
plot(1:n, x_actual, '-r', 1:n, x_estimated, '-b')
legend('actual', 'recovered', 'error');

figure
plot(x_actual - x_estimated);

figure;
plot(x_actual, x_estimated, '+b')

figure;
plot(list_norm_diff_x, '-b');

figure;
plot(list_rel_diff_x, '-r');



%If beta(or equivalently, phi) were  known exactly:

x_estimated_LS = phi_actual \ y_measured;

x_estimated_O =  zeros(n,1); %Oracle

cvx_begin
    variable x_estimated_O(n)
    minimize(norm(x_estimated_O,1))
    subject to
        norm(y_measured - phi_actual*x_estimated_O, 2) <= epsilon
cvx_end

figure;
plot(1:n, x_actual, '-r', 1:n, x_estimated, '-b', 1:n, x_estimated_O, '-g');
legend('actual', 'recovered', 'known phi');


    
norm(x_actual - x_estimated, 1)    % AA-
norm(x_actual - x_estimated_O, 1)  % O-  (=known phi)
    
%save('results/env_tol=1e-6_eps=1e-1_noise=0');
    
    
    
    
    
    
    
    
    