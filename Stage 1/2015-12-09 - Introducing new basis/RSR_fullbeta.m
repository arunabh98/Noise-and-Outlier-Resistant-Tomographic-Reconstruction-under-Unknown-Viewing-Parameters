
clear; clc; close all;
setenv('TEMP','C:\Users\Eeshan Malhotra\Temp');
rng(1);


%% PARAMETERS
m = 80;
n = 200;
k = 40; % k-sparse signal
r = 0.1;
tolerance = 1e-6; % difference in consecutive estimates of x
epsilon = 1e-2; % for optimization threshold
epsilon_noise = epsilon;
maxiter = 1000;


%% GENERATE DATA

A = normrnd(5, 5, [m, n]);
B = normrnd(5, 5, [m, n]);
U = dctmtx(n);

theta_actual = zeros(n,1);
theta_actual(randsample(n, k)) = normrnd(0, 1, [k, 1]);
beta_actual = r*(2*rand(n,1) - 1); 
delta_actual = diag(beta_actual);
phi_actual = (A + (B * delta_actual));
y_actual = phi_actual * U * theta_actual;

%noise = zeros(m, 1);
noise = normrnd(0,1,[m,1]);
noise = noise * epsilon_noise / norm(noise, 2);

y_measured = y_actual + noise;


%% ESTIMATE
tic

beta_estimated = zeros(n,1);
theta_estimated = zeros(n,1);

rel_diff_theta = Inf;
R = r * ones(n,1);
iter = 0;

list_rel_diff_theta = zeros(maxiter, 1);
list_norm_diff_theta = zeros(maxiter, 1);
list_theta_error = zeros(maxiter, 1);

while rel_diff_theta > tolerance && iter < maxiter
    iter = iter + 1
    
    %Update phi
	delta_estimated = diag(beta_estimated);
    phi_estimated = (A + (B * delta_estimated));

    %% 1. Find x
    disp('Optimizing theta');
    prev_theta_estimated = theta_estimated;
    
    cvx_begin
        variable theta_estimated(n)
        minimize(norm(theta_estimated,1))
        subject to
            norm(y_measured - phi_estimated*U*theta_estimated, 2) <= epsilon
    cvx_end
    
    norm_diff_theta = norm(theta_estimated - prev_theta_estimated, 1);
    rel_diff_theta = abs(norm(theta_estimated,1) - norm(prev_theta_estimated, 1)) / norm(prev_theta_estimated, 1);
    
    list_rel_diff_theta(iter)  = rel_diff_theta;
    list_norm_diff_theta(iter) = norm_diff_theta;
    list_theta_error(iter) = norm(theta_actual - theta_estimated, 1);
    
    %% 2. Find beta (assuming x fixed)
    disp('Optimizing beta');
    prev_beta_estimated = beta_estimated;
    
    cvx_begin
        variable beta_estimated(n)
        minimize(norm(y_measured - ((A + B*diag(beta_estimated)) * U * theta_estimated), 2))
        subject to
            -R <= beta_estimated <= R
    cvx_end
    
    diff_beta = norm(beta_estimated - prev_beta_estimated, 1);
    
    if mod(iter, 100) == 0 || iter == 1
        set = [num2str(iter) 'iters' '_m=' num2str(m) '_k=' num2str(k)];
        save(['env_AA_' set]);
        
        figure;
        plot(1:n, theta_actual, '-r', 1:n, theta_estimated, '-b')
        legend('actual', 'recovered');
        print(['actual_vs_recovered_' set], '-dpng')
        
        figure;
        plot(list_theta_error(1:iter), '-r');
        print(['thetaerror_' set], '-dpng')

    end
    

end

delta_estimated = diag(beta_estimated);
phi_estimated = (A + (B * delta_estimated));

toc

iter
save('env_full.mat');
%% EVALUATE RESULTS
% 
% close all;
% Error_theta = (theta_actual - theta_estimated);
% norm(Error_theta, 1)
% 
% figure;
% plot(1:n, theta_actual, '-r', 1:n, theta_estimated, '-b')
% legend('actual', 'recovered');
% 
% figure
% plot(theta_actual - theta_estimated);
% 
% figure;
% plot(theta_actual, theta_estimated, '+b')
% 
% figure;
% plot(list_norm_diff_theta(10:end), '-b');
% 
% figure;
% plot(list_theta_error, '-r');
% 
% % 
% 
% %If beta(or equivalently, phi) were  known exactly:
% 
% x_estimated_LS = phi_actual \ y_measured;
% 
% x_estimated_O =  zeros(n,1); %Oracle
% 
% cvx_begin
%     variable x_estimated_O(n)
%     minimize(norm(x_estimated_O,1))
%     subject to
%         norm(y_measured - phi_actual*x_estimated_O, 2) <= epsilon
% cvx_end
% 
% figure;
% plot(1:n, theta_actual, '-r', 1:n, theta_estimated, '-b', 1:n, x_estimated_O, '-g');
% legend('actual', 'recovered', 'known phi');
% 
% 
%     
% norm(theta_actual - theta_estimated, 1)    % AA-
% norm(theta_actual - x_estimated_O, 1)  % O-  (=known phi)
%     
%save('results/env_tol=1e-6_eps=1e-1_noise=0');
    
    
    
    
    
    
    
    
    