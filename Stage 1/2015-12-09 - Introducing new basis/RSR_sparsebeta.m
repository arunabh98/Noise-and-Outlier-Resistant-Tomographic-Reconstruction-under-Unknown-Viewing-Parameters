
clear; clc; close all;
setenv('TEMP','C:\Users\Eeshan Malhotra\Temp');
rng(1);


%% PARAMETERS
m = 80;
n = 200;
kX = 40; % k-sparse signal
kBeta = 60;
r = 0.1;
tolerance = 1e-6; % difference in consecutive estimates of x
epsilon = 1e-2; % for optimization threshold
epsilon_noise = epsilon;
maxiter = 100;


%% GENERATE DATA

A = normrnd(2, 2, [m, n]);
B = normrnd(2, 2, [m, n]);
U = dctmtx(n);

theta_actual = zeros(n,1);
theta_actual(randsample(n, kX)) = normrnd(0, 1, [kX, 1]);
beta_actual = zeros(n,1);
beta_actual(randsample(n, kBeta))= r*(2*rand(kBeta,1) - 1);

delta_actual = diag(beta_actual);
phi_actual = (A + (B * delta_actual));
y_actual = phi_actual * U * theta_actual;

%noise = zeros(m, 1);
noise = normrnd(0,0.5,[m,1]);
noise = noise * epsilon_noise / norm(noise, 2);

y_measured = y_actual + noise;


%% ESTIMATE
tic

beta_estimated = zeros(n,1);
beta_estimated = r*(2*rand(n,1) - 1);
theta_estimated = zeros(n,1);

rel_diff_theta = Inf;
R = r * ones(n,1);
iter = 0;

list_rel_diff_theta = zeros(maxiter, 1);
list_norm_diff_theta = zeros(maxiter, 1);
list_theta_error = zeros(maxiter, 1);
list_beta_error = zeros(maxiter, 1);
E1 = zeros(maxiter, 1);
E2 = zeros(maxiter, 1);


%while (rel_diff_theta > tolerance || rel_diff_beta > tolerance) && iter < maxiter
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
    
    E1(iter) = norm(theta_estimated, 1);
    
    %% 2. Find beta (assuming x fixed)
    disp('Optimizing beta');
    prev_beta_estimated = beta_estimated;
    
    cvx_begin
        variable beta_estimated(n)
        minimize norm(beta_estimated, 1)
        subject to
            norm(y_measured - ((A + B*diag(beta_estimated)) * U * theta_estimated), 2) <= epsilon
            -R <= beta_estimated <= R
    cvx_end

     
%     cvx_begin
%         variable beta_estimated(n)
%         minimize norm(y_measured - ((A + B*diag(beta_estimated)) * U * theta_estimated), 2)
%         subject to
%             -R <= beta_estimated <= R
%     cvx_end


    norm_diff_beta = norm(beta_estimated - prev_beta_estimated, 1);
    rel_diff_beta = abs(norm(beta_estimated,1) - norm(prev_beta_estimated, 1)) / norm(prev_beta_estimated, 1);
   
    
    list_beta_error(iter) = norm(beta_actual - beta_estimated, 1);
    
    E2(iter) = norm(beta_estimated, 1);
%     
%     if mod(iter, 100) == 0 || iter == 1
%         set = [num2str(iter) 'iters' '_m=' num2str(m) '_kX=' num2str(kX) '_kBeta=' num2str(kBeta)];
%         save(['env_AA_' set]);
%         
%         plot(1:n, theta_actual, '-r', 1:n, theta_estimated, '-b')
%         legend('actual', 'recovered');
%         print(['actual_vs_recovered_' set], '-dpng')
%         
%         plot(list_theta_error(1:iter), '-r');
%         print(['thetaerror_' set], '-dpng')
%         
%         plot(1:n, beta_actual, '-r', 1:n, beta_estimated, '-b')
%         legend('beta actual', 'beta recovered');
%         print(['beta actual_vs_recovered_' set], '-dpng')
%         
%         plot(list_theta_error(1:iter), '-r');
%         print(['thetaerror_' set], '-dpng')
%         
%         
%         
% 
%     end
    

end


        set = [num2str(iter) 'iters' '_m=' num2str(m) '_kX=' num2str(kX) '_kBeta=' num2str(kBeta)];
        save(['env_AA_' set]);
        
        plot(1:n, theta_actual, '-r', 1:n, theta_estimated, '-b')
        legend('actual', 'recovered');
        print(['actual_vs_recovered_' set], '-dpng')
        
        plot(list_theta_error(1:iter), '-r');
        print(['thetaerror_' set], '-dpng')
        
        plot(1:n, beta_actual, '-r', 1:n, beta_estimated, '-b')
        legend('beta actual', 'beta recovered');
        print(['beta actual_vs_recovered_' set], '-dpng')
        
        plot(list_theta_error(1:iter), '-r');
        print(['thetaerror_' set], '-dpng')
        




delta_estimated = diag(beta_estimated);
phi_estimated = (A + (B * delta_estimated));

toc

iter
save('env_full.mat');
%% EVALUATE RESULTS
% 
figure;
plot(1:n, theta_actual, '-r', 1:n, theta_estimated, '-b')
legend('actual', 'recovered');


figure;
plot(1:kBeta, beta_actual(beta_actual~=0), '-r', 1:kBeta, beta_estimated(beta_actual~=0), '-b')
legend('actual', 'recovered');


figure;
plot(E1 + E2);
 
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
    
    
    
    
    
    
    
    
    