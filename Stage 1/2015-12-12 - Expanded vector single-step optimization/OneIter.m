clear; clc; 
rng(5);

%% PARAMETERS

n = 200;
r = 0.1;                 
tolerance = 1e-6; % difference in consecutive estimates of x
epsilon = 1e-2; % for optimization threshold
epsilon_noise = epsilon;

numruns = 1;
m = 80;
kTheta = 10
kBeta = 200;

%% INITIALIZE

R = r * ones(n,1);
theta_error_averaged = 0;

tic
for run = 1:numruns;

    disp('====================================');
    disp(['              run ' num2str(run) '                 ']);
    theta_error = 0;

        %% GENERATE DATA

        A = normrnd(2, 2, [m, n]);
        B = normrnd(2, 2, [m, n]);
        U = dctmtx(n);

        theta_actual = zeros(n,1);
        beta_actual = zeros(n,1);
        
        noise = normrnd(0,1,[m,1]);
        noise = noise * epsilon_noise / norm(noise, 2);


        theta_actual_k = theta_actual;
        theta_actual_k(randsample(n, kTheta)) = normrnd(0, 1, [kTheta, 1]);

        beta_actual_k = beta_actual;

        beta_actual_k(randsample(n, kBeta))= r*(2*rand(kBeta,1) - 1);

        delta_actual_k = diag(beta_actual_k);
        phi_actual_k = (A + (B * delta_actual_k));

        y_actual = phi_actual_k * U * theta_actual_k;
        y_measured = y_actual + noise;

                
        %% ESTIMATE

%                beta_estimated = zeros(n,1);
        beta_estimated = r*(2*rand(n,1) - 1);
        theta_estimated = zeros(n,1);
        z_estimated = zeros(n,1);

        cvx_begin
            variables theta_estimated(n) z_estimated(n)
            minimize(norm([theta_estimated; z_estimated], 1))
            subject to
                norm(y_measured - [A*U B]*[theta_estimated; z_estimated], 2) <= epsilon
        cvx_end

        theta_error = norm(theta_estimated - theta_actual_k, 1);

    theta_error_averaged = theta_error_averaged  + theta_error;
    
end % run loop
toc
theta_error_averaged = theta_error_averaged / numruns;    


beep; pause(1); beep; pause(1); beep; 
%save('theta_error_averaged_fullbeta.mat', 'mtmat');

%% EVALUATE RESULTS

figure;
plot(1:n, theta_actual_k, '-r', 1:n, theta_estimated, '-b');
title(['m = ' num2str(m) ' | kTheta = ' num2str(kTheta) ' | kBeta = ' num2str(kBeta)]);
legend('actual', 'recovered');


    
    
figure;
plot(1:n, theta_actual_k, '-b', 1:n, U*theta_actual_k, '-r')
    
    
    
    
    
    