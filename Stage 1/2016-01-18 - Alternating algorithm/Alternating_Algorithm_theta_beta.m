clear; clc; 
rng(5);

%% PARAMETERS

n = 200;
r = 1;                 
tolerance = 1e-6; % difference in consecutive estimates of x
epsilon = 1e-2; % for optimization threshold
epsilon_noise = epsilon;
maxiter = 250;

numruns = 1;
m = 120;
kTheta = 10
kBeta = 20;

%% INITIALIZE

R = r * ones(n,1);
theta_error_averaged = 0;

list_theta_error = zeros(maxiter, 1);
list_beta_error = zeros(maxiter, 1);

tic
for run = 1:numruns;

    disp('====================================');
    disp(['              run ' num2str(run) '                 ']);
    
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
        
        beta_estimated = r*(2*rand(n,1) - 1);
        theta_estimated = zeros(n,1);
        z_estimated = zeros(n,1);

        iter = 0;
        rel_diff_theta = Inf;
        
        while rel_diff_theta > tolerance && iter < maxiter
            iter = iter + 1
             
            %Solve for theta
            prev_theta_estimated = theta_estimated;
            
            cvx_begin
                variables theta_estimated(n)
                 minimize(norm(theta_estimated, 1))
                subject to
                    norm(y_measured - ((A + B*diag(beta_estimated)) * U * theta_estimated), 2) <= epsilon
            cvx_end
            
            norm_diff_theta = norm(theta_estimated - prev_theta_estimated, 1);
            rel_diff_theta = abs(norm(theta_estimated,1) - norm(prev_theta_estimated, 1)) / norm(prev_theta_estimated, 1);
            
            list_theta_error(iter) = norm(theta_actual - theta_estimated, 1);
    
            %Solve for beta
            cvx_begin
                variable beta_estimated(n)
                minimize norm(beta_estimated, 1)
                subject to
                    norm(y_measured - ((A + B*diag(beta_estimated)) * U * theta_estimated), 2) <= epsilon
                    -R <= beta_estimated <= R
            cvx_end

%            norm_diff_beta = norm(beta_estimated - prev_beta_estimated, 1);
%            rel_diff_beta = abs(norm(beta_estimated,1) - norm(prev_beta_estimated, 1)) / norm(prev_beta_estimated, 1);
            list_beta_error(iter) = norm(beta_actual - beta_estimated, 1);
    
        end
        
theta_error_averaged = theta_error_averaged  + list_theta_error(iter);
end % run loop

toc
theta_error_averaged = theta_error_averaged / numruns;    


%% EVALUATE RESULTS

figure;
plot(1:n, theta_actual_k, '-r', 1:n, theta_estimated, '-b');
title(['m = ' num2str(m) ' | kTheta = ' num2str(kTheta) ' | kBeta = ' num2str(kBeta)]);
legend('theta actual', 'theta recovered');

figure;
plot(1:n, beta_actual_k, '-r', 1:n, beta_estimated, '-b');
title(['m = ' num2str(m) ' | kTheta = ' num2str(kTheta) ' | kBeta = ' num2str(kBeta)]);
legend('beta actual', 'beta recovered');


    
    
    
    
    
    