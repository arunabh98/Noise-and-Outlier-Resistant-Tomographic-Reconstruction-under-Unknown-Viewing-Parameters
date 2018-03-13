
clear; clc; 
rng(5);

%% PARAMETERS
n = 200;
r = 0.1;
tolerance = 1e-6; % difference in consecutive estimates of x
epsilon = 1e-2; % for optimization threshold
epsilon_noise = epsilon;
numruns = 1;
mlist = 10:20:150;
klist = 10:10:60;
maxiter = 100;

%% INITIALIZE

x_error_averaged = zeros(numel(mlist), numel(klist));

tic
for run = 1:numruns;

    disp('====================================');
    disp(['              run ' num2str(run) '                 ']);
    x_error = zeros(numel(mlist), numel(klist));

    for mcounter = 1:numel(mlist)

        m = mlist(mcounter) % number of measurements

        %% GENERATE DATA

        A = normrnd(2, 5, [m, n]);
        B = normrnd(2, 5, [m, n]);

        beta_actual = r*(2*rand(n,1) - 1); 
        delta_actual = diag(beta_actual);
        phi_actual = (A + (B * delta_actual));

        x_actual = zeros(n,1);

        for kcounter = 1:numel(klist)% k-sparse signal
            k = klist(kcounter)

            x_actual_k = x_actual;
            x_actual_k(randsample(n, k)) = normrnd(0, 1, [k, 1]);

            y_actual = phi_actual * x_actual_k;

            noise = zeros(n, 1);
            noise = normrnd(0, 1, [m,1]);
            noise = noise * epsilon_noise / norm(noise, 2);

            y_measured = y_actual + noise;


            %% ESTIMATE

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
            x_error(mcounter, kcounter) = norm(x_estimated - x_actual_k, 1);


        end % k loop

    end % m loop

    x_error_averaged = x_error_averaged  + x_error;
    
end % run loop
toc
x_error_averaged = x_error_averaged / numruns;    

%% EVALUATE RESULTS

% figure;
% plot(1:n, x_actual_k, '-r', 1:n, x_estimated, '-b');
% legend('actual', 'recovered');
% 
% norm(x_actual_k - x_estimated, 1)    % AA-
% 
%     

beep; pause(1); beep; pause(1); beep; 
save('x_error_averaged_AA.mat', 'x_error_averaged');

image(klist, mlist, repmat(x_error_averaged / max(max(x_error_averaged)),[1,1,3]))
xlabel('k (number of non-zero entries)');
ylabel('m (number of measurements)');
title('L1 norm of error vector')


