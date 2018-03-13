
clear; clc; 
rng(5);

%% PARAMETERS
n = 200;
r = 0.1;
tolerance = 1e-6; % difference in consecutive estimates of x
epsilon = 1e-2; % for optimization threshold
epsilon_noise = epsilon;
numruns = 5;
mlist = 0:10:150;
klist = 10:5:60;

%% INITIALIZE

U = dctmtx(n);
theta_error_averaged = zeros(numel(mlist), numel(klist));

for run = 1:numruns;

    disp('====================================');
    disp(['              run ' num2str(run) '                 ']);
    theta_error = zeros(numel(mlist), numel(klist));

    for mcounter = 1:numel(mlist)

        m = mlist(mcounter) % number of measurements

        %% GENERATE DATA

        A = normrnd(2, 5, [m, n]);
        B = normrnd(2, 5, [m, n]);

        beta_actual = r*(2*rand(n,1) - 1); 
        delta_actual = diag(beta_actual);
        phi_actual = (A + (B * delta_actual));

        theta_actual = zeros(n,1);

        for kcounter = 1:numel(klist)% k-sparse signal
            k = klist(kcounter)

            theta_actual_k = theta_actual;
            theta_actual_k(randsample(n, k)) = normrnd(0, 1, [k, 1]);

            y_actual = phi_actual * U * theta_actual_k;

            noise = zeros(n, 1);
            noise = normrnd(0, 1, [m,1]);
            noise = noise * epsilon_noise / norm(noise, 2);

            y_measured = y_actual + noise;


            %% ESTIMATE

            beta_estimated = zeros(n,1);
            theta_plus = zeros(n,1);
            theta_minus = zeros(n,1);
            p = zeros(n,1);
            M = [A -A B];
            R = r * ones(n,1);

            cvx_begin
                variables theta_plus(n) theta_minus(n) p(n)
                minimize(sum(theta_plus + theta_minus))

                subject to
                    norm(y_measured - M*[U*theta_plus; U*theta_minus; U*p], 2) <= epsilon

                    theta_plus >= 0
                    theta_minus >= 0
                    -r*(theta_plus + theta_minus) <= p <= r*(theta_plus + theta_minus) 
            cvx_end

            assert(all(theta_plus .* theta_minus < 1e-8), 'CHECK FAIL: x+ and x- overlap')
            theta_estimated = theta_plus - theta_minus;
            x = U * theta_estimated;
            beta_estimated(theta_estimated > 1e-8) = p(theta_estimated > 1e-8) ./ x(theta_estimated > 1e-8); % Automatically avoids div by 0 locations

            theta_error(mcounter, kcounter) = norm(theta_estimated - theta_actual_k, 1);



        end % k loop

    end % m loop

    theta_error_averaged = theta_error_averaged  + theta_error;
    
end % run loop

theta_error_averaged = theta_error_averaged / numruns;    

%% EVALUATE RESULTS

%
% figure;
% plot(1:n, theta_actual_k, '-r', 1:n, theta_estimated, '-b');
% legend('actual', 'recovered');
% norm(theta_actual_k - theta_estimated, 1)
%     

beep; pause(1); beep; pause(1); beep; 
save('theta_error_averaged.mat', 'theta_error_averaged');
save('env_plusminusmethod5.mat');
image(klist, mlist, repmat(theta_error_averaged / max(max(theta_error_averaged)),[1,1,3]))
xlabel('k (number of non-zero entries)');
ylabel('m (number of measurements)');

title('L1 norm of error vector')


