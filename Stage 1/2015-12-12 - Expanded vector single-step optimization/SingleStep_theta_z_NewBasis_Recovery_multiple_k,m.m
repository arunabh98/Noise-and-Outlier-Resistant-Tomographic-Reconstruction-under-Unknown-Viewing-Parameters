clear; clc; 
rng(5);

%% PARAMETERS

n = 200;
r = 0.1;                 
tolerance = 1e-6; % difference in consecutive estimates of x
epsilon = 1e-2; % for optimization threshold
epsilon_noise = epsilon;

numruns = 5;
mlist = 10:20:150;
kThetalist = 10:10:60;
kBetalist = 15:15:90;

%% INITIALIZE

R = r * ones(n,1);
theta_error_averaged = zeros(numel(mlist), numel(kThetalist), numel(kBetalist));

tic
for run = 1:numruns;

    disp('====================================');
    disp(['              run ' num2str(run) '                 ']);
    theta_error = zeros(numel(mlist), numel(kThetalist), numel(kBetalist));

    for mcounter = 1:numel(mlist)

        m = mlist(mcounter) % number of measurements

        %% GENERATE DATA

        A = normrnd(2, 2, [m, n]);
        B = normrnd(2, 2, [m, n]);
        U = dctmtx(n);

        theta_actual = zeros(n,1);
        beta_actual = zeros(n,1);
        
        noise = normrnd(0,1,[m,1]);
        noise = noise * epsilon_noise / norm(noise, 2);

        
        for kThetacounter = 1:numel(kThetalist)% signal sparsity
            kTheta = kThetalist(kThetacounter)

            theta_actual_k = theta_actual;
            theta_actual_k(randsample(n, kTheta)) = normrnd(0, 1, [kTheta, 1]);

            for kBetacounter = 1:numel(kBetalist)% k-sparse signal  
                [run mcounter kThetacounter kBetacounter]
                
                kBeta = kBetalist(kBetacounter)
                beta_actual_k = beta_actual;
            
                beta_actual_k(randsample(n, kBeta))= r*(2*rand(kBeta,1) - 1);
                
                delta_actual_k = diag(beta_actual_k);
                phi_actual_k = (A + (B * delta_actual_k));

                y_actual = phi_actual_k * U * theta_actual_k;
                y_measured = y_actual + noise;

                
                %% ESTIMATE

                beta_estimated = zeros(n,1);
                beta_estimated = r*(2*rand(n,1) - 1);
                theta_estimated = zeros(n,1);
                z_estimated = zeros(n,1);

                cvx_begin
                    variables theta_estimated(n) z_estimated(n)
                    minimize(norm([theta_estimated; z_estimated], 1))
                    subject to
                        norm(y_measured - [A*U B]*[theta_estimated; z_estimated], 2) <= epsilon
                cvx_end

                theta_error(mcounter, kThetacounter, kBetacounter) = norm(theta_estimated - theta_actual_k, 1);
            end % kBeta loop

        end % kTheta loop

    end % m loop

    theta_error_averaged = theta_error_averaged  + theta_error;
    
end % run loop
toc
theta_error_averaged = theta_error_averaged / numruns;    


beep; pause(1); beep; pause(1); beep; 
save('theta_error_averaged.mat', 'theta_error_averaged');

%% EVALUATE RESULTS

% figure;
% plot(1:n, x_actual_k, '-r', 1:n, x_estimated, '-b');
% legend('actual', 'recovered');
% 
% norm(x_actual_k - x_estimated, 1)    % AA-
% 
%     


f1 = figure;
for i = 1:numel(mlist) 
    subplot(2,4,i)
    btmat = squeeze(theta_error_averaged(i,:,:))
    image(kThetalist, kBetalist, repmat(btmat / max(max(btmat)),[1,1,3]))
    ylabel('kTheta');
    xlabel('kBeta');
    title(['m = ' num2str(mlist(i))]);
end
suptitle('L1 norm of error vector. Normalized separately for each individual plot');
f2 = figure;
for i = 1:numel(mlist) 
    subplot(2,4,i)
    btmat = squeeze(theta_error_averaged(i,:,:))
    image(kThetalist, kBetalist, repmat(btmat / max(max(max(theta_error_averaged))),[1,1,3]));
    ylabel('kTheta');
    xlabel('kBeta');
    title(['m = ' num2str(mlist(i))]);
end
suptitle('L1 norm of error vector. Normalized over all m');    

%%%%%%%%%%%%%%%

f3 = figure;
for i = 1:numel(kThetalist) 
    subplot(2,3,i)
    mbmat = squeeze(theta_error_averaged(:,i,:))
    image(mlist, kBetalist, repmat(mbmat / max(max(mbmat)),[1,1,3]))
    ylabel('m');
    xlabel('kBeta');
    title(['kTheta = ' num2str(kThetalist(i))]);
end
suptitle('L1 norm of error vector. Normalized separately for each individual plot');

f4 = figure;
for i = 1:numel(kThetalist) 
    subplot(2,3,i)
    mbmat = squeeze(theta_error_averaged(:,i,:))
    max(max(mbmat))
    image(mlist, kBetalist, repmat(mbmat /  max(max(max(theta_error_averaged))),[1,1,3]))
    ylabel('m');
    xlabel('kBeta');
    title(['kTheta = ' num2str(kThetalist(i))]);
end
suptitle('L1 norm of error vector. Normalized over all kThetas');


%%%%%%%%%%%%%%%%%%

f5 = figure;
for i = 1:numel(kBetalist) 
    subplot(2,3,i)
    mtmat = squeeze(theta_error_averaged(:,:,i))
    image(mlist, kThetalist, repmat(mtmat / max(max(mtmat)),[1,1,3]))
    ylabel('m');
    xlabel('kTheta');
    title(['kBeta = ' num2str(kBetalist(i))]);
end
suptitle('L1 norm of error vector. Normalized separately for each individual plot');

f6 = figure;
for i = 1:numel(kBetalist) 
    subplot(2,3,i)
    mtmat = squeeze(theta_error_averaged(:,:,i))
    image(mlist, kThetalist, repmat(mtmat / max(max(max(theta_error_averaged))),[1,1,3]))
    ylabel('m');
    xlabel('kTheta');
    title(['kBeta = ' num2str(kBetalist(i))]);
end
suptitle('L1 norm of error vector. Normalized over all kBetas');

%%%%%%%%%

figure;
imshow3D(theta_error_averaged) %scroll on kBeta
    xlabel('m');
    ylabel('kTheta');

figure;
imshow3D(permute(theta_error_averaged, [1 3 2])) %scroll on kTheta
    xlabel('m');
    ylabel('kBeta');

figure;
imshow3D(permute(theta_error_averaged, [2 3 1])) %scroll on m
    xlabel('kTheta');
    ylabel('kBeta');















