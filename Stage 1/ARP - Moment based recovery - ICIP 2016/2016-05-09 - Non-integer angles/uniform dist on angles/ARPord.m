function numcorrect = ARP(imagefile, name, sigmaNoiseFraction, numkeep, numstarts, Ord, randseed)


    % Reconstruction from projection data from unknown angles
    % In this script we'll implement angle recovery. From there, reconstruction
    % is quite directly possible

    %close all; clear; clc;

    rng(randseed);

    tic
    %% %%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%
    thetas = (0:179); % Projections to sample from
    thetas = thetas + rand(1,numel(thetas));
%     numkeep = 101; % Number of angles to recompute from
%     sigmaNoiseFraction = 0.05;

%     numstarts = 10;
%     Ord = 2;

    %% %%%%%%%%%%%%% A. GENERATE INPUT DATA %%%%%%%%%%%%%%%%%

    I= rgb2gray(imread(imagefile));
    %  figure();
    %  imshow(I,'InitialMagnification','fit');

    %%%% Calculate projections %%%%

    [P, svector] = radon(I, thetas); % P(s, theta)

    % Normalize s to a unit circle
    smax = max(abs(svector));
    svector = svector / smax;

    % Sample Uniformly
    keep = randperm(size(P,2), numkeep);

%      OR
%     Sample from a skewed distribution
%     weights = [0.2 0.3 0.12 0.03 0.35]; 
%     partitionsize = length(thetas)/length(weights);
%     
%     keep = [];
%     for i = 1:length(weights)
%         start = partitionsize * (i - 1);
%         numsample = round(numkeep * weights(i));
%         keep_partition = start + randperm(partitionsize, numsample);
%         keep = [keep  keep_partition];
%     end
%     
%     numkeep = length(keep); %to adjust for rounding off errors not balancing out

    thetaskeep = thetas(keep);
    Pgiven = P(:, keep);
    Pnonoise = Pgiven;

    expectedanswer = (thetaskeep - thetaskeep(1))';
    expectedanswer  = expectedanswer * mode(sign(expectedanswer));

    %%%% Add noise
    ref = std(Pgiven(:));
    sigmaNoise = sigmaNoiseFraction * ref;
    noise = normrnd(0, sigmaNoise, size(Pgiven));
    Pgiven = Pgiven + noise;


    %%%% Calculate Projection moments from projections %%%%

    kmax = numkeep - 1; %Number of projection moments. This equality is because this kmax is the number needed for calculating image moments
    %"The image moments of the kth order are determined by 
    % a set of projection moments of the kth order from
    % k+1 given views."

    M = ones(size(Pgiven,2), kmax+1); % each column represents kth moments for one k, for all directions 
    for k = 0:kmax
        for i = 1:size(Pgiven,2)
            M(i, k+1) = calculateProjectionMoment(Pgiven(:,i), svector, k);
        end
    end

    %The resulting M is the limited set of Projection Moments
    % M(thetanumber, momentnumber)
          M(:,1); % 0th projection moment
          M(:,2); %

    %%%% Calculate Real image moments
    IMbasetruth = [ imageMomentFromImage(I, 1, 0, size(I)/2, smax);
                    imageMomentFromImage(I, 0, 1, size(I)/2, smax); 
                  ];


    %% %%%%%%%%%%% B. DENOISING %%%%%%%%
    if sigmaNoiseFraction > -1 % Temporarily replacing 0. As such noise fraction will always be non-negative.
        
        disp('Before denoising');
        Pnoisy = Pgiven;
        norm(Pnonoise - Pgiven, 'fro')
        psnr_noisy = 10 * log10( max(max(Pnonoise))^2 / (mean(mean(( Pnonoise - Pnoisy ) .^ 2))) )

        P_denoised = denoise(Pgiven, sigmaNoise, 50, 200);
        Pgiven = P_denoised;

        disp('After denoising')
        norm(Pnonoise - P_denoised, 'fro')
        %norm(Pnonoise - P_denoised2, 'fro')
%         figure()
%         plot(Pnonoise(:,1));
%         hold on;
%         plot(P_denoised(:,1), '-r');
%         plot(Pnoisy(:,1), '-g');
%         hold off;

        psnr_denoised = 10 * log10( max(max(Pnonoise))^2 / (mean(mean(( Pnonoise - P_denoised ) .^ 2))) )

    end
    %% %%%%%% B2. Final fixes on P, M %%%%%%%%%%%%%%%
    Pgiven = max(0, Pgiven);

    % Recalculate M
    M=ones(size(Pgiven,2), kmax+1); % each column represents kth moments for one k, for all directions
    for k = 0:kmax
        for i = 1:size(Pgiven,2)
            M(i, k+1) = calculateProjectionMoment(Pgiven(:,i), svector, k);
        end
    end

    %% %%%%%%%%%%%%%%% C. RECOVER ANGLES %%%%%%%%%%%%%%%%%

    % Initialize angles randomly
    thetasestimated = randi([1 179], numkeep, 1);
    %thetasestimated(1) = 0;
    thetasfirstestimate = thetasestimated;

    PMord = reshape(M(:,2:(Ord+1)), Ord * numkeep, 1);

    numiter = 30;
    Ehlccvalues = zeros(numiter * (numkeep - 1), 1);
    
    deltas = zeros(numiter, 1);

    %% Multi-start

    thetasestimated_bystart = zeros(numkeep, numstarts);
    Ehlccvalues_bystart = zeros(numstarts, 1);

    for start = 1:numstarts

        % Initialize angles randomly
        thetasestimated = randi([1 179], numkeep, 1);

        % Starting estimate
        A = assembleA(thetasestimated, Ord); %<>%
        IMestimated = A \ PMord;

        Ehlccvec = A * IMestimated - PMord; 
        Ehlcc = norm(Ehlccvec, 2);
        ctr = 1;

        Ehlccvalues(ctr) = Ehlcc;


        for iteration = 1:numiter
            disp([num2str(start), ', ', num2str(iteration)]);

            for i = 1:numkeep
                ctr = ctr + 1;

                besttheta = thetasestimated(i);
                bestE = Ehlcc;

                for t = [-179:180]

                    thetas_iter = thetasestimated;
                    thetas_iter(i) = t;
                    A = assembleA(thetas_iter, Ord);
                    IMestimated = A \ PMord;
                    E_tvec = A * IMestimated - PMord;
                    E_t = norm(E_tvec, 2);

                    if E_t < bestE
                        besttheta = t;
                        bestE = E_t;
                    end
                end

                thetasestimated(i) = besttheta;


                A = assembleA(thetasestimated, Ord);
                IMestimated = A \ PMord;

                Ehlccvec = A * IMestimated - PMord;
                Ehlcc = norm(Ehlccvec, 2);
                Ehlccvalues(ctr) = Ehlcc;

            end

            ctr
            Ehlcciter = Ehlccvalues(ctr)
            Ehlccpreviter = Ehlccvalues(ctr - (numkeep)); % -1 iff theta0 is not updated

            delta = Ehlccpreviter - Ehlcciter;
            deltas(iteration) = delta;
            if delta < 0.0001
                break
            end

        end
        Ehlccvalues_bystart(start) = Ehlccvalues(ctr);
        thetasestimated_bystart(:, start) = thetasestimated;
    end % multistart loop


    %figure();
    %plot(Ehlccvalues_bystart, 'or');

    [~, optstart] = min(Ehlccvalues_bystart);
    thetasestimated = thetasestimated_bystart(:, optstart);

    %thetasestimated = mod(thetasestimated, 360);
    [~, sortingorder] = sort(expectedanswer);
    disp('Expected | Converged | (1) -(2) | (1) + (2)');
    %[expectedanswer(sortingorder) thetasestimated(sortingorder)
    %(mod(expectedanswer(sortingorder) - thetasestimated(sortingorder), 360)) (mod((expectedanswer(sortingorder) + thetasestimated(sortingorder)), 360))]
    [expectedanswer thetasestimated (mod(expectedanswer - thetasestimated, 360)) (mod((expectedanswer + thetasestimated), 360))]
    toc


    %Checker
    a = mean((mod(expectedanswer - thetasestimated, 360)))
    b = mean((mod(expectedanswer + thetasestimated, 360)))


    %err_a = sum(abs((mod(expectedanswer - thetasestimated, 360)) - a) <= 0);
    %err_b = sum(abs((mod(expectedanswer + thetasestimated, 360)) - b) <= 0);
    err_a = std(mod(expectedanswer - thetasestimated, 360));
    err_b = std(mod(expectedanswer + thetasestimated, 360));
    


    if err_a < err_b
        e0 = sum(abs((mod(expectedanswer - thetasestimated, 360)) - a) <= 0);
        ehalf = sum(abs((mod(expectedanswer - thetasestimated, 360)) - a) <= 0.5);
        e1 = sum(abs((mod(expectedanswer - thetasestimated, 360)) - a) <= 1);
        e3 = sum(abs((mod(expectedanswer - thetasestimated, 360)) - a) <= 3);
        e5 = sum(abs((mod(expectedanswer - thetasestimated, 360)) - a) <= 5);

    else
        e0 = sum(abs((mod(expectedanswer + thetasestimated, 360)) - b) <= 0);
        ehalf = sum(abs((mod(expectedanswer + thetasestimated, 360)) - b) <= 0.5);
        e1 = sum(abs((mod(expectedanswer + thetasestimated, 360)) - b) <= 1);
        e3 = sum(abs((mod(expectedanswer + thetasestimated, 360)) - b) <= 3);
        e5 = sum(abs((mod(expectedanswer + thetasestimated, 360)) - b) <= 5);
    end
        
    e0
    ehalf
    e1
    e3
    e5
    numcorrect = char([num2str(e0), '/', num2str(ehalf), '/', num2str(e1), '/', num2str(e3), '/', num2str(e5)]);

% %     figure();
% %     plot(Ehlccvalues(1:ctr), '-b');
% %     hold on;
% %     majors = 1:(numkeep - 1):ctr;
% %     plot(majors,Ehlccvalues(majors) , 'ro')
% %         %Plot deltas
% %         %figure();
% %         %plot(deltas(1:iteration), '-r.');
% % 
% %     figure()
% %     plot(thetasestimated, expectedanswer, '.r')


    %% %%%%%%%%%%%%%%% D. IDENTIFYING BAD THETAS %%%%%%%%%%%%%%%%%

    
    %% %%%%%%%%%%%%%%%%% IMAGE RECONSTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%

    %Using ORIGINAL angles
    I_reconstructed_original = iradon(Pgiven, thetaskeep,size(Pgiven,1));
    I_reconstructed_original = I_reconstructed_original/max(max(I_reconstructed_original));
    imwrite(I_reconstructed_original, strcat('../output/', name, '_', num2str(sigmaNoiseFraction), '_', num2str(numkeep), '_', num2str(numstarts), '_', num2str(Ord), '_', num2str(randseed), '_actualangles.jpg'));
    %figure();
    %imshow(I_reconstructed_original, 'InitialMagnification','fit');


    %Using ALL RECONSTRUCTED angles (including incorrect ones)
    I_reconstructed = iradon(Pgiven, thetasestimated,size(Pgiven,1));
    I_reconstructed = I_reconstructed/max(max(I_reconstructed));
    imwrite(I_reconstructed, strcat('../output/', name, '_', num2str(sigmaNoiseFraction), '_', num2str(numkeep), '_', num2str(numstarts), '_', num2str(Ord), '_', num2str(randseed), '_estimatedangles.jpg'));
    %figure();
    %imshow(I_reconstructed, 'InitialMagnification','fit');

    % a(logical([1,0,1,1,0,0,0,1,0,0])) = []
    beep; pause(1); beep;

end
