function numcorrect = RecoverAngles_TSP(imagefile, name, sigmaNoiseFraction, numkeep, numstarts, Ord, randseed)


    % Reconstruction from projection data from unknown angles
    % In this script we'll implement angle recovery. From there, reconstruction
    % is quite directly possible

    close all;
    %clear; clc;
    tag = 'unif';

    rng(5);

    tic
    %% %%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%
    thetas = (0:179); % Projections to sample from
%     numkeep = 170; % Number of angles to recompute from
%     sigmaNoiseFraction = 0.00;

%     imagefile ='../images/200px-mickey.jpg';

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
%    keep = sort(keep);

%      OR
%    Sample from a skewed distribution
% %     weights = [0.2 0.3 0.12 0.03 0.35]; 
% %     partitionsize = length(thetas)/length(weights);
% %     
% %     keep = [];
% %     for i = 1:length(weights)
% %         start = partitionsize * (i - 1);
% %         numsample = round(numkeep * weights(i));
% %         keep_partition = start + randperm(partitionsize, numsample);
% %         keep = [keep  keep_partition];
% %     end
    
    numkeep = length(keep); %to adjust for rounding off errors not balancing out
%%%

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
    if sigmaNoiseFraction > 0
        
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

    %% %%%%%%%%%%%%%%% C. RECOVER ANGLES FROM PROJECTIONS %%%%%%%%%%%%%%%%%
    
    [tsporder1, tsporder2, d1] = orderbytsp(Pgiven);
    [~, actualorder] = sort(thetaskeep);
    
    [stsporder1, stsporder2, d2] = orderbytsp(Pgiven(:,actualorder));
    
    uniformangles = 0:180/numkeep:179;
    tsp_predicted_angles = uniformangles(tsporder1);
    stsp_predicted_angles = uniformangles(stsporder1);
    
    figure;
    plot(expectedanswer, tsp_predicted_angles, '.r')
    saveas(gcf, ['../output_tsp/' name '_' tag '_numangles' num2str(numkeep) '_noise' num2str(sigmaNoiseFraction*100) 'pc_projections_shuffledstart.jpg'], 'jpg');
    
    figure;
    plot(expectedanswer(actualorder), stsp_predicted_angles, '.b')
    saveas(gcf, ['../output_tsp/' name '_' tag '_numangles' num2str(numkeep) '_noise' num2str(sigmaNoiseFraction*100) 'pc_projections_sortedstart,jpg'], 'jpg');
    
    save(['../output_tsp/' name '_' tag '_numangles' num2str(numkeep) '_noise' num2str(sigmaNoiseFraction*100) 'pc_projections.mat']);

    %% %%%%%%%%%%%%%%% C. RECOVER ANGLES FROM PROJECTION MOMENTS %%%%%%%%%%%%%%%%%
    
    Mmatrix = M';
    [mtsporder1, mtsporder2, md1] = orderbytsp(Mmatrix);
    [~, actualorder] = sort(thetaskeep);
    
    [mstsporder1, mstsporder2, md2] = orderbytsp(Mmatrix(:,actualorder));
    
    uniformangles = 0:180/numkeep:179;
    mtsp_predicted_angles = uniformangles(mtsporder1);
    mstsp_predicted_angles = uniformangles(mstsporder1);
    
    figure;
    plot(expectedanswer, mtsp_predicted_angles, '.r')
    saveas(gcf, ['../output_tsp/' name '_' tag '_numangles' num2str(numkeep) '_noise' num2str(sigmaNoiseFraction*100) 'pc_moments_shuffledstart.jpg'], 'jpg');
    
    figure;
    plot(expectedanswer(actualorder), mstsp_predicted_angles, '.b')
    saveas(gcf, ['../output_tsp/' name '_' tag '_numangles' num2str(numkeep) '_noise' num2str(sigmaNoiseFraction*100) 'pc_moments_sortedstart,jpg'], 'jpg');
    
    save(['../output_tsp/' name '_' tag '_numangles' num2str(numkeep) '_noise' num2str(sigmaNoiseFraction*100) 'pc_moments.mat']);
     
    
    numcorrect=0;
    


end
