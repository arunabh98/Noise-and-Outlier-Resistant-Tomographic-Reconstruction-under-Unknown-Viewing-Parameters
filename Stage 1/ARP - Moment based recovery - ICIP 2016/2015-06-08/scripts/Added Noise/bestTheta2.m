function [bestt, besterror] = bestTheta2(IM, PM, theta1)
    if ~exist('theta1', 'var')
        theta1 = 0;
    end

    %%  Closed form solution
    phi = atand(IM(1) / IM(2));
    theta =  asind( PM(2) / sqrt((IM(1)* IM(1)) + (IM(2)*IM(2))) ) - phi;
    
    %%%TODO: problem with range of asind(), since we're subtacting phi
    %%% later on!
    
    A = [cosd(theta1), sind(theta1); cosd(theta), sind(theta)];
    disp(A);
    IMcalculated = A \ PM;
    error = sum(sum(abs(IMcalculated - IM)));
    
    bestt = theta;
    besterror = error;
end
