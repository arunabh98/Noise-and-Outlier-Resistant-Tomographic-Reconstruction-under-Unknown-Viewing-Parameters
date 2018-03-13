function [bestt, besterror] = bestTheta2(IM, PM, theta1)
    if ~exist('theta1', 'var')
        theta1 = 0;
    end

    %%  Closed form solution
    
    phi = atand(IM(1) / IM(2));
    arg = PM(2) / sqrt((IM(1)* IM(1)) + (IM(2)*IM(2)));
    
    arg = max(-1, min(arg, 1));
    theta =  asind( arg ) - phi;
    
    %%%TODO: if arg of asind is not in domain, then cutting off at -1,1 is
    %%% not giving the best least squares solution
    
    A = [cosd(theta1), sind(theta1); cosd(theta), sind(theta)];
    %disp(A);
    IMcalculated = A \ PM;
    error = sum(sum(abs(IMcalculated - IM)));
    
    bestt = theta;
    besterror = error;
end
