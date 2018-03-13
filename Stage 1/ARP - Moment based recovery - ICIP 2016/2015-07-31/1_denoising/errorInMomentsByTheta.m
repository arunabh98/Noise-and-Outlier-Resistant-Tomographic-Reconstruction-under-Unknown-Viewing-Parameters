%function [IMleastsquares, IMtheta_i] = errorInMomentsByTheta(theta_i, otherthetas, PMtheta_i, PM0, PMotherthetas)
function error = errorInMomentsByTheta(theta_i, otherthetas, PMtheta_i, PM0, PMotherthetas)

    % Calculate difference in thetas as estimated by just theta_i and as by
    % a least squares solution using all the other thetas
    A1 = [cosd(otherthetas), sind(otherthetas)];
    IMleastsquares= A1 \ PMotherthetas;
    
    A2 = [cosd([0; theta_i]), sind([0; theta_i])];
    IMtheta_i = A2 \ [PM0; PMtheta_i];
    
    error = IMtheta_i - IMleastsquares;
    
end
    
