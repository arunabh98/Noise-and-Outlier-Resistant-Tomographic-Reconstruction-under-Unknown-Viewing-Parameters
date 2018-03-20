function [processed_theta] = process_theta(theta)
    
    a = (theta > 180);
    b = (theta < -180);
    c = (-180 < theta < 180);
    A = zeros(size(theta)) + 360;
    B = zeros(size(theta)) - 360;
    C = zeros(size(theta));
    A(a) = theta(a);
    B(b) = theta(b);
    C(c) = theta(c);
    A = A - 360;
    B = B + 360;
    processed_theta = A + B + C;
end
