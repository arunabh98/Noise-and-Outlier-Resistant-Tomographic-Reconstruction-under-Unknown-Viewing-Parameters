
function I = ibeam()
    % Ibeam
    I = zeros(100,100);
    x1 = 20; y1 = 10;
    x2 = 80; y2 = 30;
    x3 = 40; y3 = 30;
    x4 = 60; y4 = 70;
    x5 = 10; y5 = 70;
    x6 = 90; y6 = 90;

    I(y1:y2, x1:x2) = 1;
    I(y3:y4, x3:x4) = 1;
    I(y5:y6, x5:x6) = 1;
end

