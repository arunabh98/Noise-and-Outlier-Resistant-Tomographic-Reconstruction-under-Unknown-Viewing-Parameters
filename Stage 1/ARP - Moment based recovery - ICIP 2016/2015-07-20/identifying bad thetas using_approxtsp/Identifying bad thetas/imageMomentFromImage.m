function IMDpq = imageMomentFromImage(I, p, q, center, smax)
    %SIGN?
    sum = 0;
    for i = 1:size(I,1)
        y = ((size(I,1) - i) - center(1))/smax;
        for j = 1:size(I,2)
            x = (j - center(2))/smax;
            sum = sum + x^p * y^q * I(i,j);
        end
    end
    IMDpq = sum;
end




    