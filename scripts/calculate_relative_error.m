filename = ...
    '../results/moment_estimation/unknown_angles/5_percent_noise/';
num_angles = 80;

I1 = mat2gray(imread(strcat(filename, num2str(num_angles), '/original_image.png')));
I2 = mat2gray(imread(strcat(filename, num2str(num_angles), '/estimated_image.png')));
I3 = mat2gray(imread(strcat(filename, num2str(num_angles), '/reconstructed.png')));

disp(norm(I1-I2)/norm(I1));
disp(norm(I1-I3)/norm(I1));

