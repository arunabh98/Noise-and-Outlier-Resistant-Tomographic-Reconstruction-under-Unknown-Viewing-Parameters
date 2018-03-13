
close all;
clear; clc;
%% PARAMETERS
%imagelist = {'200px-cat.jpg', '200px-mickey.jpg', '200px-morningsun.jpg', '200px-rocket.jpg', '200px-steve.jpg', '200px-beach.jpg', '200px-woody.jpg', '256px-Earthrise.jpg'};
imagelist = {'256px-Earthrise.jpg', '200px-mickey.jpg'};
%noiselist = {0, 0.01, 0.05, 0.1};% all the values of noise we want to loop over
noiselist = {0};% all the values of noise we want to loop over
numkeep = 30; numstarts = 10; randseed = 5;
ordlist = [6];


%% %%%
imagelist'
fh = fopen('outfile.txt', 'a');
fprintf(fh, '%2s\n', '********************'); 
ctr = 0;
for Ord = ordlist
    for im = imagelist
        for noise = noiselist
            ctr = ctr + 1
            name = strsplit(char(im), '[\.\-]','DelimiterType','RegularExpression');
            name = char(name(2));

            imfile = char(['../images/', char(im)]);
            noisefrac = cell2mat(noise);

            config = sprintf('%2s,%2s,%1.2f,%2d,%2d,%2d,%2d\n', char(imfile), name, noisefrac,  numkeep, numstarts, Ord, randseed)
            numcorrect = ARPord(imfile, name, noisefrac, numkeep, numstarts, Ord, randseed)
            fprintf(fh, '%2s,%2s,%1.2f,%2d,%2d,%2d,%2d,%2s\n', char(imfile), name, noisefrac,  numkeep, numstarts, Ord, randseed, numcorrect);
            disp('----------------------------------------------')
            

        end
    end
end


fclose(fh);
beep; pause(1); beep; pause(1); beep; pause(1); beep; pause(1); beep;