function cleanbw(directory)
%
% This function finds all the .tif files in a directory and removes all but
% the largest objects. It then fills all holes in the images
%
% directory : string
%   a string specifying the relative path of the files to be segmented
%
%
% The code was developed by the Brangwynne laboratory at Princeton University. 
% If using this code (or a modified form) please cite:
%
%  W.Gilpin, S. Uppaluri, C. Brangwynne "Worms under pressure: 
%  bulk mechanical properties of C. elegans are independent of the cuticle"Äù 
%  Biophysical Journal, 2015.


addpath(directory);
pics = dir(fullfile(directory,'*.tif'));
mkdir(['clean' directory]);
N = numel(pics);

for k = 1:N
    pic = imread(pics(k).name);
    pic=im2double(pic);
    
%     se = strel('disk',8);    
%     pic = imdilate(pic,se);
    
    pic = ~imfill(~pic,'holes');
    
    imwrite(pic, fullfile(['clean' directory],pics(k).name),'Compression','None');
end