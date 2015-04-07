function skelmovie(directory)
%
% This function finds all the segmented .tif files in a directory and
% generates new images with a best-fit skeleton superimposed. Do not use
% this module to actually take measurements.
%
% directory : string
%   a string specifying the relative path of the files to be skeletonized
%
%
% The code was developed by the Brangwynne laboratory at Princeton University. 
% If using this code (or a modified form) please cite:
%
%  W.Gilpin, S. Uppaluri, C. Brangwynne "Worms under pressure: 
%  bulk mechanical properties of C. elegans are independent of the cuticle"€ 
%  Biophysical Journal, 2015.

addpath(directory);
pics = dir(fullfile (directory,'*.tif')); % get .tif files in directory
mkdir(['skel' directory]);                 % make an output directory
N = numel(pics);
for k = 1:N
    pic = imread(pics(k).name);
    pic=im2bw(pic);
%     bw=pic;
    kk=~pic;
    mm=bwmorph(skeleton(kk)>35,'skel',Inf); % tune this parameter
    se = strel('disk',3);                   % tune this parameter
    mm = imdilate(mm,se);
    % mm=~mm;
    RGBmm = imcomplement(cat(3, 0.278431*mm, 0.788235*mm, 0.478431*mm));
    RGBmm(find(RGBmm==0))=1;
    RGBkk=cat(3, ~kk, ~kk, ~kk);
    newpic=RGBkk+imcomplement(RGBmm);
%     newpic = mm; % turn on this line to get just skeleton
    imwrite(newpic, fullfile(['skel' directory],pics(k).name),'Compression','None');
end