% directory = '/Users/williamgilpin/Desktop/SoftLiving/Time lapse images/12_7 N2 small';
directory='test_images';
addpath(directory);
pics = dir(fullfile(directory,'*.tif')); % 10/29/2014 Edited this line to include a space after fullfile''
% get all .tif files in the specified absolute path
mkdir(['seg' directory]);
N = numel(pics);
for k = 1:N
    pic = imread(pics(k).name);
    pic=im2double(pic);
    pic=pic(:,:,1);
    picref=pic;
    pic=imadjust(pic);
    thres = graythresh(pic);
    %     thres=thres-.2;
    thres=.8;
    pic=im2bw(pic,thres); % use ~ if you want white on black
    imwrite(pic, fullfile(['seg' directory],pics(k).name),'Compression','None');
end