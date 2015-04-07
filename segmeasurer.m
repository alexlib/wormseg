close all; %clear allsegrate1worm3
clearvars pics;
directory='segdevice1worm1';
addpath(directory);
pics=dir(fullfile(directory,'*.tif')); % get all .tif files in the specified MATLAB relative path
N = numel(pics);
data=zeros(N,4); % each column is different measurement
for k = 1:N
    pic = ~imread(pics(k).name);
    comp = bwconncomp(pic);
    if comp.NumObjects==0
        data(k,1) = 0;
        data(k,2) = 0;
        data(k,3) = 0;
        data(k,4) = 0;
    else
        props = regionprops(comp, pic, 'MajorAxisLength','MinorAxisLength','Area','Perimeter');
        areas=cell2mat({props(:).Area});
        [maxsize,windex]=max(areas);
        data(k,1)=sum(areas);
        data(k,2)=props(windex).MajorAxisLength;
        data(k,3)=props(windex).MinorAxisLength;
        data(k,4)=maxsize; % size of just he biggest object
    end
end
plot(data(:,1));
% save(directory,directory)
% plot((data(:,1)-data(1,1))/data(1,1))
% hold all
%plot((data(:,2)-data(1,2))/data(1,2))
% plot(data(:,2)./data(:,1))