function finaldata = seglength(directory)
% enter a string representing a directory to probe
% >> data = finaldata('segdevice3rate5worm1');
close all;
clearvars pics;
addpath(directory);
pics=dir(fullfile(directory,'*.tif')); % get all .tif files in the specified MATLAB relative path
N = numel(pics);
data=zeros(N,4);                % each column is different measurement
for k = 1:N
    pic = ~imread(pics(k).name);
    pic2=pic;
    pic2=bwmorph(skeleton(pic2)>35,'skel',Inf);
    [L,num] = bwlabel(pic);
    if num==0
        data(k,:)=0;
        disp('no objects found in frame ');
        disp(k);
    else
        areas=zeros(num,1);
        for ii=1:num
            areas(ii) = bwarea(L==ii);
        end
        data(k,1)=sum(areas);
        data(k,2)=max([areas; 0]);
        [L,num] = bwlabel(pic2);
        areas=zeros(num,1);
        for ii=1:num
            areas(ii) = bwarea(L==ii);
        end
        data(k,3)=sum(areas);
        data(k,4)=max([areas; 0]);
    end
end

finaldata=data;


%     comp = bwconncomp(pic);
%     comp2=bwconncomp(pic2);
%     
%     if comp.NumObjects==0
%         data(k,1) = 0;
%         data(k,2) = 0;
%         data(k,3) = 0;
%         data(k,4) = 0;
%     else



    %             props = regionprops(comp, pic,'Area');
    %             areas=cell2mat({props(:).Area});
    %             [maxsize,windex]=max(areas);
    %             data(k,1)=sum(areas); % sum of just the biggest object
    %             data(k,2)=maxsize; % size of just the biggest object
    %             props = regionprops(comp2, pic2,'Area');
    %             areas=cell2mat({props(:).Area});
    %             [maxsize,windex]=max(areas);
    %             data(k,3)=maxsize;


