% wormseg
%
% a set of analysis routines for segmentation and skeletonization of 
% large numbers of bright field microscope images
%
% The code was developed by the Brangwynne laboratory at Princeton University. 
% If using this code (or a modified form) please cite:
%
%  W.Gilpin, S. Uppaluri, C. Brangwynne "Worms under pressure: 
%  bulk mechanical properties of C. elegans are independent of the cuticle"Äù 
%  Biophysical Journal, 2015.
%
%

%% Segment all the images in a specified directory
thres = .8; % the brightness threshold for segmentation
dir_str = 'sample_raw_images'; % the directory in which the images are found
segmovie(dir_str, thres)

%% Clean up all the images 
cleanbw('segsample_raw_images');

%% Find skeleton of all images and save images with skeletons overlaid
% this code is optional for measurements because the measurement code calls
% the skeletonization separately. This cell is purely for visualization.
skelmovie('cleansegsample_raw_images')

%% define structure and get rid of anything that isn't a folder

% get all files at current directory level starting with wildcard
folds=dir('clean*');

N = numel(folds);
ii=1;
while ii<=N      
    if ~folds(ii).isdir
        folds(ii)=[];
        ii=ii-1;
    end
    ii=ii+1;
end
disp(N)



%% Find the lengths of everything in the directory
N=numel(folds);
for ii=1:N
    folds(ii).data=seglength(folds(ii).name);
    disp(folds(ii).name);
end
disp('loop complete')



%% OPTIONAL: median filter the measurements to discard any outliers
% this cell generally isn't necessary if you have clean microscopy data.
windex=4;
N=numel(folds);
for ii=9:9
    data=folds(ii).data;
    L=length(data);
    data=medfilt1(data,10);
    data=data(2:end-1,:);
    folds(ii).data=data;
end


%% plot results
close all
aindex=2;
lindex=4;
N=numel(folds);
data=folds(1).data;
figure()
plot(linspace(0,1,length(data)),(pi/4)*(data(:,2).^2./data(:,4))/(1.1039e5)^1.5);
hold all;
for ii=1:N
    data=folds(ii).data;
    plot(linspace(0,1,length(data)),(pi/4)*(data(:,2).^2./data(:,4))/(1.1039e5)^1.5)
end 


%% plot results by volume versus pressure
close all
donormalize=0;
half_periods=3;
aindex=2;
lindex=4;
minvol=.5; % minimum fraction of volume in syringe at max pressure
N=numel(folds);
data=folds(1).data;
figure()
for ii=1:N
    data=folds(ii).data;
    sel=1:floor((1/half_periods)*length(data));
     ps=1./linspace(1,minvol,length(sel))-1;
    if donormalize
        cor_f=(pi/4)*((data(1,aindex).^2)./data(1,lindex)/(1.1039e5)^1.5);
    else
        cor_f=1;
    end
    plot(101.325*ps,(pi/4)*((data(sel,aindex).^2)./data(sel,lindex)/(1.1039e5)^1.5)/cor_f,'.'); hold all
end
xlim([0,101.325]);
% set bottom y limit to zero but keep automatically set top limit
yL = ylim(gca);
ylim([0, yL(2)]);




%% extract measurements from data structure
allcolors=[0.59375 0.257813 0.886719;0.800781 0.0625 0.460938;0.278431 0.788235 0.478431;0.917647 0.682353 0.105882;0.372549 0.596078 1.; 0.8 0.8 .8;1.0 .3882 .2784];
repcolors=[allcolors; allcolors;allcolors;allcolors;allcolors];
N=numel(folds);
lindex=4;
aindex=2;
rates = [.4167 .25 .125 .08333 .0625 .41666 .81888]; % double check these last two numbers (mL/min)
rates = (rates/1.25)*101.325/60; % convert into kPa/s
frac_elas=1/3;
half_periods=3; % select how many times it went down then up
selrange=frac_elas/half_periods; % select fraction of initial data that you think is relevant
minvol=.5; % minimum fraction of volume in syringe at max pressure
maxpress=1;
slopes=zeros(N,1);
slopes2=zeros(N,1);
slopes3=zeros(N,1);
rate=zeros(N,1);
inits=zeros(N,1);
fins=zeros(N,1);
mids=zeros(N,1);
figure();
for ii=1:N
    % initial size
    data=folds(ii).data;
    disp(folds(ii).name)
    lengths=data(:,lindex);
    areas=data(:,aindex);
    inits(ii)=((areas(1)^2)/lengths(1))/(1.1039e5)^1.5;
    fins(ii)=lengths(end);
    mids(ii)=lengths(floor(length(data)/half_periods));
    % pressurization rate
    nam=folds(ii).name;
    locr=findstr(nam, 'rate');
    if ~isempty(locr)                       % only fill if rate is in the name
        rat_ind=str2num(nam(locr+4));
        rate(ii)=rates(rat_ind);
    end
    
    % fit an elastic modulus
    sel = 1:floor(selrange*length(data)); % select some of data to fit
    areas_short=areas(sel);
    lengths_short=lengths(sel);
    ps=(frac_elas*1)./linspace(1,minvol,length(areas_short))-(frac_elas*1);
    ps=1./linspace(1,minvol,floor(length(areas)/half_periods))-1;
    ps_short=ps(sel);
    %ps = linspace(0,frac_elas,length(sel)); % a linear pressure range 
    pts_to_fit=(1./(lengths_short/lengths_short(1))).*( areas_short/areas_short(1) ).^2;
    pts_to_fit=1-pts_to_fit;
    
    % model worm as cylinder
    pts_to_fit(isinf(pts_to_fit))=NaN;
    pts_to_fit(isnan(pts_to_fit))=nanmean(pts_to_fit);
    f0 = fit(pts_to_fit,ps_short','poly1'); % fit in opposite order to get true modulus
    plot(f0,'k',pts_to_fit,ps_short','.');  hold all;
    set(findobj(gca, 'Type', 'Line', 'Color', 'k'), 'Color',repcolors(ii,:),'LineWidth',4,'MarkerSize',60);
    ci = confint(f0,0.95);
    disp(ci)
    slopes(ii) = (ci(2,1)+ci(1,1))/2;
    folds(ii).bulkmod = slopes(ii)*101.325;
    %if half_periods~=1
    if false
        %sel=floor(length(data)/half_periods):floor(length(data)/half_periods+selrange*length(data));
        sel=floor(2*length(data)/half_periods-selrange*length(data)):floor(2*length(data)/half_periods);
        areas_short=areas(sel);
        lengths_short=lengths(sel);
       
        %ps = linspace(0,frac_elas,length(sel)); % a linear pressure range
        pts_to_fit=(1./(lengths_short/lengths_short(1))).*( areas_short/areas_short(1) ).^2;
        pts_to_fit=1-pts_to_fit;
       
        % model worm as cylinder
        pts_to_fit(isinf(pts_to_fit))=NaN;
        pts_to_fit(isnan(pts_to_fit))=nanmean(pts_to_fit);
        dlen=length(sel)-length(ps_short);
        if dlen>0
            pts_to_fit(end-dlen+1:end)=[];
        elseif dlen<0
            ps_short(end-dlen+1:end)=[];
        end
        f0 = fit(pts_to_fit,ps_short','poly1'); % fit in opposite order to get true modulus
        ci = confint(f0,0.95);
        slopes2(ii) = (ci(2,1)+ci(1,1))/2;
        sel=floor(2*length(data)/half_periods):floor(2*length(data)/half_periods+selrange*length(data));
        areas_short=areas(sel);
        lengths_short=lengths(sel);
        pts_to_fit=(1./(lengths_short/lengths_short(1))).*( areas_short/areas_short(1) ).^2;
        pts_to_fit=1-pts_to_fit;
        
        % model worm as cylinder
        pts_to_fit(isinf(pts_to_fit))=NaN;
        pts_to_fit(isnan(pts_to_fit))=nanmean(pts_to_fit);
        dlen=length(sel)-length(ps_short);
        if dlen>0
            pts_to_fit(end-dlen+1:end)=[];
        elseif dlen<0
            ps_short(end-dlen+1:end)=[];
        end
        f0 = fit(pts_to_fit,ps_short','poly1'); % fit in opposite order to get true modulus   
        ci = confint(f0,0.95);
        slopes3(ii) = (ci(2,1)+ci(1,1))/2;
    end
end
% if no rate info, the experiment was done over 15 min
rate(rate==0)=.1126
legend off;legend off;
ylim([min(ps_short) max(ps_short) ])
% xlim([-.03 .25 ])        % have to set this manually
slopes=slopes*101.325; % convert to kPa
slopes2=-slopes2*101.325;
slopes3=slopes3*101.325;
% inits=inits/sqrt(1.1039e5); % convert to mm
figure();
plot(inits(1),slopes(1),'.','MarkerSize',30)
hold all
for ii=2:N
    plot(inits(ii),slopes(ii),'.','MarkerSize',30)
end 
disp(mean(slopes));disp(std(slopes));
disp(mean(slopes2));disp(std(slopes2));
disp(mean(slopes3));disp(std(slopes3));
legend off



%% write the slopes to a text file
fileID = fopen('slopes.txt','w');
fprintf(fileID,'%4.4f\n',slopes);
fclose(fileID);


%% plot hysteresis
allslopes=abs([slopes slopes2 slopes3]);
slopestats=[mean(allslopes); std(allslopes)];
errorbar([1:3],mm(1,:),mm(2,:),'.');hold on;plot(mm(1,:),'.k','MarkerSize',30);


%% look at isotropy
lindex=4;
aindex=2;
N=numel(folds);
figure()
for ii=1:N
    data=folds(ii).data;
    data=data(1:floor(1/3*length(data)),:);
    plot(101.325*linspace(0,1,length(data)),(data(:,lindex)/data(1,lindex))./((data(:,aindex)./data(:,lindex))/(data(1,aindex)/data(1,lindex))),'.'); hold all;
end
xlim([0,101.325]);
ylim([.82 1.18])
pbaspect([3 1 1])


%% look at rates
figure()
for ii=1:N
    plot(rate(ii),slopes(ii),'.k','MarkerSize',30); hold all;
end
xlabel('');
ylabel('');
[f0,gof] = fit(rate(1:end-1),slopes(1:end-1),'poly1');
p11 = predint(f0,rate(1:end-1),0.95,'functional','off');
plot(f0,'r');hold on;plot(rate(1:end-1),p11,'m--');
xlim([0, .7])






%% plot results
% figure();plot(inits,slopes,'.')
%figure();plot(rate,slopes,'.')
 scatter3(inits,slopes,rate,100,rate,'.')
 colormap(jet);
 zlabel('rate');
xlabel('initial size');
ylabel('inverse modulus')
%plot([min1; min2; max1; max2])


%% plot divergence of modulus
% dependency: sgfilter.m on FileExchange, by YangQuan Chen
close all
donormalize=1;
half_periods=3;
aindex=2;
lindex=4;
minvol=.5; % minimum fraction of volume in syringe at max pressure
N=numel(folds);
for ii=1:N
    data=folds(ii).data;
    sel=1:floor((1/half_periods)*length(data)); % take only first pressure stroke
    ps=100./linspace(1,minvol,length(sel))-100; % volume goes to .5 in sel
    vols=(pi/4)*((data(sel,aindex).^2)./data(sel,lindex)/(1.1039e5)^1.5);
    vols=medfilt1(vols,10);
    %vols=data(sel,aindex);
    %bulks=diff(vols)./vols(1:end-1);
    % gauss_ker = gausswin(50)
    % gauss_ker = gauss_ker/sum(gauss_ker); % Normalize.
    % bulks=conv(gauss_ker,bulks);

    % ll=10;
    % b=ones(ll,1)/ll;
    % bulks2=filter(1,b,bulks);
    lh=floor(.15*length(vols)); rh=lh;
    kernel=sgfilter(lh,rh,1,1); % take derivative and fit a line
    dv=conv(kernel,vols);    dv(1:lh)=[];    dv(end-rh+1:end)=[];
    dp=diff(ps)';    dp(end+1)=dp(end);
    bulks=(dp./dv).*vols;
%     bulks=bulks.^(3/2); % convert to true bulk modulus
    [~,win]=min(bulks);
    win=win(1); % only care about first index
    vfrac=1./vols;
    [sorvals,~] = sort(vfrac,'descend');
    sorvals=sorvals(1:10);  % guess asymptotic volume fraction by finding median of final few points
    maxvfrac=median(sorvals); 
    vfrac=vfrac/maxvfrac;       % normalize to final volume fraction
    %vfrac=vols/vols(2);    %plot versus volume, not volume fraction
    vfrac(1:win)=[]; bulks(1:win)=[];
    
    crr=0;
    cll=0;
    loglog(vfrac(lh:end-rh),bulks(lh:end-rh),'.'); hold all
    disp(['data set: ' num2str(ii)]);
    disp([num2str(max(size(vols))) '->' num2str(max(size(bulks(lh:end-rh))))]);
    disp(lh)
    %f=fit(log(vfrac(lh+cll:end-rh-crr)),log(bulks(lh+cll:end-rh-crr)),'poly1')
end


