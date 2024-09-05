clear all;
% stats of the tracks in the selected region only

% Define image dimensions and load oxygen image
SX=25091; SY=18749;

oxy= (abs(imread('oxy_stack_11_may_small_20_gauss0006.tif')));
oxy_gray=im2gray(oxy);

% smoothen the oxygen level
oxy_gauss=imgaussfilt(oxy_gray,20);

% Make the mask
mask0=oxy_gauss>0.035;
xx=1:size(mask0,2);
yy=1:size(mask0,1);
[XX,YY]=meshgrid(xx,yy);
mask0(XX>950)=0;
mask_large = imresize(mask0,[SY nan]);
oxy_large = imresize(oxy_gauss,[SY nan]);
[gmx,gmy]=gradient(mask_large);
gm=sqrt(gmx.^2+gmy.^2);
cx=[]; cy=[];
for i=1:size(gm,1)
    for j=1:size(gm,2)
        if gm(i,j)>0.001
            cx=[cx,i];
            cy=[cy,j];
        end
    end
end

% Load all cell tracks
load('cellStats.mat');

% Find the tracks within the gradient mask
idM=[]; % id of the tracks within the mask
for i=1:length(xT)
    yP=min(size(mask_large,1),max(1,round(yT(i))));
    xP=min(size(mask_large,2),max(1,round(xT(i))));
    if mask_large(yP,xP)
        idM=[idM,i];
    end
end

% Extract ss (cell size) and llInt (integrated track length) for cells within gradient mask
ss=ss(idM);
llInt=llInt(idM);
xT=xT(idM);
yT=yT(idM);

% Define bin edges based on ss
binEdges = linspace(min(ss), max(ss), 50);  % Adjust the number of bins as needed
mean_llInt = zeros(1, length(binEdges)-1);
std_llInt = zeros(1, length(binEdges)-1);

for i = 1:length(binEdges)-1
    indices = ss >= binEdges(i) & ss < binEdges(i+1); % indices of tracks within the mask
    counts(i) = sum(indices); % N added
    if any(indices)
        mean_llInt(i) = mean(llInt(indices));
        std_llInt(i) = std(llInt(indices))/sqrt(sum(indices));  % Standard error   
    end
end

% Plotting mean cell size vs integrated trajectory length with binning
figure;
scatter(binEdges(1:end-1), mean_llInt, 100, counts, 'filled'); % N added

xlabel('Mean Cell Size (ss)');
ylabel('Integrated Trajectory Length (llInt)');
title('Mean Cell Size vs Integrated Trajectory Length, Binned');
set(gca, 'fontsize', 21); 

colormap('jet'); % N added

%  Error Bars N added
hold on;
errorbar(binEdges(1:end-1), mean_llInt, std_llInt, 'vertical', 'k', 'LineStyle', 'none');

%Appearance - N added
colorbar;
ylabel(colorbar, 'Count');