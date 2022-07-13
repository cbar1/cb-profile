close all 
clear
clc 

% Define the filter size we will use in step 2:
filtsize = 85;

% Creating test image 'im' by splicing together two built in images.
% Also zero-padding (adding zeros around the border) with half the 
% filter size (filtsize) we will use so that the filter could be 
% centered on any actual image pixel, including those near the border.
% 'coins.png' contains bright nickels and dimes on a dark background
% 'eight.tif' contains dark quarters on a bright background, so we invert it
% to match 'coins.png'
im1 = imread('coins.png');
[r,c] = size(im1);
im2 = imread('eight.tif');
[r2,c2] = size(im2);
filtsizeh = floor(filtsize/2);
im = zeros(r+r2+filtsize,c+filtsize);
im(filtsizeh+1:filtsizeh+r+r2,filtsizeh+1:filtsizeh+c) = [im1;255-im2(:,1:c)];
[r,c] = size(im);

% imagesc(im);colormap(gray);title('test image');axis equal;

%% 1. Localize the centroid of each coin
% Otsu threshold
msk = OtsuThreshold(im);
% figure; imagesc(msk); colormap(gray); title('Otsu'); axis equal;

% Dilate 9x9
msk_dil = imdilate(msk,ones(9,9));
% figure; imagesc(msk_dil); colormap(gray); title('Dilated'); axis equal;

% Erode 23x23
msk_dil_erd = imerode(msk_dil,ones(23,23));
% figure; imagesc(msk_dil_erd); colormap(gray); title('Eroded'); axis equal;

% Connected components to get centroids of coins:
cc = bwconncomp(msk_dil_erd);
props_struct = regionprops(cc);
centroid = zeros(length(props_struct),2);
component_size = zeros(length(props_struct),1);
for i=1:length(props_struct)
    centroid(i,:) = round(props_struct(i).Centroid);
    component_size(i) = props_struct(i).Area;
end


%% 2. Measure features for each coin using a bank of matching filters
% make matching filters to create features
% Define diameters to use for filters
dimediameter = 31;
quarterdiameter = 51;
nickeldiameter = 41;

% Use the MakeCircleMatchingFilter function to create matching filters for dimes, nickels, and quarters
% (This is in a separate Matlab grader problem. Save your work, 
%       complete the corresponding grader problem and embed the solution 
%       in the helper function list below.)
nickelfilter = MakeCircleMatchingFilter(nickeldiameter,filtsize);
dimefilter = MakeCircleMatchingFilter(dimediameter,filtsize);
quarterfilter = MakeCircleMatchingFilter(quarterdiameter,filtsize);

% Evaluate each of the 3 matching filters on each coin to serve as 3 feature measurements 
D = zeros(length(centroid),3);
for i=1:length(centroid)
    D(i,1) = corr(dimefilter(:),reshape(msk_dil_erd(centroid(i,2)-filtsizeh:...
        centroid(i,2)+filtsizeh,centroid(i,1)-filtsizeh:centroid(i,1)+filtsizeh),[filtsize*filtsize,1]));
    D(i,2) = corr(nickelfilter(:),reshape(msk_dil_erd(centroid(i,2)-filtsizeh:...
        centroid(i,2)+filtsizeh,centroid(i,1)-filtsizeh:centroid(i,1)+filtsizeh),[filtsize*filtsize,1]));
    D(i,3) = corr(quarterfilter(:),reshape(msk_dil_erd(centroid(i,2)-filtsizeh:...
        centroid(i,2)+filtsizeh,centroid(i,1)-filtsizeh:centroid(i,1)+filtsizeh),[filtsize*filtsize,1]));    
end

%% 3. Perform k-means clustering of features for unsupervised learning classifier, get correct classifications
rng(0);
totcount=[];
cls_init = kmeans(D, 3);

% relabel centroid classes based on average size of the objects in each class. smallest will be dime, next nickel, and largest quarter
class_ave_object_size = zeros(3, 1);
numCoinCls = length(unique(cls_init));
for i = 1:numCoinCls
    class_ave_object_size(i, 1) = mean(component_size(cls_init == i));
end 

%get indices to sort cls_init by size 
[~, classmap] = sort(class_ave_object_size, 'ascend');


%first indice in classmap holds the class number label corresponding to the smallest
%average coin size, etc. 
%redefine the smallest class number, given by the number in classmap(1), to be 1
%the second smallest coin to be 2, and the largest to 3 in the class vector
cls = zeros(size(cls_init));
for i = 1:numCoinCls
    cls(cls_init == classmap(i)) = i;
end

%% What index the classes have been relocated to 
classmap_indx = [find(classmap==1); find(classmap==2); find(classmap==3)];
%% Visualize mapping 
table(unique(cls_init), class_ave_object_size, classmap_indx, classmap, 'VariableNames',...
    {'init_cls', 'ave_sizes', 'mapping', 'classmap'})
table(cls_init, cls, 'VariableNames', {'Cls_init', 'Cls_mapped'})

%% Visualize the result & sum coin values
values = zeros(height(centroid), 1);
figure; imagesc(im);colormap(gray);title('test image');hold on;axis equal;
% plot circles around each coin with different color/diameter unique to each type and count the change
for i = 1:height(centroid)
    values(i) = AddCoinToPlotAndCount(centroid(i, 1), centroid(i, 2), cls(i));
end
title([num2str(totcount),' cents'])

total_sum = sum(values)/100; %dollars 
fprintf('Change: $%.2f', total_sum)








%% 
%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%
function [coinvalue,x_plot,y_plot,col] = AddCoinToPlotAndCount(x,y,cls)
% initialize radians for defining x_plot and y_plot using cos and sin functions
rads = 0:2*pi/32:2*pi;
% initialize parameters for radius and color of circle for each type of coin
radius = [22, 30, 40];
colorStrs = ['r', 'g', 'm'];
coinValues = [10, 5, 25];


x_plot = cos(rads).*radius(cls)+x;
y_plot = sin(rads).*radius(cls)+y;

coinvalue = coinValues(cls);
col = colorStrs(cls);
plot(x_plot,y_plot,col);
end 

function filter = MakeCircleMatchingFilter(diameter,filtsize)
filter = zeros(filtsize,filtsize);
radius = diameter/2;
c = (filtsize+1)/2;
for i=1:filtsize
    for j=1:filtsize
        if (i-c)*(i-c) + (j-c)*(j-c) <= radius*radius
            filter(i,j) = 1;
        end
    end
end
end

function [msk,thrsh] = OtsuThreshold(im)
hst = imhist(im);
res = otsuthresh(hst);
thrsh = res*255;
msk = im>thrsh;
end
