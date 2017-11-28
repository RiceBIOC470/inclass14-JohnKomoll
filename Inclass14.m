%Inclass 14
%GB comments
1 100 
2a 100 
2b 100

%Work with the image stemcells_dapi.tif in this folder

% Read in file
file = 'stemcells_dapi.tif';
img = imread(file);
img = im2double(img);
imshow(img, [])

% (1) Make a binary mask by thresholding as best you can

% Define mask with threshold
threshold = 0.0048;
mask = img > threshold;
figure
imshow(mask)

% (2) Try to separate touching objects using watershed. Use two different
% ways to define the basins. (A) With erosion of the mask (B) with a
% distance transform. Which works better in this case?

% A. Use erosion of the mask

% Figure out which centroids represent multiple cells
CC = bwconncomp(mask);
cell_data = regionprops(CC, 'Area');
area = [cell_data.Area];
fused_cells = area > mean(area) + 0.1*std(area);  % Many clusters exist, so take many larger chunks
sublist = CC.PixelIdxList(fused_cells);
sublist = cat(1, sublist{:});
fusedMask = false(size(mask));
fusedMask(sublist) = 1;
figure
imshow(fusedMask, 'InitialMagnification', 'fit')

% Erode the multi-cell clusters to get centers and outside region
nucmin = imerode(fusedMask, strel('disk', 7));
figure
imshow(nucmin, 'InitialMagnification', 'fit')
outside = ~imdilate(fusedMask, strel('disk', 1));

% Define the basins for watershed
basin = imcomplement(bwdist(outside));
basin = imimposemin(basin, nucmin | outside);
pcolor(basin);
shading flat;

% Perform watershed and combine nuclei with nonproblematic ones
L = watershed(basin);
newmask = L > 1 | (mask - fusedMask);
figure
imshow(newmask, 'InitialMagnification', 'fit')

% B. Use a distance transformation
D = -bwdist(~mask);
D(~mask) = -Inf;
L = watershed(D);
L(~mask) = 0;
rgb = label2rgb(L, 'jet', [.5 .5 .5]);
figure
imshow(rgb,'InitialMagnification', 'fit')
figure
imshow(L > 1, 'InitialMagnification', 'fit')

% When comparing Figure 5 (Part A) and Figure 7 (Part B), it is clear that
% Part A, using the erosion mask, was much more successful. Part B
% oversegments the nuclei, while Part A only leaves a few clusters of cells
% undivided.

