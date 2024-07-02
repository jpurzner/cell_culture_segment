function connectResStruct = nuc_cyt_quant_list(rawdir, varargin)
%%% handle inputs args %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display = 1;
if ismember('dispoff',varargin)
    display = 0; 
end 

%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO DO: 
% - paramaterize the cytoplasmic channel naming

% The second channel will be used for cyto quants

names = {'DAPI', 'GFP'};
%names = {'DAPI', 'GFP', 'TexasRed'};
%names = {'DAPI' 'FITC', 'Cy5'};
%names = {'DAPI' 'GFP' 'TexasRed', 'Cy5', 'Cy3'};
%names = {'DAPI' 'TexasRed'};
file_names = strcat(names, '.tif');
file_names = fullfile(rawdir, file_names);

rawIm = cellfun(@imread, file_names, 'UniformOutput', false);
rawIm = cellfun(@double, rawIm, 'UniformOutput', false);
% have to normalize so adaptthresh works
%rawIm = cellfun(@(rawIm) rawIm./65535, rawIm, 'UniformOutput', false);

%%% segmentation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(rawIm{1},1) > 2000
    nucr=10*2;
    debrisarea=40*(2^2);
    boulderarea=1000*(2^2);
    num_thresh=6;
    nuc_filter = 40;
    
else
    nucr=7;
    debrisarea=40;
    boulderarea=400;
    num_thresh=6;
    nuc_filter = 20; 
end 

%%% remove smears %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foregroundthresh=1000;
areathresh=20000*(2^2);
bgprctile=20;
[dapi,~]=removesmears_2(rawIm{1},foregroundthresh,areathresh,bgprctile);
%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_maskA=threshmask(dapi);
nuc_maskA=markershed(nuc_maskA,round(nucr*2/3));

nuc_maskA=imfill(nuc_maskA,'holes');
nuc_maskA=bwareaopen(nuc_maskA,debrisarea);
nuc_maskA=segmentdeflections_bwboundaries(nuc_maskA,nucr,debrisarea);
nuc_maskA=logical(nuc_maskA);
%%% Note border/large/warped objects for later removal %%%%%%%%%%%%%%%%%%%%
antiborder_mask=imclearborder(nuc_maskA);
border_mask=nuc_maskA-antiborder_mask;
antilargewarped_mask=excludelargeandwarped_3(nuc_maskA,boulderarea,0.8);
largewarped_mask=nuc_maskA-antilargewarped_mask;
badmask=border_mask | largewarped_mask;
goodmask=logical(nuc_maskA-badmask);


[numberedNuc1,numOfNucs]=bwlabel(goodmask);
nucData = regionprops(numberedNuc1,'PixelIdxList','PixelList', 'Centroid','Area', 'Eccentricity','Solidity', 'MajorAxisLength', 'MinorAxisLength');
nuc_sd = arrayfun(@(x) (std(nucData(x).PixelList(:,2))), 1:size(nucData,1));
nuc_range = arrayfun(@(x) (max(nucData(x).PixelList(:,2)) - min(nucData(x).PixelList(:,2))), 1:size(nucData,1));
nuc_sd = num2cell(nuc_sd);
nuc_range = num2cell(nuc_range);

[nucData.nuc_sd]  = nuc_sd{:};
[nucData.nuc_range] = nuc_range{:};

% filter for area 
keep_nuc = [nucData.Area] > nuc_filter;
nucData = nucData(keep_nuc);
numOfNucs = size(nucData,1);

keep_region_idx = find(keep_nuc); 
keep_label_mask = ismember(numberedNuc1, keep_region_idx);
keep_label = immultiply(keep_label_mask, numberedNuc1);

%%% segment cytoplasm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cytr = 20;
cytMarker = rawIm{2};

% 65535 is the max value for a 16-bit image. adaptthresh requires values
% between 0-1, so I normalized w/ this value. This is fixed downstream to 
% get accurate pixel intensities. Sorry, I'll fix this later, just wanted 
% to get it done !! 
cyt_mask=threshmask_adaptive(cytMarker, 0.5, 101, 65535);
cyt_mask=markershed(cyt_mask,round(100));
cyt_mask=imfill(cyt_mask,'holes');
cyt_mask=bwareaopen(cyt_mask,debrisarea);
cyt_mask=segmentdeflections_bwboundaries(cyt_mask,cytr,debrisarea);
cyt_mask=logical(cyt_mask);

cyt_mask = bwareaopen(cyt_mask, 300); 
cyt_mask = imfill(cyt_mask, 'holes'); 
se = strel('disk', 2);
cyt_mask = imopen(cyt_mask, se); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% voronoi%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imageSize = size(rawIm{2});

centroids = cat(1, nucData.Centroid);
[V, C] = voronoin(centroids);

labelVoronoi = zeros(imageSize(1), imageSize(2)); 

for k = 1:length(C)
    if any(C{k} == 1)
        vertices = V(C{k}(C{k} ~= 1), :); 
        vertices(:,1) = min(max(vertices(:,1), 1), imageSize(2));
        vertices(:,2) = min(max(vertices(:,2), 1), imageSize(1));
    else
        vertices = V(C{k}, :);
    end
    mask = poly2mask(vertices(:,1), vertices(:,2), imageSize(1), imageSize(2));
    labelVoronoi(mask) = k;
end

cytoMaskA = labelVoronoi;
logicalMask = (~(labelVoronoi > 0 & ~nuc_maskA & cyt_mask)); %labels retained
cytoMaskA(logicalMask) = 0;

%wholeCellMask = labelVoronoi;
%logicalMaskB = (~(labelVoronoi > 0 & cyt_mask | nuc_maskA)); %labels retained
%wholeCellMask(logicalMaskB) = 0;


cytData = regionprops(cytoMaskA, 'PixelIdxList', 'PixelList');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% saving data to arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_im_cell = repmat({zeros(1,numOfNucs)},1,size(rawIm,2));
median_im_cell = repmat({zeros(1,numOfNucs)},1,size(rawIm,2));


% need to vectorize this nested loop 
for nucInd=1:numOfNucs  %for each cell in the image
    for i=1:size(rawIm, 2)
        sum_im_cell{i}(nucInd) = sum(rawIm{i}(nucData(nucInd).PixelIdxList));
        median_im_cell{i}(nucInd) = median(rawIm{i}(nucData(nucInd).PixelIdxList));
    end 
end

sum_im_cyt = zeros(1, numOfNucs);
median_im_cyt = zeros(1, numOfNucs);

for cytInd=1:numOfNucs
    if any(cytoMaskA(:) == cytInd)
        cytPixels = cytData(cytInd).PixelIdxList;
        sum_im_cyt(cytInd) = sum(rawIm{2}(cytPixels));
        median_im_cyt(cytInd) = median(rawIm{2}(cytPixels));
    else
        sum_im_cyt(cytInd) = 0;
        median_im_cyt(cytInd) = 0;
    end
end
%%%%% final structs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_struct = cell2struct(sum_im_cell, strcat(names, '_sum'),2);
median_struct = cell2struct(median_im_cell, strcat(names, '_median'),2);

sum_struct.(strcat(names{2}, '_cyt_sum')) = sum_im_cyt;
median_struct.(strcat(names{2}, '_cyt_median')) = median_im_cyt;
%sum_struct.TexasRed_whole_sum = sum_struct.TexasRed_sum + sum_im_cyt;
%median_struct.TexasRed_whole_median = (median_struct.TexasRed_median + median_im_cyt) / 2;

im_struct = catstruct(sum_struct, median_struct);
connectResStruct.nucData = nucData;
connectResStruct.imData = im_struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%sp_vals = catstruct(sp_vals, im_struct);
%process_length(goodmask, rawIm{3}, connectResStruct.nucData, keep_label)

% labeling connected shapes on the nuclear mask
[L, num] = bwlabel(nuc_maskA);

%%% segmentation visualizations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if display 
  
    %nuclear segmentations 
    coloredLabels = label2rgb (keep_label, 'hsv', 'k', 'shuffle');
    figure, imshow(coloredLabels)

    %coloredVoronoi = label2rgb (labelVoronoi, 'hsv', 'k', 'shuffle');
    %figure, imshow(coloredVoronoi)

    coloredcytoplasm = label2rgb (cytoMaskA, 'hsv', 'k', 'shuffle');
    figure, imshow(coloredcytoplasm);

    %wholeCell = label2rgb (wholeCellMask, 'hsv', 'k', 'shuffle');
    %figure, imshow(wholeCell);

  %  outlined DAPI segments
  %  extractmask=bwmorph(goodmask > 0,'remove');
  %  RGB=imadjust(mat2gray(dapi));
  %  RGB(:,:,2)=extractmask;
  %  RGB(:,:,3)=0;
  %  figure,imshow(RGB);
   
    %cytoplasmic segmentations
    cyto_image = mat2gray(rawIm{2}); 
    cyto_mask_RGB = repmat(cyto_image, [1, 1, 3]); 
    cyto_mask_RGB(:,:,2) = cyto_mask_RGB(:,:,2) + (cyt_mask * 0.5);
    cyto_mask_RGB(cyto_mask_RGB > 1) = 1;
    figure, imshow(cyto_mask_RGB);
    
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% smear removal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [correctedimage,stainflag]=removesmears_2(rawimage,foregroundthresh,areathresh,bgprctile)
smear_mask=rawimage>foregroundthresh;
smear_mask=bwareaopen(smear_mask,areathresh);
if max(smear_mask(:))>0
    stainflag=1;
else
    stainflag=0;
end
smear_mask=imdilate(smear_mask,strel('disk',100,0));
foreground=rawimage(~smear_mask);
background=prctile(foreground(:),bgprctile);
correctedimage=rawimage;
correctedimage(smear_mask)=background;
end

function [correctedimage,stainflag]=removesmears_3(rawimage,foregroundthresh_2,areathresh_2,bgprctile_2)
smear_mask=rawimage>foregroundthresh_2;
smear_mask=bwareaopen(smear_mask,areathresh_2);
if max(smear_mask(:))>0
    stainflag=1;
else
    stainflag=0;
end
smear_mask=imdilate(smear_mask,strel('disk',100,0));
foreground=rawimage(~smear_mask);
background=prctile(foreground(:),bgprctile_2);
correctedimage=rawimage;
correctedimage(smear_mask)=background;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mask=threshmask(image)
normlog=mat2gray(log(image));
thresh=graythresh(normlog);
mask=im2bw(normlog,thresh);
mask=imfill(mask,'holes');
end

function cytmask = threshmask_adaptive(image, sensitivity, neighborhoodSize, normFactor)
if size(image, 3) == 3
    image = rgb2gray(image);
end
imageNormalized = image ./ normFactor;
imageNormalized = im2double(imageNormalized);
T = adaptthresh(imageNormalized, sensitivity, 'NeighborhoodSize', neighborhoodSize);
cytmask = imbinarize(imageNormalized, T);
cytmask = imfill(cytmask, 'holes');
end

function maskwatershed=markershed(mask,eroderadius)
valleys=-bwdist(~mask);
basins=imerode(mask,strel('disk',eroderadius,0));
bigvalleys=bwdist(mask);
outerridges=watershed(bigvalleys);
outerridges=outerridges==0;
finalvalleys=imimposemin(valleys,basins | outerridges);
finalridges=watershed(finalvalleys);
mask(finalridges==0)=0;
maskwatershed=mask;
end

function mask=segmentdeflections_bwboundaries(mask,nucr,debrisarea)
[B,L]=bwboundaries(mask,'noholes');
obnum=max(L(:));
bordermask=zeros(size(mask));
for ci=1:obnum
    orderedset=B{ci};
    % If using bwboundaries, reverse the order.
    orderedset=[orderedset(end:-1:1,2) orderedset(end:-1:1,1)];
    bordermask=splitdeflections_4_bwboundaries(orderedset,bordermask,nucr);
end
mask=mask & ~bordermask;
mask=~bwmorph(~mask,'diag');
mask=bwareaopen(mask,debrisarea);
end

function [bordermask,bridgeflag]=splitdeflections_4_bwboundaries(orderedset,bordermask,nucr)
bridgeflag=0; %returned as 1 if deflections are bridged
nucr=round(nucr/4)*4; %make sure nucr is a multiple of 4
perilength=size(orderedset,1);
%%% detect deflection vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vIdx=getdeflections(orderedset,nucr); %returns boundary indices
%%% count vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vnum=length(vIdx);
if vnum<2
    return; %if less than two vertices are detected, exit function
end
%%% calculate perimeter distance between adjacent vertices %%%%%%%%%%%%%%%%
periIdx=vIdx;
periIdxadj1=[periIdx(2:end);perilength+periIdx(1)];
pairperi1=periIdxadj1-periIdx;
%%% pair and bridge vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while vnum>=2
    vpos=orderedset(vIdx,:);
    %%% Determine adjacent vertices that define the highest curvature %%%%%
    vposadj1=[vpos(2:end,:);vpos(1,:)];
    pair1=vposadj1-vpos;
    pairdist1=sqrt(sum(pair1.^2,2));
    curvature1=pairperi1./pairdist1;
    [bestcurve,curve1idx]=sort(curvature1);
    % If best curvature is too low, stop segmenting this object.
    if bestcurve(end)<2
        break
    end
    bestcurveidx=curve1idx(end);
    if bestcurveidx==vnum
        bestcurveidxadj=1;
    else
        bestcurveidxadj=bestcurveidx+1;
    end
    % If this point is reached, a split will be performed, so mark it.
    bridgeflag=1;
    %%% Bridge the vertices defining the best curvature %%%%%%%%%%%%%%%%%%%
    [bx,by]=bridge(vpos(bestcurveidx,:),vpos(bestcurveidxadj,:));
    % [bx,by] can equal NaN if the vertices were the same position (occurs
    % when two objects share a single coordinate. If this is the case,
    % skip mapping the bridge and proceed with perimeter update and vertex
    % removal.
    if ~isnan(bx)
        for bci=1:length(bx)
            %bridgeflag=1;
            bordermask(by(bci),bx(bci))=1;
        end
    end
    %%% assign new perimeter distances & remove old vertices %%%%%%%%%%%%%%
    previdx=bestcurveidx-1;
    if previdx==0
        previdx=vnum;
    end
    % Given vertex 3-4 gave best curvature and is now bridged, define the
    % perimeter from vertex 2 to 5: p(2-5)=p(2-3)+bridge+p(4-5).
    pairperi1(previdx)=pairperi1(previdx)+length(bx)-1+pairperi1(bestcurveidxadj);
    % Remove the vertices and perimeters of the vertices defining the best
    % curve.
    vIdx([bestcurveidx,bestcurveidxadj])=[];
    pairperi1([bestcurveidx,bestcurveidxadj])=[];
    vnum=length(vIdx);
end
end

function [bx,by]=bridge(vpos1,vpos2)
lengthx=vpos2(1)-vpos1(1);
lengthy=vpos2(2)-vpos1(2);
longerside=max([abs(lengthx) abs(lengthy)]);
stepx=lengthx/longerside;
stepy=lengthy/longerside;
bx=zeros(longerside+1,1);by=zeros(longerside+1,1);
for bs=0:longerside
    bx(1+bs)=vpos1(1)+round(bs*stepx);
    by(1+bs)=vpos1(2)+round(bs*stepy);
end
end

function vIdx=getdeflections(orderedset,nucr)
%%% set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
offsetshort=nucr/4;
if offsetshort==1
    offsetshort=2;
end
offsetlong=2*offsetshort;
gradientoffsetshort=offsetshort;
gradientoffsetlong=offsetlong;
shortgradthresh = pi/6;
longgradthresh = pi/6;
%%% calculate angular deflection at each point of the boundary %%%%%%%%%%%%
perilength=size(orderedset,1);
%%%%%%% short steps %%%%%%%%%%%%%%%%%%
orderedsetoffsetshort=[orderedset(offsetshort+1:end,:);orderedset(1:offsetshort,:)];
shortdiff=orderedsetoffsetshort-orderedset;
shortgrad=atan2(shortdiff(:,2),shortdiff(:,1));   %angle in radians
shortgradoffset=[shortgrad(gradientoffsetshort+1:end,:);shortgrad(1:gradientoffsetshort,:)];
shortgraddiff=shortgradoffset-shortgrad;
shortgraddiff=shortgraddiff+2*pi*(shortgraddiff<0);  %account for 4 quadrants
%%%%%%% long steps %%%%%%%%%%%%%%%%%%%
orderedsetoffsetlong=[orderedset(offsetlong+1:end,:);orderedset(1:offsetlong,:)];
longdiff=orderedsetoffsetlong-orderedset;
longgrad=atan2(longdiff(:,2),longdiff(:,1));
longgradoffset=[longgrad(gradientoffsetlong+1:end,:);longgrad(1:gradientoffsetlong,:)];
longgraddiff=longgradoffset-longgrad;
longgraddiff=longgraddiff+2*pi*(longgraddiff<0);
%%% find deflections above threshold %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% short steps %%%%%%%%%%%%%%%%%%
shortgraddiff(shortgraddiff>=pi)=0;  %exclude
vIdxmaskshort=shortgraddiff>shortgradthresh;
%%%%%%% long steps %%%%%%%%%%%%%%%%%%%
vIdxmasklong=longgraddiff>longgradthresh & longgraddiff<pi;
vIdxmasklong=[zeros(offsetlong,1);vIdxmasklong(1:end-offsetlong)];
vIdxmasklong=imdilate(vIdxmasklong,strel('square',1+nucr/2));
%%% find local maxima of short steps %%%%%%%%%%%%%%%%%%%%%%%
vIdxmaskshort=imclose(vIdxmaskshort,strel('square',3));  %join proximal deflection islands
vIdxobs=regionprops(bwlabel(vIdxmaskshort),'PixelIdxList');
maxmask=zeros(size(vIdxmaskshort));
for rpc=1:length(vIdxobs)
    pix=vIdxobs(rpc).PixelIdxList;
    [~,index]=max(shortgraddiff(pix));
    maxmask(pix(index)+offsetshort)=1;
end
maxmask=maxmask(1:perilength);  %remove any overhang
%%% find coincidence of long mask & local maxima of short mask %%%%%%%%%%%%
vIdxmask=vIdxmasklong & maxmask;
vIdx=find(vIdxmask);
end

function nuc_mask=excludelargeandwarped_3(nuc_mask,boulderarea,soliditythresh)
antimask=bwareaopen(nuc_mask,boulderarea);
nuc_mask=nuc_mask-antimask;
nuc_label=bwlabel(nuc_mask);
nuc_solidity=cell2mat(struct2cell(regionprops(nuc_label,'Solidity')));
warpedobjects=find(nuc_solidity<soliditythresh); %default 10xbin1:0.9
nuc_label(ismember(nuc_label,warpedobjects))=0;
nuc_mask=nuc_label>0;
end

function biascorr_nuc(nucraw, nucr)
    
    
    tilenum=31; %number of tiles (e.g. 10 --> 10x10) for calculating background within each tile
    prctilethresh=50; %the intensity percentile to take within each tile

    
    %%% subtract CMOS offset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nucraw=nucraw-cmosoffset;
    %%% generate dilated foreground mask to block out cells %%%
    nucmask=threshmask(nucraw,3);
    nucmask=imdilate(nucmask,strel('disk',round(nucr),0));
    blurnan=imfilter(nucraw,fspecial('disk',3),'symmetric');
    blurnan(nucmask)=NaN;
    %%% calculate background in each tile %%%%%%%%%%%%%%%%%%%%%
    bgblock=blockpercentile_blockimage(blurnan,tilenum,prctilethresh);
    %%% normalize each tile against the center tile %%%%%%%%%%%
    midrc=ceil(tilenum/2);
    refval=bgblock(midrc,midrc);
    if ~isnan(refval)
        bgblocknorm=bgblock/refval;
        biasstack=cat(3,biasstack,bgblocknorm);
    end
end 

end