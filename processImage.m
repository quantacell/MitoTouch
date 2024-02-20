% ----------------------------------------
% Copyright Clariant, CNRS, Université de Montpellier, IRD, EPHE
% Contributors: Abdel Aouacheria, Victor Racine
% Contacts : abdel.aouacheria@umontpellier.fr, victor.racine@quantacell.com 

%  
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License v3.0 only.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
%  
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see https://spdx.org/licenses/GPL-3.0-only.html
% ----------------------------------------


function [data, dataHeaders, T, masks, maskNames]=processImage(mergeFilename, outputFilename, nucleiThreshold, cellThreshold, mitoThrehsold, erodeDist, filterBorder, isCellBased)

if ~exist('mergeFilename','var')
    mergeFilename='13-1- MERGE.tif';
    outputFilename=  'fibro\out';
end
if numel(outputFilename)>0
    [~, condition, ~]=fileparts(outputFilename);
else
    condition='condition';
end


if strcmpi(mergeFilename(end-4:end), '.tiff')
    mito=imread(mergeFilename);
    filenameChan2=strrep(mergeFilename, 'ch1', 'ch2');
    nuclei=imread(filenameChan2);
    filenameChan3=strrep(filenameChan2, 'ch2', 'ch3');
    cell=imread(filenameChan3);

else
    nuclei=double(imread(mergeFilename, 3));
    cell=double(imread(mergeFilename, 2));
    mito=double(imread(mergeFilename, 1));
end
im={};
imNames={};
im{1}=nuclei;
imNames{1}='nuclei';

im{2}=cell;
imNames{2}='cell';

masks{1}=segmentNuclei(nuclei,nucleiThreshold, 1, 3000, 1);


maskNames{1}='Nuclei';
[masks{2}, masks{1}]=segmentMembraneIntensity(masks{1}, cell, 3, cellThreshold, 512);
maskNames{2}='CellsMask';

if erodeDist>0
    masks{2}=imerode(masks{2}, strel('disk', erodeDist));
    masks{2}=removeSmallIntensity(bwlabel(masks{2}>0), masks{1}, 0);
    masks{1}(masks{2}==0)=0;
    masks{1}=(masks{1}>0).*double(masks{2});
end

if filterBorder>0
    masks=filterBorderProcessing(masks, filterBorder);
    
end



im{3}=mito;
imNames{3}='mito';
[masks{3}]=segmentMitochondrion(masks{2}, mito, mitoThrehsold/2,1, 50);
maskNames{3}='MitoClusterMask';
T=table();
if numel(outputFilename)>0
    [data, dataHeaders, masks, maskNames]=measureMito(masks, maskNames, im, imNames, isCellBased);
else
    masks{4}=double(bwskel(getImageInList(masks, maskNames, 'MitoClusterMask')>0)).*getImageInList(masks, maskNames, 'CellsMask');
    maskNames{4}='MitoSkelMask';
    data=[];
    dataHeaders={};
end

T.Well=repmat({condition}, size(data, 1), 1);
T.Field=repmat({''}, size(data, 1), 1);
T.WellField=repmat({''}, size(data, 1), 1);
 T.SubFolder1=repmat({''}, size(data, 1), 1);
 T.SubFolder2=repmat({''}, size(data, 1), 1);
for t=1:numel(dataHeaders)
    name=dataHeaders{t};
    if (name(1)=='1') || (name(1)=='2') || (name(1)=='3') 
        name=['ch', name];
    end
    T.(name)=data(:,t);
end
if numel(outputFilename)>0
    imagesc(cell);
    drawRegions(masks{3}, [], false, [], false);
    colormap gray;
    daspect([1 1 1]);
    drawnow;

    for m=1:numel(maskNames)
        imwrite(uint16(masks{m}), [outputFilename, '_mask', maskNames{m},'.tif']);
    end

    for m=1:numel(imNames)
        imwrite(uint16(im{m}), [outputFilename, '_im', imNames{m},'.tif']);
    end

    imwrite(uint8(double(masks{3}>0)+1), [outputFilename, '_mask', maskNames{3},'.png']);
else
    data=[];
    dataHeaders={};
end

end

function masks=filterBorderProcessing(masks, filterBorder)

masks2B=masks{2};
masks2B(2:end-1, 2:end-1)=0;
stats=regionprops(masks{2}, 'perimeter');
statsB=regionprops(masks2B, 'area');
toKeep=[];
for i=1:numel(stats)
    if i<=numel(statsB)
        if (stats(i).Perimeter-statsB(i).Area)/stats(i).Perimeter>=filterBorder
            toKeep=cat(1, toKeep, i);
        end
    end
end
for i=1:numel(masks)
    masks{i} = bwlabel(ismember(masks{i},toKeep));
end
end

function imOut=removeJonctions(imIn)
imInd=imdilate(imIn, ones(3));
jonct=(imInd~=imIn);
jonct(imIn==0)=0;
imOut=imIn;
imOut(jonct)=0;
end


function [data, dataHeaders]=addValue(data, dataHeaders, val, header, cellIndex)
indexHeader=-1;
for i=1:numel(dataHeaders)
    if strcmpi(dataHeaders{i}, header)
        indexHeader=i;
    end
end
if indexHeader==-1
    dataHeaders=cat(2, dataHeaders, {header});
    indexHeader=numel(dataHeaders);
end

data(cellIndex, indexHeader)=val;

end

function [data, dataHeaders, masks, maskNames]=measureMito(masks, maskNames, im, imName, isCellBased)
dataHeaders={};%'CellId', 'Cell_Area', 'Nucleus_Area', 'X', 'Y', 'Count', 'MeanArea', 'Circularity', };

data=[];

maskCells=getImageInList(masks, maskNames, 'CellsMask');
maskNuclei=getImageInList(masks, maskNames, 'Nuclei');
maskMitoCluster=getImageInList(masks, maskNames, 'mitoClusterMask');
cellCount=max(max(maskCells));
imMito=getImageInList(im, imName, 'mito');
skelFinal=double(bwskel(maskMitoCluster>0)).*maskCells;
if isCellBased
    

    for i=1:cellCount
        %[data, dataHeaders]=addValue(data, dataHeaders, i, 'condition', i);
        [data, dataHeaders]=addValue(data, dataHeaders, i, 'CellId', i);
        maskCell0=(maskCells==i);
        
        stats=regionprops(maskCell0, imMito, 'Area', 'MinorAxisLength', ...
            'MajorAxisLength', 'MeanIntensity', ...
            'EquivDiameter', 'MaxIntensity', 'Perimeter');
        
        
        [data, dataHeaders]=addValue(data, dataHeaders, stats(1).Area, 'Cell_Area', i);
        [data, dataHeaders]=addValue(data, dataHeaders, stats(1).MeanIntensity, 'Cell_MeanIntensity', i);
        [data, dataHeaders]=addValue(data, dataHeaders, stats(1).MaxIntensity, 'Cell_MaxIntensity', i);
        [data, dataHeaders]=addValue(data, dataHeaders, stats(1).Perimeter, 'Cell_Perimeter', i);
        [data, dataHeaders]=addValue(data, dataHeaders, stats(1).MinorAxisLength/stats(1).MajorAxisLength, 'Cell_Compaction', i);
        [data, dataHeaders]=addValue(data, dataHeaders, stats(1).EquivDiameter/stats(1).Perimeter/4*pi, 'Cell_Roundness', i);
        
        
        [n,r] = boxcount(maskCell0);
        s=-gradient(log(n))./gradient(log(r));
        %[data, dataHeaders]=addValue(data, dataHeaders, s(2), 'MitoCluster_Fractal2', i);
        [data, dataHeaders]=addValue(data, dataHeaders, s(4), 'MitoCluster_Fractal8', i);
        [data, dataHeaders]=addValue(data, dataHeaders, s(6), 'MitoCluster_Fractal32', i);
        [data, dataHeaders]=addValue(data, dataHeaders, s(8), 'MitoCluster_Fractal64', i);
        maskNucleus0=(maskNuclei==i);
        maskMitoCluster0=bwlabel(maskMitoCluster.*maskCell0);
        stats=regionprops(maskMitoCluster0, imMito, 'Area', 'MinorAxisLength', ...
            'MajorAxisLength', 'EulerNumber', 'MeanIntensity', ...
            'EquivDiameter', 'MaxIntensity', 'Perimeter', 'Solidity');
        [data, dataHeaders]=addValue(data, dataHeaders, numel([stats.Area]), 'MitoCluster_Count', i);
        [data, dataHeaders]=addValue(data, dataHeaders, mean([stats.Area]), 'MitoCluster_Area', i);
        [data, dataHeaders]=addValue(data, dataHeaders, mean([stats.MinorAxisLength]./[stats.MajorAxisLength]), 'MitoCluster_Compaction', i);
        [data, dataHeaders]=addValue(data, dataHeaders, mean([stats.MajorAxisLength]./[stats.MinorAxisLength]), 'MitoCluster_Elongation', i);
        [data, dataHeaders]=addValue(data, dataHeaders, mean([stats.EquivDiameter]./[stats.Perimeter]/4*pi), 'MitoCluster_Roundness', i);
        [data, dataHeaders]=addValue(data, dataHeaders, mean([stats.EulerNumber]), 'MitoCluster_EulerNumber', i);
        [data, dataHeaders]=addValue(data, dataHeaders, mean([stats.EulerNumber]), 'EulerNumber', i);
        [data, dataHeaders]=addValue(data, dataHeaders, mean([stats.MeanIntensity]), 'MitoCluster_MeanIntensity', i);
        [data, dataHeaders]=addValue(data, dataHeaders, mean([stats.MaxIntensity]), 'MitoCluster_MaxIntensity', i);

        [data, dataHeaders]=addValue(data, dataHeaders, mean([stats.Perimeter]), 'MitoCluster_Perimeter', i);
        [data, dataHeaders]=addValue(data, dataHeaders, mean([stats.Solidity]), 'MitoCluster_Solidity', i);

        skel=bwskel(maskMitoCluster0>0);
        dist=bwdist(maskMitoCluster0==0);

        [data, dataHeaders]=addValue(data, dataHeaders, mean(dist(skel)), 'Skel_Width', i);
        [data, dataHeaders]=addValue(data, dataHeaders, sum(sum(skel)), 'Skel_Length', i);
        branchpoints=bwmorph(skel>0, 'branchpoints');
        [data, dataHeaders]=addValue(data, dataHeaders, sum(sum(branchpoints)), 'Skel_BranchPointsCount', i);
        endpoints=bwmorph(skel>0, 'endpoints');
        [data, dataHeaders]=addValue(data, dataHeaders, sum(sum(endpoints)), 'Skel_EndPointsCount', i);

        [data, dataHeaders]=addValue(data, dataHeaders, sum(sum(branchpoints))/(sum(sum(endpoints))+sum(sum(branchpoints))), 'Skel_BranchPointsEndPointsRatio', i);

        %branchpoints=imdilate(branchpoints, ones(3));
        skel2=skel;
        skel2(imdilate(branchpoints>0, ones(3)))=0;
        skel2=bwlabel(skel2);

        stats2=regionprops(skel2, 'Area', 'MinorAxisLength', 'MajorAxisLength', 'EquivDiameter', 'Perimeter');
        [data, dataHeaders]=addValue(data, dataHeaders, mean([stats2.MinorAxisLength]./[stats2.MajorAxisLength]), 'Mito_Compaction', i);
        [data, dataHeaders]=addValue(data, dataHeaders, mean([stats2.MajorAxisLength]./[stats2.MinorAxisLength]), 'Mito_Elongation', i);
        ed=[stats2.EquivDiameter];
        pe=[stats2.Perimeter];
        [data, dataHeaders]=addValue(data, dataHeaders, mean(ed(pe>0)./pe(pe>0)/4*pi), 'Mito_Roundness', i);
        [data, dataHeaders]=addValue(data, dataHeaders, mean([stats2.Area]), 'Mito_Length', i);

        dist=bwdist(maskCell0==0);
        [data, dataHeaders]=addValue(data, dataHeaders, mean(dist(skel2>0)), 'DistToCellMembrane', i);

        dist2=bwdist(maskNucleus0>0);

        [data, dataHeaders]=addValue(data, dataHeaders, mean(dist2(skel>0)), 'DistToNuclei', i);

        dist3=dist./(dist+dist2);
        [data, dataHeaders]=addValue(data, dataHeaders, mean(dist3(skel>0)), 'RatioDistMemb0Nucl1', i);
    end
    masks{4}=skelFinal;
    maskNames{4}='MitoSkelMask';
else
    maskCells=maskCells>0;
    maskNuclei=maskNuclei>0;
    maskMitoCluster=maskMitoCluster>0;
    branchpoints=bwmorph(skelFinal>0, 'branchpoints');
    branchpoints=imdilate(branchpoints, ones(3));
    skelFinal(branchpoints>0)=0;
    skelFinal=bwlabel(skelFinal);
    dist=bwdist(skelFinal>0);
    dist(maskMitoCluster==0)=-Inf;
    wat=watershed(dist);
    wat(maskMitoCluster==0)=0;
    maskMito=bwlabel(removeSmall(wat, 10));
    mitoCount=max(max(maskMito));
    
    stats=regionprops(maskMito, imMito, 'Area', 'MinorAxisLength', ...
            'MajorAxisLength', 'EulerNumber', 'MeanIntensity', ...
            'EquivDiameter', 'MaxIntensity', 'Perimeter');
    [data, dataHeaders]=addValue(data, dataHeaders, 1:mitoCount, 'MitoId', 1:mitoCount);
   
    [data, dataHeaders]=addValue(data, dataHeaders, [stats.Area], 'Mito_Area', 1:mitoCount);
    [data, dataHeaders]=addValue(data, dataHeaders, [stats.MinorAxisLength]./[stats.MajorAxisLength], 'Mito_Compaction', 1:mitoCount);
    [data, dataHeaders]=addValue(data, dataHeaders,[stats.EquivDiameter]./[stats.Perimeter]/4*pi, 'Mito_Roundness', 1:mitoCount);

    [data, dataHeaders]=addValue(data, dataHeaders,[stats.MeanIntensity], 'Mito_MeanIntensity', 1:mitoCount);
    [data, dataHeaders]=addValue(data, dataHeaders, [stats.MaxIntensity], 'Mito_MaxIntensity', 1:mitoCount);

    [data, dataHeaders]=addValue(data, dataHeaders, [stats.Perimeter], 'Mito_Perimeter', 1:mitoCount);


    dist=bwdist(maskCells==0);
    dist2=bwdist(maskNuclei>0);
    dist3=dist./(dist+dist2);
    stats=regionprops(maskMito, dist, 'MeanIntensity');
    stats2=regionprops(maskMito, dist2, 'MeanIntensity');
    stats3=regionprops(maskMito, dist3, 'MeanIntensity');
    [data, dataHeaders]=addValue(data, dataHeaders, [stats.MeanIntensity], 'DistToCellMembrane', 1:mitoCount);
    [data, dataHeaders]=addValue(data, dataHeaders, [stats2.MeanIntensity], 'DistToNuclei', 1:mitoCount);
    [data, dataHeaders]=addValue(data, dataHeaders, [stats3.MeanIntensity], 'RatioDistMemb0Nucl1', 1:mitoCount);
    
    masks{4}=skelFinal;
    maskNames{4}='MitoSkelMask';
    
    masks{5}=maskMito;
    maskNames{5}='MitoMask';
end


end

function [data, dataHeaders]=measureMasks(masks, maskNames, im, imName)
dataHeaders={'CellId', 'Cell_Area', 'Nucleus_Area', 'X', 'Y'};%, ...
    %'MeanInt_Chan2', 'MaxInt_Chan2', 'SumInt_Chan2', 'MeanInt_ChanAlpha', 'MaxInt_ChanAlpha','SumInt_ChanAlpha', ...
    %'%age in MaskAlpha'
    %};

maskCells=masks{2};
maskNuclei=masks{1};
cellCount=max(max(maskCells));
nbFeatures=numel(dataHeaders)+6+numel(imName);

data=zeros(cellCount, nbFeatures);
data(:,1)=(1:cellCount)';
if cellCount==0
    return;
end
statsReg=regionprops(maskNuclei, 'Area', 'Centroid');
statsRegCell=regionprops(maskCells, 'Area');
data(:, 2)=[statsRegCell.Area]';
data(:, 3)=[statsReg.Area]';
centroids=reshape([statsReg.Centroid]', 2, numel(statsReg));
data(:, 4)=centroids(1, :)';
data(:, 5)=centroids(2, :)';

cpt=numel(dataHeaders);
for c=1:numel(imName)
    
    statsReg=regionprops(maskCells, im{c}, 'Area', 'MeanIntensity', 'MaxIntensity');
    cpt=cpt+1;
    dataHeaders{cpt}=[imName{c}, '_mean'];
    data(:, cpt)=[statsReg.MeanIntensity]';
    if c<=3
        cpt=cpt+1;
        dataHeaders{cpt}=[imName{c}, '_max'];
        data(:, cpt)=[statsReg.MaxIntensity]';
        cpt=cpt+1;
        dataHeaders{cpt}=[imName{c}, '_sum'];
        data(:, cpt)=[statsReg.Area]'.*[statsReg.MeanIntensity]';
    end
end

end

function [maskCells, maskNuclei]=segmentMembraneIntensity(maskNuclei, intMap, sigma, threshold, sizeLimit)
originalSize=size(maskNuclei);
if ~exist('sizeLimit', 'var')
    sizeLimit=256;
end
if max(originalSize)>sizeLimit
    ratio=sizeLimit/max(originalSize);
    maskNuclei=imresize(maskNuclei, ratio, 'nearest');
    intMap=imresize(intMap, ratio);
    sigma=sigma*ratio;
else
    ratio=1;
end
intMap0=intMap;
intMap=imfilter(double(intMap), fspecial('gaussian', round((sigma)*6+1), sigma), 'symmetric', 'same'); 
intMap=intMap-max(max(intMap));
intensities = intMap;
intensities(maskNuclei>0)=-Inf;
 
wat = watershed(intensities);

bg=prctile(intMap(:), 5); 
th=bg+threshold; 
maskPos=intMap>th;
touchZero=unique(wat(maskPos==0));
wat2=~ismember(wat,touchZero);
wat=double(wat2).*double(wat);
%wat(maskPos==0)=0;
%wat=bwlabel(wat);
maskCells=agreagator(maskNuclei, wat, intMap, true, false);
maskCells=fillRemainingCells(maskCells, maskPos);
maskNuclei=(maskNuclei>0).*double(maskCells);




if ratio~=1
    maskNuclei=imresize(maskNuclei, originalSize,'nearest');
    maskCells=imresize(maskCells, originalSize,'nearest');
end

end

function maskCells3=fillRemainingCells(maskCells, maskPos)
maskCells=removeJonctions(maskCells);
D = bwdist(maskCells>0);
%D(maskPos==0)=-Inf;

maskCells2=uint16(watershed(D));
maskCells2(maskPos==0)=0;
maskCells2=bwlabel(maskCells2>0);
maskCells3=removeSmallIntensity(maskCells2, maskCells, 0);
maskCells3=bwlabel(maskCells3>0);
end

function [maskMito]=segmentMitochondrion(maskCells, mito, mitoThrehsold, sigma, minSize)

mitoF=imfilter(double(mito), fspecial('gaussian', (sigma)*6+1, sigma), 'symmetric', 'same'); 
mitoBG=imfilter(mitoF, fspecial('gaussian', round((sigma*4)*6+1), sigma*4),    'symmetric', 'same'); 
mitoF=max(0, mitoF-mitoBG);
sd=mitoThrehsold*std(mitoF(:));
maskMito=removeSmall((mitoF>sd)&(maskCells>0), minSize);
maskMito=maskCells.*(maskMito>0);
end

function maskCells=agreagator(maskNuclei, wat, intensities, takeMin, isDisplay)
stop=0;
maskCells=maskNuclei;

watd=imdilate(wat, ones(3));
wat(wat==0)=watd(wat==0);

wat(maskCells>0)=0;
intensities=intensities-min(min(intensities))+1;

while ~stop
    maskCellsd=imdilate(maskCells, ones(3));

    v=find(((maskCellsd>0)~=(maskCells>0))&(wat>0));
    if numel(v)==0
        stop=true;
    else
        if takeMin
            [extremaInt, indexPixExtrema]=min(intensities(v));
        else
            [extremaInt, indexPixExtrema]=max(intensities(v));
        end
        indexMin=maskCellsd(v(indexPixExtrema));
        indexStructToAdd=wat(v(indexPixExtrema));
        if indexStructToAdd==0
            error(['Error : ']);
        end
        newMask=double(wat==indexStructToAdd)*indexMin;
        maskCells=max(maskCells, newMask);
        wat(maskCells>0)=0;
        if isDisplay
            imagesc(intensities);
            drawRegions(removeJonctions(maskCells), [], true,  [], false);
            daspect([1 1 1]);
            drawnow
        end
    end
end
end

function [maskCells, maskNuclei]=segmentDonuts(maskNuclei, distMax)
D = bwdist(maskNuclei>0);
D(maskNuclei>0) = -Inf;
D(D>distMax)=-Inf;
maskCells = watershed(D);
maskCells=removeSmallIntensity(maskCells, maskNuclei>0, 0);
maskNuclei=(maskNuclei>0).*double(maskCells);
end

function im=getImageInList(imList, imName, key)

im=[];
for i=1:numel(imList),
    if strcmpi(imName{i}, key)
        im=imList{i};
        return;
    end
end
end



function mask=getMask(masks, maskNames, keyList)
if ~iscell(keyList)
    keyList={keyList};
end
lastMask=[];
for i=1:numel(keyList),
    mask=getImageInList(masks, maskNames, keyList{i})>0;
    if numel(lastMask)>0
        mask=logical(mask.*lastMask);
    end
    lastMask=mask;
end
end



function count=cellCount(masks, maskNames, keyList)
mask=getMask(masks, maskNames, keyList)>0;
maskTag=bwlabel(mask);
count=max(maskTag(:));
end

function area=getArea(masks, maskNames, keyList)
mask=getMask(masks, maskNames, keyList)>0;
area=sum(sum(mask>0));
end

function int=getIntensity(masks, maskNames, keyListMask, ims, imNamse, keyIm)
mask=getMask(masks, maskNames, keyListMask)>0;
im=double(getImageInList(ims, imNamse, keyIm));
int=mean(im(mask));
end

function int=getIntensityMed(masks, maskNames, keyListMask, ims, imNamse, keyIm)
mask=getMask(masks, maskNames, keyListMask)>0;
im=double(getImageInList(ims, imNamse, keyIm));
int=mean(im(mask));
end

function mask=detectNuclei(im)

imE=imerode(im, strel('disk', 17));
imF=double(im-imE);
%th=graythresh(imF);
%th=th*(max(imF(:))-min(imF(:)))+min(imF(:));
mask=imF>prctile(imF(:), 80);
mask=removeSmall(mask, 30)>0;
end

function [maskCells2]=segmentNuclei(im,offset, sigma, minSize, ratio)
if ~exist('ratio','var')
    ratio=1/8;
end
originalSize=size(im);

im=imresize(im, ratio);
minSize=minSize*ratio*ratio;


imF=imfilter(double(im), fspecial('gaussian', round((sigma)*6+1), sigma),    'symmetric', 'same'); 
imBG=imopen(imF, strel('disk', 25*8/ratio)); %imfilter(imF, fspecial('gaussian', round((sigma*4)*6+1), sigma*4),    'symmetric', 'same'); 
imF=max(imF-imBG, 0);

med=median(imF(:));

maskCells=(imF>offset*50);
maskCells=imfill(maskCells, 'holes');
maskCells=removeSmall(maskCells, minSize);
maskCells2=imerode(maskCells, ones(3));


if ratio~=1
    maskCells2=imresize(maskCells2, originalSize,'nearest');
end
end


function mask=detectObject(im, minSize)

backgroundTH=prctile(im(:), 5)+0.5;
mask=im>backgroundTH;
mask=imclose(mask,strel('disk', 5));
mask=gaussianFiltering(double(mask), 5)>0.5;
mask=removeSmall(mask, minSize)>0;
maskF=imfill(mask, 'holes');
maskF=maskF-mask;
maskF=removeBig(maskF, minSize)>0;
mask=mask+maskF;
mask=bwlabel(mask);
end

function [data, dataHeaders]=measureObject(im2D, mask, tagChannel)
R = regionprops(mask,double(im2D),  ...
    'MaxIntensity', 'MeanIntensity', 'Area');
RLoc=R(1);
fn=fieldnames(RLoc);
featureCount=ones(numel(fn), 1);
for i=1:numel(fn),
    featureCount(i)=numel(RLoc.(fn{i}));
end
nbDataCol=sum(featureCount)+1; % +1 for the cellId
data=[];
dataHeaders={'CellId'};
for i=1:numel(R)
    RLoc=R(i);
    dataLoc=zeros(1, nbDataCol);
    dataLoc(1)=i;
    cpt=1; % 1 for CellId
    for f=1:numel(fn),
        vals=RLoc.(fn{f});
        for fi=1:numel(vals),
            cpt=cpt+1;
            dataLoc(1, cpt)=vals(fi);
            if i==1,
                if numel(vals)==1,
                    featureName=fn{f};
                else
                    featureName=[fn{f}, num2str(fi)];
                end
                if exist('tagChannel', 'var')
                    featureName=[featureName, '_', tagChannel];
                end
                dataHeaders{1, cpt}=featureName;
            end
        end
    end
    data=cat(1, data, dataLoc);
end

end

function [data, dataHeaders]=measureObjectMask(mask)

R = regionprops(mask,...
    'Area', 'BoundingBox', 'Centroid', 'ConvexArea', 'Eccentricity', ...
    'EquivDiameter', 'EulerNumber', 'Extent', 'FilledArea', ...
    'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Perimeter', ...
    'Solidity');
RLoc=R(1);
fn=fieldnames(RLoc);
featureCount=ones(numel(fn), 1);
for i=1:numel(fn),
    featureCount(i)=numel(RLoc.(fn{i}));
end
nbDataCol=sum(featureCount)+1; % +1 for the cellId
data=[];
dataHeaders={'CellId'};
for i=1:numel(R)
    RLoc=R(i);
    dataLoc=zeros(1, nbDataCol);
    dataLoc(1)=i;
    cpt=1; % 1 for CellId
    for f=1:numel(fn),
        vals=RLoc.(fn{f});
        for fi=1:numel(vals),
            cpt=cpt+1;
            dataLoc(1, cpt)=vals(fi);
            if i==1,
                if numel(vals)==1,
                    featureName=fn{f};
                else
                    featureName=[fn{f},  num2str(fi)];
                end
                dataHeaders{1, cpt}=featureName;
            end
        end
    end
    data=cat(1, data, dataLoc);
end


end

function bwimCentroid=pointOnCentroid(bwim)
stats=regionprops(bwim, 'Centroid');
bwimCentroid=zeros(size(bwim));
for i=1:numel(stats),
    if max(isnan(stats(i).Centroid))==0
        bwimCentroid(round(stats(i).Centroid(2)), round(stats(i).Centroid(1)))=1;
    end
end
end

function maskCells=detectCellsGray(im, sigma, minSize)

imF=255-imfilter(double(im), fspecial('gaussian', round((sigma)*6+1), sigma),    'symmetric', 'same'); 
backgroundTH1=prctile(imF(:), 5)+0.13;
maskCells=imF>backgroundTH1;

intesnities=-imF;
intesnities=intesnities-max(max(intesnities));
intesnities(maskCells==0)=-Inf;
maskCells2=uint16(watershed(intesnities));
maskCells2(maskCells==0)=0;
maskCells=removeSmall(maskCells2, minSize);


end

