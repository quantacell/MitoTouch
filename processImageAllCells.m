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

function [im, imNames,masks, maskNames, data, dataHeaders]=processImageAllCells(inputFilename, inputFilenameMask, ...
                    outputFilename)


im={};
im1=double(imread(inputFilename));


imNames={};
for i=1:size(im1,3)
    im{i}=im1(:,:,i);
    imNames{i}=[num2str(i)];
end
sample=double(imread(inputFilenameMask));

blue=getImageInList(im, imNames, '3');
masks{1}=segmentNuclei(blue,10, 0.3, 8);
masks{1}=bwlabel((masks{1}).*(sample>0));
maskNames{1}='Nuclei';
[masks{2}, masks{1}]=segmentDonuts(masks{1}, 7);
maskNames{2}='CellsMask';


[data, dataHeaders]=measureMasks(masks, maskNames, im, imNames);;
T=table();

for t=1:numel(dataHeaders),
    name=dataHeaders{t};
    if (name(1)=='1') || (name(1)=='2') || (name(1)=='3') 
        name=['ch', name];
    end
    T.(name)=data(:,t);
end
csvFilename=[outputFilename, '.csv'];
writetable(T, csvFilename);

for m=1:numel(maskNames)
    imwrite(uint16(masks{m}), [outputFilename, '_mask', maskNames{m},'.tif']);
end

end

function analyseThresholds(dataHeaders, data, im1, im, imNames, masks, maskNames, sample, outputFilename)
[path, file, ext]=fileparts(outputFilename);
path=strrep(path, 'results\', 'resultsThreshold\');
path=strrep(path, 'mask', '');
if ~exist(path, 'dir'),
    mkdir(path);
end
outputFilename=fullfile(path, file);
im1=incrustateRegions(im1, sample, [0.8,0.8,0.8]);
im1=incrustateRegions(im1, getImageInList(im, imNames, 'epiderm'), [1,1,0], 2);

redThreshold=[];
greenThreshold=[];
cellCountEpi=[];
cellCountDerm=[];
cellCountEpiGreenPos=[];
cellCountEpiRedNeg=[];
cellCountEpiRedNegGreenPos=[];
cellCountDermGreenPos=[];
cellCountDermRedNeg=[];
cellCountDermRedNegGreenPos=[];
areaSample=sum(sum(sample>0));
areaSampleTab=[];

areaEpi=sum(sum(getImageInList(im, imNames, 'epiderm')>0));
areaEpiTab=[];


for redTh=70:10:120
    for greenTh=100:10:150
        areaSampleTab=cat(1, areaSampleTab, areaSample);
        areaEpiTab=cat(1, areaEpiTab, areaEpi);
        
        redThreshold=cat(1, redThreshold, redTh);
        greenThreshold=cat(1, greenThreshold, greenTh);


        selEpi=data(:,findInCells(dataHeaders, 'epiderm_mean'))>0.5;
        selDerm=data(:,findInCells(dataHeaders, 'epiderm_mean'))<=0.5;
        selRedNeg=data(:,findInCells(dataHeaders, '3_mean'))<redTh;
        selGreenPos=data(:,findInCells(dataHeaders, '2_mean'))>greenTh;
        %selRedGreenPos=selRedPos&selGreenPos;

        cellCountEpi=cat(1, cellCountEpi, sum(selEpi));
        cellCountDerm=cat(1, cellCountDerm, sum(selDerm));
        cellCountEpiGreenPos=cat(1, cellCountEpiGreenPos, sum(selEpi & selGreenPos));
        cellCountEpiRedNeg=cat(1, cellCountEpiRedNeg, sum(selEpi & selRedNeg));
        cellCountEpiRedNegGreenPos=cat(1, cellCountEpiRedNegGreenPos, sum(selEpi & selRedNeg & selGreenPos));

        cellCountDermGreenPos=cat(1, cellCountDermGreenPos, sum(selDerm & selGreenPos));
        cellCountDermRedNeg=cat(1, cellCountDermRedNeg, sum(selDerm & selRedNeg));
        cellCountDermRedNegGreenPos=cat(1, cellCountDermRedNegGreenPos, sum(selDerm & selRedNeg & selGreenPos));

        im2=incrustateRegions(im1, masks{2}, [0,1,1], 1, selGreenPos & ~selRedNeg);
        im2=incrustateRegions(im2, masks{2}, [1,0,1], 1, selGreenPos & selRedNeg);
        im2=incrustateRegions(im2, masks{2}, [0.5,0.5,0.5], 1, ~selGreenPos);
        
        image(im2)
        imwrite(im2, [outputFilename, '_cellMontage_thR', num2str(redTh),'_thG', num2str(greenTh), '.jpg']);
    end
end


for redTh2=0:10:50
    for greenTh2=50:10:100
        areaSampleTab=cat(1, areaSampleTab, areaSample);
        areaEpiTab=cat(1, areaEpiTab, areaEpi);
        
        redTh=median(data(:,findInCells(dataHeaders, '3_mean')))+redTh2;
        greenTh=median(data(:,findInCells(dataHeaders, '2_mean')))+greenTh2;
        redThreshold=cat(1, redThreshold, redTh);
        greenThreshold=cat(1, greenThreshold, greenTh);


        selEpi=data(:,findInCells(dataHeaders, 'epiderm_mean'))>0.5;
        selDerm=data(:,findInCells(dataHeaders, 'epiderm_mean'))<=0.5;
        selRedNeg=data(:,findInCells(dataHeaders, '3_mean'))<redTh;
        selGreenPos=data(:,findInCells(dataHeaders, '2_mean'))>greenTh;
        %selRedGreenPos=selRedPos&selGreenPos;

        cellCountEpi=cat(1, cellCountEpi, sum(selEpi));
        cellCountDerm=cat(1, cellCountDerm, sum(selDerm));
        cellCountEpiGreenPos=cat(1, cellCountEpiGreenPos, sum(selEpi & selGreenPos));
        cellCountEpiRedNeg=cat(1, cellCountEpiRedNeg, sum(selEpi & selRedNeg));
        cellCountEpiRedNegGreenPos=cat(1, cellCountEpiRedNegGreenPos, sum(selEpi & selRedNeg & selGreenPos));

        cellCountDermGreenPos=cat(1, cellCountDermGreenPos, sum(selDerm & selGreenPos));
        cellCountDermRedNeg=cat(1, cellCountDermRedNeg, sum(selDerm & selRedNeg));
        cellCountDermRedNegGreenPos=cat(1, cellCountDermRedNegGreenPos, sum(selDerm & selRedNeg & selGreenPos));

        im2=incrustateRegions(im1, masks{2}, [0,1,1], 1, selGreenPos & ~selRedNeg);
        im2=incrustateRegions(im2, masks{2}, [1,0,1], 1, selGreenPos & selRedNeg);
        im2=incrustateRegions(im2, masks{2}, [0.5,0.5,0.5], 1, ~selGreenPos);
        
        image(im2)
        imwrite(im2, [outputFilename, '_cellMontage_thRMed+', num2str(redTh2),'_thGMed+', num2str(greenTh2), '.jpg']);
    end
end

tab=table(areaSampleTab, areaEpiTab, redThreshold, greenThreshold, cellCountEpi, cellCountDerm, cellCountEpiGreenPos, cellCountEpiRedNeg, cellCountEpiRedNegGreenPos, ...
    cellCountDermGreenPos, cellCountDermRedNeg, cellCountDermRedNegGreenPos);

writetable(tab, [outputFilename, '_thresholdAnalysis.xlsx']);
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

function [maskCells2]=segmentNuclei(im,offset, sigma, minSize)



imF=imfilter(double(im), fspecial('gaussian', round((sigma)*6+1), sigma),    'symmetric', 'same'); 
imBG=imfilter(imF, fspecial('gaussian', round((sigma*4)*6+1), sigma*4),    'symmetric', 'same'); 
imF=imF-imBG;
%imFe=imerode(imF, strel('disk', 15));
%level=median(imF(:))+offset;
med=median(imF(:));
%offset=max(offset, med+(prctile(imF(:), 99)-med)/5);
%offset=med+(prctile(imF(:), 99.9)-med)/5;


maskCells=(imF>offset);
maskCells=removeSmall(maskCells, minSize);
maskCells=imdilate(maskCells, ones(3));
intesnities=-imF;
%intesnities=imfilter(double(intesnities), fspecial('gaussian', (sigma)*6+1, sigma),    'symmetric', 'same'); 
intesnities=intesnities-max(max(intesnities));
intesnities(maskCells==0)=-Inf;
maskCells2=uint16(watershed(intesnities));
maskCells2(maskCells==0)=0;
maskCells2=removeSmall(maskCells2, minSize);
%maskCells2=keepRound(maskCells2, 0.7);

%imagesc(im);
%drawRegions(maskCells2, [], false);
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
%intesnities=imfilter(double(intesnities), fspecial('gaussian', (sigma)*6+1, sigma),    'symmetric', 'same'); 
intesnities=intesnities-max(max(intesnities));
intesnities(maskCells==0)=-Inf;
%intesnities(mask==0)=-Inf;
maskCells2=uint16(watershed(intesnities));
maskCells2(maskCells==0)=0;
%maskCells2(mask==0)=0;
maskCells=removeSmall(maskCells2, minSize);


end

