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

function [T, cptImpage]=processAllImages(mainFolder, ouputFolder, T, cptImpage, nucleiThreshold, cellThreshold, mitoThreshold, erodeDist, filterBorder, isCellBased, isParallel)

if numel(mainFolder)==0
    mainFolder=uigetdir();
    if numel(mainFolder)<=1
        return;
    end
    
    
    ouputFolder=fullfile(mainFolder, 'results');
    if ~exist(ouputFolder, 'dir')
        mkdir(ouputFolder);
    end
    tic;
    listOfImages={};
    listOfImages=scanFolder(mainFolder, listOfImages);
    
    
    cptImpage=numel(listOfImages);
    
    if cptImpage==0
        msgbox(['Unable to find any image in folder ', mainFolder],'modal');
        T=table();
        return;
    end
    
    try
        T=processAllImages2(listOfImages, ouputFolder, nucleiThreshold, cellThreshold, mitoThreshold, erodeDist, filterBorder, isCellBased, isParallel);
        csvFilename=fullfile(ouputFolder, 'final.csv');
        writetable(T, csvFilename);

        t=toc;
        writeParameters(listOfImages, ouputFolder, nucleiThreshold, cellThreshold, mitoThreshold, erodeDist, filterBorder, isCellBased, isParallel, t);
    catch ex
        msgException(ex);
        if ~exist('T', 'var')
            T=table;
        end
        if ~exist('cptImpage', 'var')
            cptImpage=0;
        end
    end
    return;
end

list=dir(mainFolder);
try
for i=1:numel(list)
    if list(i).isdir
        if list(i).name(1)~='.'
            [T, cptImpage]=processAllImages(fullfile(mainFolder,list(i).name), ouputFolder, T, cptImpage, nucleiThreshold, cellThreshold, mitoThreshold, erodeDist, filterBorder, isCellBased, isParallel);
        end
    end
end
catch ex
    msgException(ex);
    if ~exist('T', 'var')
        T=table;
    end
    if ~exist('cptImpage', 'var')
        cptImpage=0;
    end
end

end

function msgException(ex)
msgbox(['Error: ', ex.message], 'modal');
for i=1:numel(ex.cause)
    msgbox(['Cause: ', ex.cause{i}], 'modal');
end
end

function writeParameters(listOfImages, ouputFolder, nucleiThreshold, cellThreshold, mitoThreshold, erodeDist, filterBorder, isCellBased, isParallel, processingTime)
t = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH-mm-ss'));
t=strrep(t, ':','-');
filename=['parameters ', t, '.txt'];
str={};
str{numel(str)+1}='List of parameters';
str{numel(str)+1}=['Date: ', t];
str{numel(str)+1}=['ouputFolder: ', ouputFolder];
str{numel(str)+1}=['Filename: ', filename];
str{numel(str)+1}='List Of Images:';
for i=1:numel(listOfImages)
    str{numel(str)+1}=listOfImages{i};
end
str{numel(str)+1}=['nucleiThreshold: ', num2str(nucleiThreshold)];
str{numel(str)+1}=['cellThreshold: ', num2str(cellThreshold)];
str{numel(str)+1}=['mitoThreshold: ', num2str(mitoThreshold)];
str{numel(str)+1}=['erodeDist: ', num2str(erodeDist)];
str{numel(str)+1}=['filterBorder: ', num2str(filterBorder)];
str{numel(str)+1}=['isCellBased: ', num2str(isCellBased)];
str{numel(str)+1}=['isParallel: ', num2str(isParallel)];
str{numel(str)+1}=['processingTime (sec): ', num2str(processingTime)];

fid=fopen(fullfile(ouputFolder, filename), 'w');
for i=1:numel(str)
    fprintf(fid, '%s\n', str{i});
end
fclose(fid);
end

function T=processAllImages2(listOfImages, ouputFolder, nucleiThreshold, cellThreshold, mitoThreshold, erodeDist, filterBorder, isCellBased, isParallel)

TList={};
filter='MERGE.tif';
if isParallel
    parfor i=1:numel(listOfImages)
        [well, field, well_field, subFolder1, subFolder2]=processFilename(listOfImages{i}, ouputFolder, filter);
        ouputFolderLoc=fullfile(ouputFolder,[well,'_', field]);
        [~, ~, Tloc]=processImage(listOfImages{i}, ouputFolderLoc,nucleiThreshold, cellThreshold, mitoThreshold, erodeDist, filterBorder, isCellBased);
        Tloc.Well=repmat({well},size(Tloc, 1),1);
        Tloc.Field=repmat({field},size(Tloc, 1),1);
        Tloc.WellField=repmat({well_field},size(Tloc, 1),1);
        Tloc.SubFolder1=repmat({subFolder1},size(Tloc, 1),1);
        Tloc.SubFolder2=repmat({subFolder2},size(Tloc, 1),1);

        TList{i}=Tloc;
    end
else
    for i=1:numel(listOfImages)
        [well, field, well_field, subFolder1, subFolder2]=processFilename(listOfImages{i}, ouputFolder, filter);
        ouputFolderLoc=fullfile(ouputFolder,[well,'_', field]);
        [~, ~, Tloc]=processImage(listOfImages{i}, ouputFolderLoc,nucleiThreshold, cellThreshold, mitoThreshold, erodeDist, filterBorder, isCellBased);
        Tloc.Well=repmat({well},size(Tloc, 1),1);
        Tloc.Field=repmat({field},size(Tloc, 1),1);
        Tloc.WellField=repmat({well_field},size(Tloc, 1),1);
        Tloc.SubFolder1=repmat({subFolder1},size(Tloc, 1),1);
        Tloc.SubFolder2=repmat({subFolder2},size(Tloc, 1),1);

        TList{i}=Tloc;
    end
end
T=[];
for i=1:numel(TList)
    if numel(TList{i})>0
        if numel(T)==0
            T=TList{i};
        else
            try
                T=cat(1, T, TList{i});
            catch
                
            end
        end
    end
end
end



function listOfImages=scanFolder(ouputFolder, listOfImages)

list=dir(ouputFolder);
for i=1:numel(list)
    if list(i).isdir
        if list(i).name(1)~='.'
            newFolder=fullfile(ouputFolder, list(i).name);
            listOfImages=scanFolder(newFolder, listOfImages);
        end
    end
end

filter='ch1sk1fk1fl1.tiff';


list=dir(fullfile(ouputFolder,['*', filter]));
for i=1:numel(list)

    file=fullfile(ouputFolder,list(i).name);
    listOfImages=cat(1, listOfImages, file);
end

end

function [well, field, well_field, subFolder1, subFolder2]=processFilename(file, outputFolder, filter)
di=find(file(1:numel(outputFolder))~=outputFolder);
    last=strtrim(file(di(1):end-numel(filter)));
    if last(end)=='-'
        last=strtrim(last(1:end-1));
    end
    well=lower(last);
    field=well;

    campPos=strfind(well, '-');
    well=well(1:campPos(1)-1);
    field=field(campPos(1)+1:end);

    wellSplit=strsplit(well, '\');
    if numel(wellSplit)>1
        subFolder1=wellSplit{end-1};
    else
        subFolder1='missing';
    end

    if numel(wellSplit)>2
        subFolder2=wellSplit{end-2};
    else
        subFolder2='missing';
    end
    
    well=wellSplit{end};
    well_field=[well, '_', field];
end

