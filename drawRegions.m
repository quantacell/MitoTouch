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

function B=drawRegions(maskCells, handles, displayTag, color, isHole)
if ~exist('displayTag', 'var')
    displayTag=true;
end

if ~exist('handles', 'var')
    handles=[];
end

if ~exist('color', 'var')
    color=[];
end
if displayTag
    maskCells=imerode(maskCells, [0,1,0; 1, 1, 1; 0, 1, 0]);
end
%B=bwboundaries(maskCells>0,4, 'noholes');
if ~exist('isHole', 'var') || ~isHole
    B=bwboundaries(maskCells>0,4, 'noholes');
else
    B=bwboundaries(maskCells>0,4);
end
if numel(handles)==0
   hold('on');
else
    hold(handles, 'on');
end
for k = 1:length(B)
    boundary = B{k};
    if displayTag
        index=maskCells(boundary(1,1), boundary(1,2));
    else
        index=k;
    end
    if numel(color)==0,
        col=getColor(index);
    else
        col=color;
    end
    if numel(handles)==0
        plot(boundary(:,2), boundary(:,1), 'w', 'Color', col)
    else
        plot(handles, boundary(:,2), boundary(:,1), 'w', 'Color', col)
    end
    if displayTag
        y=min(boundary(:,1))-20;
        x=mean(boundary(:,2));
        if numel(handles)==0
            text(x, y, num2str(index), 'Color', col);
        else
            text(x, y, num2str(index), 'Color', col, 'Parent', handles.axes1);
        end
    end
end
if numel(handles)==0
   hold('off');
else
    hold(handles, 'off');
end

end

function color=getColor(index);
%colors=[0, 0, 0.5625; 0.4375, 1, 0.5625; 0.6875, 0, 0; 0, 0, 0.875; 0, 0.1875, 1; 0, 0.5, 1; 0, 0.8125, 1; 0.125, 1, 0.875; 0.75, 1, 0.25; 1, 0.9375, 0; 1, 0.625, 0; 1, 0.3125, 0; 1, 0, 0];
colors=[0, 0, 128; 0, 255, 255; 34, 139, 34; 127, 255, 0; 255, 215, 0; 255, 127, 80; 252, 15, 192;
    181, 126, 220; 128, 0, 0; 132, 49, 121]/255;
color=colors(mod(index-1, size(colors, 1))+1, :);

end