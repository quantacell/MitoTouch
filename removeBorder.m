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


function maskOut=removeBorder(mask) 
if islogical(mask)
    mask=bwlabel(mask);
end
R = regionprops(mask,'BoundingBox');
ind=[];
for i=1:numel(R),
    if (ceil(R(i).BoundingBox(1))>1 ) && (ceil(R(i).BoundingBox(2))>1 ) && ...
            (ceil(R(i).BoundingBox(1))+R(i).BoundingBox(3)-1<size(mask, 2)) && ...
            (ceil(R(i).BoundingBox(2))+R(i).BoundingBox(4)-1<size(mask, 1))
        ind=[ind, i];
    end
end
maskOut = bwlabel(ismember(mask,ind));
end