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

function [ maskOut ] = removeAllExceptBrigthest( mask,im, N )

mask=bwlabel(mask>0);
R = regionprops(mask,im, 'MeanIntensity');
[~, ind]=sort([R.MeanIntensity], 'descend');
if N>numel(ind)
    maskOut=mask;
else
    maskOut = bwlabel(ismember(mask,ind(1:N)));
end

end

