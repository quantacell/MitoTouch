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

function im=gaussianFiltering(im, sigma, scale)
    if ~exist('scale', 'var')
        scale=round(3*sigma);
    end
    psfEs=fspecial('gaussian', [1, 2*scale+1], sigma)';
    for z=1:size(im,3)
        imLoc=imfilter(im(:,:,z)', psfEs,    'symmetric', 'same'); 
        im(:,:,z)=imfilter(imLoc', psfEs,    'symmetric', 'same'); 
    end
end

