function scaled_out = calib_factor(xdim , ydim , VCT_sens , theta_m , normalization , BIN_NUM , TotalCounts)
% function scaled_out = calib_fact(xdim , ydim , zdim , VCT_sens , theta)

% This function calculates the scaling factor for the PET imaging system.
% The image has the dimension of (xdim,ydim), and only "VCT_sens" percent of
% the counts are successfully being detected. In order to find the correct
% scaling, we generate an image of rods (in this case, just a slice of a
% rod which is a centered voxel(s)), then forward project it. The scaling
% factor is evaluated from the ratio of the number of counts in the image
% and in sinogram with respect to the VCT_sens.
% -------------------------------------------------------------------------
% % AUTHOR: 
% - Saeed Ashrafinia
% -------------------------------------------------------------------------
% **** If you use this code in a study, please cite the following paper ***
% S. Ashrafinia, H. Mohy-ud-din, N. Karakatsanis, A. Jha, M. Casey, D.
% Kadrmas and A. Rahmim, “Generalized PSF modeling for optimized
% quantitative-task performance”, Phys. Med. Biol., vol. 62, pp. 5149-5179,
% 2017.   
% 
% -------------------------------------------------------------------------
% STATEMENT:
% This file is part of PET Recon Package by Saeed Ashrafinia, Rahmimlab.com
% --> Copyright (C) 2013-2017  Saeed Ashrafinia, Johns Hopkins University
% 
%    This package is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
% -------------------------------------------------------------------------

% Generating the rod
PSF_width = 0;
rod_img = zeros(xdim,ydim);
rod_img(ceil(xdim/2):ceil(xdim/2)+PSF_width , ceil(ydim/2):ceil(ydim/2)+PSF_width) = TotalCounts;

% Calculate scaling factor for each z direction
% rod_sino = radon(rod_img , theta , BIN_NUM);
rod_sino =zeros(size(normalization));
for m=1:size(theta_m ,2);
    rod_sino (:,:,m) = radon(rod_img , theta_m(:,m) , BIN_NUM);
end


rod_sino = rod_sino .* normalization;
scaled_out = sum(rod_img(:)) * VCT_sens / sum(rod_sino(:));



% **** 3D ****
% rod_img = zeros(xdim,ydim,zdim);
% rod_img(ceil(xdim):ceil(xdim)+1 , ceil(ydim):ceil(ydim)+1 , ceil(zdim):ceil(zdim)+1) = 1e9;

% % Calculate scaling factor for each z direction
% scaled_out = zeros(1,zdim);
% for z = 1 : zdim
%     rod_sino = radon(rod_img , theta);
%     scaled_out(z) = sum(rod_img(:)) * percnt / sum(rod_sino(:));
% end