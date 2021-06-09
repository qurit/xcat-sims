function scale_factor = scaling_factor(BqDim , d_x , d_y , d_z , t1 , t2 , lambda , calib_fact , IMAGE_DECAYED)

% scale_factor = scaling_factor(BqDim , d_x , d_y , d_z , t1 , t2 , lambda , calib_fact)
% This function returns the image corrected for the following factors:
% Bq oder of magnitude, voxel size , time , decay and calibration factor.
% volume of a voxel is needed, and calculated via d_x, d_y.
% t1, t2 and lambda are needed for calucating the decay factor.
% calibration factor can be evaluated from calib_fact.m
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
% --> Copyright (C) 2013-2018  Saeed Ashrafinia, Johns Hopkins University
% 
%    This package is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
% -------------------------------------------------------------------------

% Correct for Bq unit
if BqDim == 'K'
    ratio = 1e3;
elseif BqDim == 'M'
    ratio = 1e6;
else
    ratio = 1;
end

% find voxel volume
vol = d_x * d_y * d_z;

% find time
delta_t = t2 - t1; 

% decay
if  IMAGE_DECAYED == 0
    decay_factor = exp(-lambda * t1) * (1 / (lambda * delta_t)) * (1 - exp(-lambda * delta_t));
else
    decay_factor = 1;
end

% scaling factor
scale_factor = ratio * vol * (60 * delta_t) * decay_factor * calib_fact;
