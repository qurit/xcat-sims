function [sensit_image_all] = Gen_SensSino(Kernel_Set , atten , norm_image , NUMVAR , SUBSETS , RECONST_RM , theta_m , xdim , ydim , filter_tr , iRadonInterp)
% This function creates the sensitivity sinogram for reconstruction
% simulatiion. It applies the RM kernel and has the option to
% reconstruct with RM.
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



[binsdim , thdim_m , ~] = size(atten);
sensit_image_all = zeros(xdim , ydim , SUBSETS , NUMVAR);

for count_var = 1:NUMVAR
    % Load kernel
    convolFINAL_tr = transpose(squeeze(Kernel_Set(:,:,count_var)));
    
    % Creating sensitivity image
    % Projection-space blurring (transposed)
    for m=1:SUBSETS
        % Projection-space transpose operation
        proj=atten(:,:,m).*norm_image(:,:,m);
        if RECONST_RM
            proj2=zeros(binsdim,thdim_m);
            for th=1:thdim_m
                for bin=1:binsdim
                    proj2(:,th)=proj2(:,th)+proj(bin,th)*convolFINAL_tr(bin,:)';
                end
            end
        else
            proj2=proj;
        end
        
        % Geometric projection
        aaa = iradon(proj2,theta_m(:,m),iRadonInterp,'none');
        aaa2=aaa(2:xdim+1,2:ydim+1);
        
        % Image-space blurring (transposed)
        if RECONST_RM
            sensit_image_all(:,:,m,count_var)=conv2(aaa2,filter_tr,'same');
        else
            sensit_image_all(:,:,m,count_var)=aaa2;
        end
    end
end