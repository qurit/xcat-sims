% Please read all the following information carefully. 
% 
% This framework performs a 2D simulation of PET data with and without
% noise, using analytical modeling of various image degradation effects.   
% It also performs image reconstruction, and has an option for point-spread
% function (PSF) modeling, also knows as PSF modeling, as described below.
% (Please note RM in this code refers to resolution- or PSF modeling.)
% 
% This file requires an input true image "ImgTrue", and in case of modeling
% attenuation, it requires an attenumation map "Mu_Map.mat", and in case of
% modeling detector normalization, it requires a normalization map
% "norm_map". Please put these three files in the same directory as this
% main code. There are already 2 normalization maps included in this
% package for 128x128 and 256x256 images that are for GE Discovery RX. 
% 
% The true image should be called from a function or loaded from a 
% variable at the following section below in this code: 
% "%% Load Reference True Image" in line 99 or 
% "%% Generate True Image" in line 243. Some examples are included in this
% line, such as a NEMA phantom generator. 
% 
% The size of Mu_Map should be in correspondence with the reconstruction 
% image dimension. The size of normalization_map should be the same size as
% the image sinogram (no subsets). Refer to number of angles and bins in 
% "Gen_PSF_Kernels" to check that. Both mu_map and norm_map should be
% located and are loaded from the current directory.
% 
% The uptake value of the tumor and its background tissue may be required
% for quantitation analysis, and are set in the variable "ROIs".
% 
% THIS CODE RUNS "Gen_PSF_Kernels.m" TO GENERATE PSF KERNELS. 
% Scanner parameters are determined in that file.
% 
% This code has been developed for reconstructing a single tumour. For
% multiple tumors this code has to be run separately and the variable
% "TmrCount" should be set to track the tumor index.
% THIS CODE STARTS WITH A HIGHER RESOLUTION IMAGE AND SCALES THE RADON T
% DOWN TO MATCH THE OUTPUT IMAGE.
% 
% This code saves a cropped version of the reconstruction image to preserve
% RAM and hard drive space. The cropped image dimensions can be set in
% variable "BkgndBox" with respect to the image size.
% 
% The output variables depend on whether the reconstruction is
% noisy/noise-free and signal-present/signal-absent and are SIG_PRS,
% SIG_PRS_NF, SIG_ABS, and SIG_ABS_NF, respectively. 
% Each of these variables, the program saves the reconstructed images of
% every iteration, every PSF kernel, and every noise realization. The
% noise-free variables have one dimension less (just 1 realization).
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
    
%%

clear;  warning('Clearing all variables. Turn it off if needed.'); 

tic

%% Defining Variables
SLICES = 301;    % RF. Number of slices in image.
CYCLES = 3;    % RF. NUMBER OF NOISE REALIZATIONS

for c = 1:CYCLES
    
PSF_Kernel  = 0;        % 0 for no PSF modeling, 1 for true PSF modeling, 2 for no PSF + true PSF modeling, 3 for generalized PSF modeling (no PSF + true PSF + under- and over-estimated PSFs). Generalized PSF are defined in Gen_PSF_Kernels.m.
NOISE_REALZ = 1;        % Number of noise realizations. Please modify "CYCLES" variable instead and leave this as one. % RF.
ITERATIONS  = 2;
SUBSETS     = 24;       % Number of OSEM Subsets
xdim        = 256;    
input_xdim = 512;       % Size of the input image (can be the same, x2, x4 or x8 of xdim), max 1024

IMG_ABS_PRS = 0;        % 0: run for signal absent; 1: present; 2: absent and present
RECON_NF_NOISY = 1;     % 0 for Noise-Free, 1 for Noisy, 2 for Noisy and Noise-free.

RECONST_RM = 1;         % Using RM in reconstruction
SIMULATE_RM = 1;        % Using RM in forward projection
IMAGE_DECAYED = 0;      % Normaly 0; set to 1 if image is not decay corrected.
HIGH_RES_TRUE = 1;      % Use a higher resolution image as the true input and then scale down it projection
LOAD_ATTENUATION = 0;   % Load attenuation map or consider a matrix of zeros
CORRECT_IMG_UPTKE = 0;  % Depends if the initial image is in the units of activity of concentration. If it's already concent., that usually has to be, there's no need to do it.
LOAD_NORMALIZATION = 1; % Load normalization image or consider a matrix of ones

SAVE        = 1;        % If save the final result      
DIRECTORY = char(pwd);  % Saving Directory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder=char('../input/')      % FOLDER TO RETRIEVE XCAT IMAGES FROM
file=char('snm1_act_1')
%file=char('example')                % NAME OF FILE TO SIMULATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate Partial PSF kernels
% ***** To select kernels refer to "Gen_PSF_Kernels" *****
Gen_PSF_Kernels         % Generate partial PSF kernels

%% Load Reference True Image
% More information on the reference image is right after beginning of the
% first for loop. 
%  Gen_Image_True;       % Load the XCAT phantom reference image (loads "ImgTrue")

name=append(folder,file)
n = strcat(name,'.bin')

fileID = fopen(n, 'r');
ImgTrue = fread(fileID, 512*512*SLICES,'single');
ImgTrue = reshape(ImgTrue, 512, 512, SLICES);
MSA_obj = []

BKGND_SUV       = 10;   % Uptake value of the background tissue (of the tumor)
TUMOR_SUV       = 24;   % Uptake value of the tumor tissue

%% Save tumor and background values
TmrCount        = 6;    % Select from 1 to 8 different tumors, Refer to Gen_Image_True for indices. This variable is only used for naming the output image.   

%% Continue Initializing
ScanDuration    = 3;   % Scan duration in minutes, *** If changed here, change in FDG_TAC_all.m too"
VCT_sensitivity = 7300/1e6; % cps/MBq = cps/Mdps MBq->KBq; 
                        %The number is 1700 for 2D, and 7300 for 3D for GE Discovery RX scanner.

FOV = 384; % RF added this variable.

ydim            = xdim;
d_z             = 1.5; % originally 3.27mm. Changed by RF. 
voxel_size      = FOV/xdim; % Voxel size
TotalCounts     = 1e9;
rng('shuffle', 'twister');
iRadonInterp    = 'linear';

% If image requires correction. If it's activity, this changes it to concentration.
% For a new image with different resolution, this should be corrected accordingly.
if CORRECT_IMG_UPTKE                        
    image_correction = 1 ./ (1000 * voxel_size/10 * voxel_size/10 * d_z/10);
else
    image_correction = 1;                   %#ok<*UNRCH>
end
ROIs            = [BKGND_SUV TUMOR_SUV];  % [background_tissue_uptake Tumour_tissue_uptake]  

% offset correction for higher input resolution down sampling. Should be -6
% for HiResXdimFactor = 4.
% Apparently the best value for HiResXdimFactor=8 is -14 and -13 (the
% latter is slightly better), a factor of 2 is not a good choice.
if exist('ImgTrue' , 'var') %|| size(ImgTrue , 1) ~= in_xdim
    input_xdim                 = size(ImgTrue , 1);
    disp('The size of ref image is now determined from itself; ignoring predefined "in_xdim".');
end
% How much scaling up the true high resolution image
HiResScale         = input_xdim / xdim;

if ~HIGH_RES_TRUE
    if HiResScale ~= 1
        error('Fix HiRes sizes!');
    end
    HiResScale = 1;
end
if xdim ==128
    try
        if HiResScale  == 4
            OffsetCorr      = -6;
        elseif HiResScale == 8
            OffsetCorr      = -13;
        end
        
    catch
        error('Sorry; The offset correction for your desired input image resolution has not been previously set. Please use the above few lines to define it. You may do some try and errors to find the best offset value.');
    end
    BkgndBox            = [41 91 32 95]; % [x1 x2 y1 y2]
elseif xdim == 256
    try
        if HiResScale  == 4
            OffsetCorr      = -8;           % -8 works the best
        elseif HiResScale == 8
            OffsetCorr      = -13;
        end
    catch
        error('Sorry; The offset correction for your desired input image resolution has not been previously set. Please use the above few lines to define it. You may do some try and errors to find the best offset value.');
    end
    BkgndBox            = [1 256 1 256]; %[70 185 50 205]; % [79+10 184-24-15 60+30+40 195-5 ];% [x1 x2 y1 y2]
end


%% Input Attenuation Map and Normalization Map
% Input attenuation map. Refer to Mid_simul_code_KM for saving this var.
% For different slice of the phantom and different phantom size, this
% should change.
if LOAD_ATTENUATION
    if xdim     == 128
        load Mu_Map_128.mat
    elseif xdim ==256
        load Mu_Map_256.mat
    end
else
    mu_map = zeros(xdim , ydim);
end

% Input Normalization
if LOAD_NORMALIZATION
    if xdim == 128
        load 'norm_map_128' norm_original
    elseif xdim == 256
        load 'norm_map_256' norm_original
    end
    % Check Mid_simul_code_KM for how to save this file.
    % We're correcting for the number of bins in the normalization image.
    % The normalization sinogram contains 323 bins. Our image is 128*128, with a width of 3.47 and 3.47.
    % Therefore, one side of the image would be 128*3.47 = 444.16, and the diagonal would be that number times sqrt(2) = 628.1371.
    % Now, knowing an average bin size after arc correction is 2.249, the number of bins would be 628.1371/2.249 = 279.2962.
    % Let's take it 279. then we'll have to cut the normalization map as following:
else
    norm_original = ones(binsdim,thdim);
end
if numel(norm_original) ~= prod([binsdim,thdim])
    norm_original = imresize(norm_original,[binsdim thdim]);
end
disp(' ');
disp(['Phantom, masks and maps are loaded @ ' , num2str(toc,5) , ' sec.']);
disp('Make sure to verify the total counts (~0.5e+5) and the attenuation sinogram (~1/20).');
% disp(' ');


%% Positron range
[Y2, X2] = meshgrid(1:7,1:7);
filter = exp( -0.56*sqrt((X2-4).^2+(Y2-4).^2)*bin_size);     % coefficient should be 0.56 for Rb-82
filter = filter/sum(filter(:)); %#ok<NASGU> %disp(num2str(filter));
filter = 1;             % No positron range modeling:
filter_tr = filter;     % because filter is symmetric; it's transpose (vectorized format) is equivalent


switch IMG_ABS_PRS
    case 0
        ABS_PRS = 1;        % Reconstruct Signal Absent
    case 1
        ABS_PRS = 2;        % Reconstruct Signal Present
    case 2
        ABS_PRS = 1 : 2;    % Reconstruct Signal Absent and Present
end


for i = 1:SLICES

    disp(i);
    
for img_prs_abs = ABS_PRS      % 1: Signal Absent, 2: Signal present
       
    %% Generate reference truth image for signal abs and signal present
    % Note: if importing from an XCAT phantom, the file is usually a 3D
    % matrix of activity. Its first slice is the background image
    % containing all of the tissues. Each of the next slices represent only
    % a "tumor tissue". So for tumor_absent, we just pass the first slice.
    % For tumor_present, we have to add one of the slices to the first
    % slice. We also have to be carefule since the tumor activity is being
    % "added" to its background; not replaced.
    % Here we also make some simple examples of an image of all zeros and a
    % small square in the middle.
    
    %% Generate True Image
    
    if img_prs_abs      ==  1
        image_true      = squeeze(ImgTrue(:,:,i)); % SELECT FOR IMAGE TRUE: WHOLE IMAGE
%        image_true      = zeros(size(squeeze(ImgTrue(:,:,1))));% SELECT FOR IMAGE TRUE: SINGLE TUMOR
%        image_true      = zeros(in_xdim);% SELECT FOR IMAGE TRUE: SINGLE TUMOR
        %image_true = GenNEMA(xdim , voxel_size/HiResScale , [ROIs(1) ,  repmat(ROIs(1),1,6)]);
        T1 = 'Sig Absent';
    elseif img_prs_abs  ==  2
        image_true      = squeeze(ImgTrue(:,:,1)+ImgTrue(:,:,TmrCount+1)); %      WHOLE IMAGE
%         image_true      = squeeze(ImgTrue(:,:,TmrCount+1));  %                    SINGLE TUMOR of Whole image 
%         image_true      = ones(size(squeeze(ImgTrue(:,:,1)))); %                    All ones, values from ImgTrue
%         image_true      = ones(in_xdim); %                                          All ones
%         image_true(in_xdim/2-2*HiResScale:in_xdim/2+2*HiResScale , in_xdim/2-2*HiResScale              :              in_xdim/2+2*HiResScale) = sum(ROIs(:)); % SINGLE TUMOR @ center
%         image_true(in_xdim/2-2*HiResScale:in_xdim/2+2*HiResScale , in_xdim/2-2*HiResScale+50*HiResScale:+50*HiResScale+in_xdim/2+2*HiResScale) = sum(ROIs(:)); % SINGLE TUMOR @ off center
%         image_true(in_xdim/2-2*HiResScale:in_xdim/2+2*HiResScale , in_xdim/2-2*HiResScale-50*HiResScale:-50*HiResScale+in_xdim/2+2*HiResScale) = sum(ROIs(:)); % SINGLE TUMOR @ off center
%         image_true(in_xdim/2-2*HiResScale+50*HiResScale:+50*HiResScale+in_xdim/2+2*HiResScale , in_xdim/2-2*HiResScale:in_xdim/2+2*HiResScale) = sum(ROIs(:)); % SINGLE TUMOR @ off center
%         image_true(in_xdim/2-2*HiResScale-50*HiResScale:-50*HiResScale+in_xdim/2+2*HiResScale , in_xdim/2-2*HiResScale:in_xdim/2+2*HiResScale) = sum(ROIs(:)); % SINGLE TUMOR @ off center
%         image_true(in_xdim/2-2*HiResScale-50*HiResScale:-50*HiResScale+in_xdim/2+2*HiResScale , in_xdim/2-2*HiResScale+50*HiResScale:+50*HiResScale+in_xdim/2+2*HiResScale) = sum(ROIs(:)); % SINGLE TUMOR @ off center
%         image_true(in_xdim/2-2*HiResScale-50*HiResScale:-50*HiResScale+in_xdim/2+2*HiResScale , in_xdim/2-2*HiResScale-50*HiResScale:-50*HiResScale+in_xdim/2+2*HiResScale) = sum(ROIs(:)); % SINGLE TUMOR @ off center
%         image_true(in_xdim/2-2*HiResScale+50*HiResScale:+50*HiResScale+in_xdim/2+2*HiResScale , in_xdim/2-2*HiResScale+50*HiResScale:+50*HiResScale+in_xdim/2+2*HiResScale) = sum(ROIs(:)); % SINGLE TUMOR @ off center
%         image_true(in_xdim/2-2*HiResScale+50*HiResScale:+50*HiResScale+in_xdim/2+2*HiResScale , in_xdim/2-2*HiResScale-50*HiResScale:-50*HiResScale+in_xdim/2+2*HiResScale) = sum(ROIs(:)); % SINGLE TUMOR @ off center
%         image_true = image_true .* C; %                                             Mask with giant circle

%         image_true = GenNEMA(input_xdim , voxel_size/HiResScale , [ROIs(1) ,  repmat(ROIs(2),1,6)]);
        T1 = 'Sig Prsent';
    end
    
    %image_true = image_true.b * image_correction; % RF edited
    image_true = image_true * image_correction;
    %figure; imagesc(image_true); axis tight; colorbar; colormap(jet); title(['True image, ',int2str(size(image_true,1)),'\times',int2str(size(image_true,2))]); 

    %% Initial projection-space size and counts estimaton
    % Initialize Y_bar0 size using the information provided in "edit radon"
    % binsdim = 2*ceil(norm([xdim ydim]-floor(([xdim ydim]-1)/2)-1))+3; 
    Y_bar0  = zeros(2*ceil(norm(size(image_true)-floor((size(image_true)-1)/2)-1))+3 , angles_m , SUBSETS);
%     for m = 1:SUBSETS
%         Y_bar0(:,:,m) = radon(image_true,theta_m(:,m));
%     end
    
    atten = ones(binsdim,thdim_m,SUBSETS);
    for m = 1:SUBSETS
%         atten(:,:,m)  = exp(-radon(mu_map , theta_m(:,m))*bin_size);
        atten(:,:,m)  = exp(-radon(mu_map , theta_m(:,m))); % Because mu map saved by XCAT is also per voxel and we need per mm
    end
    
    %% Normalization map , or input the true one
    norm_image = reshape(norm_original,binsdim,thdim_m,m);
    
    
    %% Calculate sensitivity sinogram
    [sensit_image_all] = Gen_SensSino(KernelsSet , atten , norm_image , NUMVAR , SUBSETS , RECONST_RM , theta_m , xdim , ydim , filter_tr , iRadonInterp);
    disp(['Sensitivity Sinogram is generated @ ' , num2str(toc,5) , ' sec.']);
%     disp(' ');
    
    %% Calculate Calibration factor and scaling factor
    lambda = 0.0063;
    start_time = 60;
    end___time = start_time + ScanDuration;
    
    calib_fact = calib_factor(xdim , ydim , VCT_sensitivity , theta_m , norm_image , NUM_BINS , TotalCounts) ;
    scale_factor = scaling_factor('K' , voxel_size/HiResScale/10 , voxel_size/HiResScale/10 , d_z / 10 , ...
        start_time , end___time , lambda , calib_fact , IMAGE_DECAYED);
    
    %% Create Empty Matrices to save reconstructed images
    if      img_prs_abs == 1
        SIG_ABS         = single(zeros(numel(BkgndBox(1):BkgndBox(2)) , numel(BkgndBox(3):BkgndBox(4)) , ITERATIONS , NOISE_REALZ , NUMVAR));
        SIG_ABS_NF      = single(zeros(numel(BkgndBox(1):BkgndBox(2)) , numel(BkgndBox(3):BkgndBox(4)) , ITERATIONS , NUMVAR));
    elseif  img_prs_abs == 2
        SIG_PRS         = single(zeros(numel(BkgndBox(1):BkgndBox(2)) , numel(BkgndBox(3):BkgndBox(4)) , ITERATIONS , NOISE_REALZ , NUMVAR));
        SIG_PRS_NF      = single(zeros(numel(BkgndBox(1):BkgndBox(2)) , numel(BkgndBox(3):BkgndBox(4)) , ITERATIONS , NUMVAR));
    end
    
    %% Forward Projection
    % Image-space blurring
    if SIMULATE_RM
        image_true2 = conv2(image_true,filter,'same');  
    else
        image_true2 = image_true;
    end
    
    % Geometric projection
    for m = 1:SUBSETS
        [Y_bar0(:,:,m),xp] = radon(image_true2,theta_m(:,m)); 
    end
    
    
    OffsetCorr = -8;
    
    % If image_true is high res, scale the sinogram down 
    % Other factors we don't know. To know, run
    % with different C and plot using the commands in Mid_simul_codes.m to
    % see if C doesn't make any streaks.
    if HiResScale ~= 1
        Y_summed = zeros(binsdim , angles_m , SUBSETS);
        for k = 3 : binsdim  
            TmpInd = ((k-1) * HiResScale + OffsetCorr +1 ): (k * HiResScale + OffsetCorr);
            curr_bins = TmpInd(TmpInd <= size(Y_bar0,1) & TmpInd > 0); % To make sure indicies are between 1 and numbins
            Y_summed(k,:,:) = sum(Y_bar0(curr_bins , : , :) , 1);
        end
        Y_bar0 = Y_summed;
    end
    
    
    % Projection-space blurring
    if SIMULATE_RM
        Y_bar = zeros(binsdim,thdim_m,SUBSETS);
        for m = 1:SUBSETS
            for th = 1:thdim_m
                for bin = 1:binsdim
                    Y_bar(:,th,m) = Y_bar(:,th,m) + Y_bar0(bin,th,m) * KernelFull(bin,:)';
                end
            end
        end
    else
        Y_bar = Y_bar0;
    end
    % Attenuation: we assume this is done after next step; or that it's
    % smooth; to allow taking it out of forward and back-projections in EM
    % but we properly simulate here
    Y_bar = Y_bar .* atten .* norm_image;
    
    
    %% Reconstruct Noisy / Noise-Free 
    Y_bar_realization = zeros(binsdim,thdim_m,SUBSETS,NOISE_REALZ);
    
    switch RECON_NF_NOISY 
        case 0
            NNF = 1;        % Reconstruct Noise-Free
        case 1
            NNF = 2;        % Reconstruct Noisy
        case 2
            NNF = 1 : 2;    % Reconstruct Noisy and Noise-Free
    end
    for IsNoisy = NNF
        if IsNoisy == 2     % 2: add noise  1: noise-free
            %% Add noise
            Y_bar = Y_bar * scale_factor;
            disp(['The total number of counts is: ' , num2str(sum(Y_bar(:)),6)]);
            for count_rlz = 1:NOISE_REALZ
                Y_bar_realization(:,:,:,count_rlz) = poissrnd(Y_bar); % renders the data noisy
            end
            Y_bar = Y_bar / scale_factor;
            Y_bar_realization = Y_bar_realization / scale_factor/ HiResScale^2; %dividing by factor^2 to compensate finer resolution in d_x and d_y
            T2 = 'Noisy';
            tmp_NOISE_REALZ = NOISE_REALZ;
            
            
        else
            % If no noise, added, we only have single run
            Y_bar_realization(:,:,:,1) = Y_bar / HiResScale^2; %/ scale_factor;
            T2 = 'Noise Free'; 
            tmp_NOISE_REALZ = 1;
             
        end
        MeanImg = zeros(xdim,ydim,ITERATIONS,NUMVAR);
        
        
        %% RECONSTRUCTION USING DIFFERENT KERNELS
        %% Noisy Reconstructions
        delta2 = 0.00000001;
        unit_vector_image = ones(xdim,ydim);
        
            
        for count_var = 1:NUMVAR
            convolFINAL = squeeze(KernelsSet(:,:,count_var));
            convolFINAL_tr = transpose(convolFINAL);
            
            for count_rlz = 1:tmp_NOISE_REALZ
                %disp(' I ');
                disp(['Reconstructing ' , T1 , ', ' , T2 , ', Realization # ',num2str(count_rlz),'  of ',num2str(tmp_NOISE_REALZ ),...
                    ', for PSF # ',num2str(count_var),'  of ',num2str(NUMVAR),' @ ' num2str(toc,6),' sec.']);
                %% Initialization
                image_old = unit_vector_image;
%                 image_old = image_true; warning('Starting recon from the true image. Do Not get surprized by good results!!!');
                
                for count_it = 1:ITERATIONS
                    %                 disp(['Iter = ',num2str(count_it)]);
                    for m = 1:SUBSETS
                        sensit_image = sensit_image_all(:,:,m,count_var);
                        
                        %% Forward model
%                         % Image-space blurring % I'm removing this if
%                         condition because it's of no use for FDG. 
%                         if RECONST_RM
%                             image_old2 = conv2(image_old,filter,'same');
%                         else
                            image_old2 = image_old;
%                         end
                        
                        % Geometric projection
                        [data,xp] = radon(image_old2,theta_m(:,m));
                        
                        % Projection-space blurring
                        if RECONST_RM
                            expected_data = zeros(binsdim,thdim_m);
                            for th = 1:thdim_m
                                for bin = 1:binsdim
                                    expected_data(:,th) = expected_data(:,th) + data(bin,th)*convolFINAL(bin,:)';
                                end
                            end
                        else
                            expected_data = data; %#ok<*UNRCH>
                        end
                        
                        non_zero_exp_data = (expected_data>0);
                        zero_exp_data = 1-non_zero_exp_data;
                        expected_data = expected_data+delta2*zero_exp_data; % just to make sure we're not deviding by zero
                        %ratio = non_zero_exp_data.*Y./expected_data; % for a randomly generated image
                        ratio = non_zero_exp_data.*Y_bar_realization(:,:,m,count_rlz)./expected_data; % to calculate mean image
                        
                        
                        %% Backward model
                        % Projection-space blurring (transposed)
                        if RECONST_RM
                            ratio2 = zeros(binsdim,thdim_m);
                            for th = 1:thdim_m
                                for bin = 1:binsdim
                                    ratio2(:,th) = ratio2(:,th)+ratio(bin,th)*convolFINAL_tr(bin,:)';
                                end
                            end
                        else
                            ratio2 = ratio;
                        end
                        
                        % Geometric projection
                        aaa = iradon(ratio2,theta_m(:,m),'none');
                        correction = aaa(2:xdim+1,2:ydim+1);
                        
                        % Image-space blurring (transposed)
                        if RECONST_RM
                            correction = conv2(correction,filter_tr,'same');
                        end
                        
                        image_new = image_old.*correction./sensit_image;
                        image_new(isnan(image_new)) = 0; %Make sure there's no NaN
                        image_oldd = image_old;
                        image_old = image_new;
                    end
                    
                    
                    
                    %% Saving the boxed images
                    if IsNoisy  ==  2               % Save Noisy Images
                        if     img_prs_abs == 1     % Save Noisy Signal Absent
                            SIG_ABS(:,:,count_it,count_rlz,count_var) = image_new(BkgndBox(1):BkgndBox(2),BkgndBox(3):BkgndBox(4));
                        elseif img_prs_abs == 2    % Save Noisy signal Present
                            SIG_PRS(:,:,count_it,count_rlz,count_var) = image_new(BkgndBox(1):BkgndBox(2),BkgndBox(3):BkgndBox(4));
                        end
                    else                            % Save Noise-Free Images
                        if     img_prs_abs == 1     % Save Noise-Free Signal Absent
                            SIG_ABS_NF(:,:,count_it,count_var) = image_new(BkgndBox(1):BkgndBox(2),BkgndBox(3):BkgndBox(4));
                        elseif img_prs_abs == 2     % Save Noise-Free Signal Present
                            SIG_PRS_NF(:,:,count_it,count_var) = image_new(BkgndBox(1):BkgndBox(2),BkgndBox(3):BkgndBox(4));
                        end
                    end
                    
                    % Save for the mean image
                    MeanImg(:,:,count_it,count_var) = MeanImg(:,:,count_it,count_var) + image_new;
                    
                end
            end
            
            %figure; imagesc(image_new); axis tight; colorbar; colormap(jet); title({['Last Reconstructed Image, Last Iter, Last Realz, PSF #',int2str(count_var),' of ',int2str(NUMVAR)];...
            %    ['Iter #',num2str(count_it),' | Realz #',int2str(count_rlz),' | PSF #',int2str(count_var),' | ',int2str(size(image_new,1)),'\times',int2str(size(image_new,2))]});
        end
        
        %% Save signal present | signal absent | noisy | noise-free
        if IsNoisy  ==  2                   % Save Noisy Images
            if img_prs_abs  ==  1           % Save Noisy Signal Absent
                Mean_SIG_ABS = MeanImg / NOISE_REALZ;
            elseif img_prs_abs  ==  2       % Save Noisy signal Present 
                Mean_SIG_PRS = MeanImg / NOISE_REALZ;
            end
        elseif IsNoisy == 1
            if img_prs_abs == 2
                Mean_SIG_PRS = MeanImg / NOISE_REALZ;
                FileName = ['SIG_PRS_NF__',int2str(xdim),'_R',int2str(NOISE_REALZ),'_I',int2str(ITERATIONS),'_S',int2str(SUBSETS),'_K',int2str(NUMVAR),'_Tmr',int2str(TmrCount),'.mat'];
                Path_str = [DIRECTORY,FileName];
                save (Path_str ,'SIG_PRS_NF' , 'image_true' , '-v7.3');
            elseif img_prs_abs == 1
                Mean_SIG_ABS = MeanImg / NOISE_REALZ;
                FileName = ['SIG_ABS_NF__',int2str(xdim),'_R',int2str(NOISE_REALZ),'_I',int2str(ITERATIONS),'_S',int2str(SUBSETS),'_K',int2str(NUMVAR),'_Tmr',int2str(TmrCount),'.mat'];
                Path_str = [DIRECTORY,FileName];
                save (Path_str ,'SIG_ABS_NF' , 'image_true' , '-v7.3');
            end
         dicomwrite(Mean_SIG_ABS, 'name')
        end
    end
end
clear MeanImg
disp(['Reconstructing ' , T1 , ', ' , T2 , ' is over @ ' num2str(toc,6),' sec.']);

%% Save SIG_PRS, SIG_ABS and data
if SAVE
    if (IMG_ABS_PRS == 1) || (IMG_ABS_PRS == 2)
        FileName = ['SIG_PRS_ALL__',int2str(xdim),'_R',int2str(NOISE_REALZ),'_I',int2str(ITERATIONS),'_S',int2str(SUBSETS),'_K',int2str(NUMVAR),'_Tmr',int2str(TmrCount),'.mat'];
        Path_str = [DIRECTORY,FileName];
        save (Path_str , 'SIG_PRS' , 'SIG_PRS_NF' , 'Mean_SIG_PRS' , 'image_true' , 'TmrCount' , '-v7.3');
    end
    
    if (IMG_ABS_PRS == 0) || (IMG_ABS_PRS == 2)
        FileName = ['SIG_ABS_ALL__',int2str(xdim),'_R',int2str(NOISE_REALZ),'_I',int2str(ITERATIONS),'_S',int2str(SUBSETS),'_K',int2str(NUMVAR),'.mat'];
        Path_str = [DIRECTORY,FileName];
        save (Path_str , 'SIG_ABS','SIG_ABS_NF','Mean_SIG_ABS', '-v7.3');
    end
    
    
    if 0 % Run the following lines to save only the noise-free
%         DIRECTORY = 'F:\Saeed_Shared\SimulationData\XCATsimulation\';
        %DIRECTORY = 'S:\SimulationData\XCATsimulation\';
        %DIRECTORY = 'C:\Users\Roberto\Documents\Research Positions\Qurit Lab\Lymphoma\PET_sim_recon-master';
        DIRECTORY = '';
        % DIRECTORY = 'C:\Users\Saeed\OneDrive_JHMI\OneDrive\Projects\XCAT\';
        FileName = ['SIG_PRS__',int2str(xdim),'_R',int2str(NOISE_REALZ),'_I',int2str(ITERATIONS),'_S',int2str(SUBSETS),'_K',int2str(NUMVAR),'_Tmr',int2str(TmrCount),'_NF.mat'];
        Path_str = [DIRECTORY,FileName];
        save (Path_str ,'SIG_PRS_NF' , 'image_true' , 'TmrCount' , '-v7.3');
    end
end

MSA_obj = cat(3, MSA_obj, Mean_SIG_ABS(:,:,1));

end

s = strcat('../dicom-conversion/input/', file, '_R', num2str(c), '.bin')
n = num2str(s)

fileD = fopen(n, 'w')
fwrite(fileD, MSA_obj, 'single')
fclose(fileD)

clear MSA_obj;

end
