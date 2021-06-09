function [P] = GenNEMA(ImDim , VoxSize , ROIs)
% ===========================
% This code generates NEMA phantom for a given dimension and voxel size.
% The values of each ROI also can be assigned.
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

%% Initialization
% ImDim = 128;
% VoxSize = 3.47;
% ROIs = [1 2 3 4 5 6 7]; %seven numbers to assign to background and six ROIs [BK , ROI_1 , ... , ROI_6]
ROIs(2:end) = ROIs(2:end) - ROIs(1);

%% Draw Upper semi-circle
shapeInserter=vision.ShapeInserter('Fill' , 10 , 'Shape' , 'Circles' , 'FillColor' , 'Black');
C1 = int32([ImDim/2 , (ImDim/2 + round(35/VoxSize)) , round(147/VoxSize) ]);
J1=step(shapeInserter , ones(ImDim) , C1);

shapeInserter=vision.ShapeInserter('Fill' , 10 , 'Shape' , 'Rectangles'  , 'FillColor' , 'White');
R1 = int32([1 , (ImDim/2 + round(35/VoxSize))+1 , ImDim-1 , ImDim/2 ]);
J2=step(shapeInserter , J1 , R1);

P1=J2;
P1(J2==J2((ImDim/2 + round(35/VoxSize))+2 , (ImDim/2 + round(35/VoxSize))+2)) = J2(1 , 1);
P1(P1~=1)=0;
P1 = ~P1;

%% Draw lower two circles and rectangle
shapeInserter=vision.ShapeInserter('Fill' , 10 , 'Shape' , 'Circles' , 'FillColor' , 'Black');
C2 = int32([(ImDim/2 - round(70/VoxSize)) , (ImDim/2 + round(35/VoxSize)) , round(77/VoxSize) ; (ImDim/2 + round(70/VoxSize)) , (ImDim/2 + round(35/VoxSize)) , round(77/VoxSize)]);
J = step(shapeInserter , ones(ImDim) , C2);

shapeInserter = vision.ShapeInserter('Fill' , 10 , 'Shape' , 'Rectangles'  , 'FillColor' , 'Black');
R2 = int32([(ImDim/2 - round(70/VoxSize)) , (ImDim/2 + round(35/VoxSize)) , ceil(140/VoxSize) , ceil(77/VoxSize)  ]);
J3 = step(shapeInserter , J , R2);

P2 = J3;
P2(J3~=1) = 0;
P2 = ~P2;
BK = P2|P1;

%% Draw six circles
shapeInserter=vision.ShapeInserter('Fill' , 10 , 'Shape' , 'Circles' , 'FillColor' , 'Black');
C2 = int32([(ImDim/2 + round(114.4/2/VoxSize)) , ImDim/2 , ceil(37/2/VoxSize)]);
J_1 = step(shapeInserter , ones(ImDim) , C2);
J_1(J_1 ~= 1) = 0;
J_1 = ~J_1;

shapeInserter=vision.ShapeInserter('Fill' , 10 , 'Shape' , 'Circles' , 'FillColor' , 'Black');
C2 = int32([(ImDim/2 - round(114.4/2/VoxSize)) , ImDim/2 , ceil(17/2/VoxSize)]);
J_4 = step(shapeInserter , ones(ImDim) , C2);
J_4(J_4 ~= 1) = 0;
J_4 = ~J_4;

shapeInserter=vision.ShapeInserter('Fill' , 10 , 'Shape' , 'Circles' , 'FillColor' , 'Black');
C2 = int32([ImDim/2 + round((114.4/2/VoxSize) * cos(60 * pi/180)) , ImDim/2 - round((114.4/2/VoxSize) * sin(-60 * pi/180)) , ceil(28/2/VoxSize)]);
J_2 = step(shapeInserter , ones(ImDim) , C2);
J_2(J_2 ~= 1) = 0;
J_2 = ~J_2;

shapeInserter=vision.ShapeInserter('Fill' , 10 , 'Shape' , 'Circles' , 'FillColor' , 'Black');
C2 = int32([ImDim/2 - round((114.4/2/VoxSize) * cos(60 * pi/180)) , ImDim/2 + round((114.4/2/VoxSize) * sin(60 * pi/180)) , ceil(22/2/VoxSize)]);
J_3 = step(shapeInserter , ones(ImDim) , C2);
J_3(J_3 ~= 1) = 0;
J_3 = ~J_3;

shapeInserter=vision.ShapeInserter('Fill' , 10 , 'Shape' , 'Circles' , 'FillColor' , 'Black');
C2 = int32([ImDim/2 - round((114.4/2/VoxSize) * cos(60 * pi/180)) , ImDim/2 - round((114.4/2/VoxSize) * sin(60 * pi/180)) , ceil(13/2/VoxSize)]);
J_5 = step(shapeInserter , ones(ImDim) , C2);
J_5(J_5 ~= 1) = 0;
J_5 = ~J_5;

shapeInserter=vision.ShapeInserter('Fill' , 10 , 'Shape' , 'Circles' , 'FillColor' , 'Black');
C3 = int32([ImDim/2 + round((114.4/2/VoxSize) * cos(60 * pi/180)) , ImDim/2 - round((114.4/2/VoxSize) * sin(60 * pi/180)) , ceil(10/2/VoxSize)]);
if C2(end) <= C3(end)
    shapeInserter = vision.ShapeInserter('Fill' , 10 , 'Shape' , 'Rectangles'  , 'FillColor' , 'Black');
    C3 = int32([ImDim/2 + round((114.4/2/VoxSize) * cos(60 * pi/180))-1 , ImDim/2 - round((114.4/2/VoxSize) * sin(60 * pi/180))-1 ,...
         -1+2*ceil(10.5/2/VoxSize) , -1+2*ceil(10.5/2/VoxSize)]);
end
J_6 = step(shapeInserter , ones(ImDim) , C3);
J_6(J_6 ~= 1) = 0;
J_6 = ~J_6;
% imagesc(BK+J_1+J_4+J_6+J_2+J_5+J_3);

%% Assigning values to each ROI
P = ROIs(1) * BK + ROIs(2) * J_1 + ROIs(3) * J_2 + ROIs(4) * J_3 + ROIs(5) * J_4 + ROIs(6) * J_5 + ROIs(7) * J_6 ;
% imagesc(P)


