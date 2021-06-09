![img](https://user-images.githubusercontent.com/54952340/121406995-b00aee00-c913-11eb-9b93-c6782fcb61f6.png)
# xcat-sims

Summary:

•	Novel lymphatic system was defined for the 4D-extended cardiac torso (XCAT) phantom

•	A pipeline was developed to simulate and reconstruct PET images from the XCAT phantom, applying header information that is compatible with clinical radiology software

•	As example application, simulated lymphoma patients were modelled using XCAT phantom and can be freely accessed

Reference:

Please use the following references if you publish results with help from this software tool:

R. Fedrigo, et al., “Development of the Lymphatic System in the 4D XCAT Phantom for Improved Multimodality Imaging Research”, J. Nuc. Med., vol. 62, publication 113, 2021.

The Matlab component of this framework is adapted from the PET simulation and image reconstruction tool, as described by:

S. Ashrafinia, et al., “Generalized PSF modeling for optimized quantitative-task performance”, Phys. Med. Biol., vol. 62, pp. 5149-5179, 2017.

Technical Description:

The novel lymphatic system for the 4D-extended cardiac torso (XCAT) phantom enhances the ability to model diseases, such as lymphoma. The lymphatic system (nodes, vessels) was defined using non-uniform rational basis spline (NURBS) surfaces. Multichannel large deformation diffeomorphic metric mapping (MC-LDDMM) method was used to propagate from the template phantom to different XCAT anatomies. The XCAT general parameter script was used to generate files that define the ground truth radioactivity and attenuation for a simulated patient.

Example files are provided, which can be used to model patients with non-Hodgkin’s lymphoma.  Bulky tumours were created by altering lymph node morphology and function (nodes were expanded, stretched, converged, and had increased tracer uptake). This methodology is highly flexible and may be used in virtual clinical trials to simulate nodes with different pathology.

A framework was developed in Matlab and Python to simulate and reconstruct PET images using patients modelled with the XCAT phantom. Output files are converted to dicom and necessary header information is applied, such that the images can be viewed using clinical radiology software.

Compatibility note:

The MATLAB scripts were executed using the following toolboxes and versions.

MATLAB, version                           9.8

Simulink                                  10.1

Computer Vision Toolbox                   9.2 

Curve Fitting Toolbox                     3.5.11

DSP System Toolbox                        9.10

Deep Learning Toolbox                     14.0

Fuzzy Logic Toolbox                       2.7

Image Processing Toolbox                  11.1

Signal Processing Toolbox                 8.4

Statistics and Machine Learning Toolbox   11.7

Symbolic Math Toolbox                     8.5

Contact:

For questions, recommendations, or feedback, please contact: rfedrigo@student.ubc.ca
