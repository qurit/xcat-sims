Summary:

• Analytic simulation of PET data (with and without noise)

• Begin with your own image, or you can use our included NEMA NU-2 phantom generator

• Enables noisy and noise-free reconstructions. 

• Statistical PET image reconstruction framework using 2D ordered-subset expectation maximization (OS-EM) algorithm

• Incorporate attenuation and normalization modeling and correction using a CT sinogram and detector normalization map, respectively.

• Integrated generalized PSF modeling (resolution modeling) to model true PSF, no-PSF, as well as under- and overestimated PSF in reconstruction

• Can vary a range of data acquisition and image reconstruction parameters

• Capable of simulating and reconstructing signal (tumor) absent, i.e. replacing the signal value with the background, as well as signal present, to study the effects 
of model parameters on the background in the exact same location as the tumor.



Reference:

Please use the following reference if you publish results obtained using this software tool: 

S. Ashrafinia, et al., “Generalized PSF modeling for optimized quantitative-task performance”, Phys. Med. Biol., vol. 62, pp. 5149-5179, 2017.



Technical Description:

This Matlab®-based framework incorporates PET forward projection and reconstruction steps. 

The framework models many realistic PET reconstruction processes, such as decay of radioactivity, ability to perform attenuation correction using a CT sinogram as well as detector sensitivity normalization using normalization sinogram. The framework is capable of reconstructing noisy and/or noise-free, as well as signal present/absent. For the noisy reconstruction it incorporates realistic random Poisson noise to the data to generate different noise realizations.

A single slice can be loaded as the true image. The PSF modeling has different settings including no PSF, true PSF, and our proposed generalized PSF modeling that uses under-and overestimated PSF kernels in the reconstruction as described in the paper above. The true PSF is derived analytically from modeling PET resolution degradation effects including photon non-collinearity, inter-crystal scattering and inter-crystal penetration, and is being used in the forward projection. In the reconstruction, either the true PSF can be used to model PSF reconstruction. Other PSF kernels are under- or overestimated version of this analytically derived true PSF that was used in the forward projection.

Images are saved at the end of every iteration, for every noise realization, and every PSF setting. The framework provides beginning with (forward projecting) a higher resolution image, while reconstructing a lower resolution image, which is a more realistic representation of imaging objects.

We also include our digital NEMA phantom generation code that generates a NEMA NU-2 image quality phantom according to NEMA standards to use as a true image. 

Please being with "Main_Recon.m". More explanation in given in the comments section of each file.


Contact

For support, contributions and questions, please contact s.ashrafinia (at) gmail (dot) com

