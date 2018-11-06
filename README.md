# Noise-and-Outlier-Resistant-Tomographic-Reconstruction-under-Unknown-Viewing-Parameters

This repository contains the implementation of the novel algorithm presented in this [paper][1].

## Abstract
The paper presents an algorithm for effectively reconstructing an object from a set of its tomographic projections without any knowledge of the viewing directions or any prior structural information, in the presence of pathological amounts of noise, unknown shifts in the projections, and outliers among the projections. The outliers are mainly in the form of a number of projections of a completely different object, as compared to the object of interest. We introduce a novel approach of first processing the projections, then obtaining an initial estimate for the orientations and the shifts, and then define a refinement procedure to obtain the final reconstruction. Even in the presence of high noise variance (up to 50% of the average value of the (noiseless) projections) and presence of outliers, we are able to successfully reconstruct the object. We also provide interesting empirical comparisons of our method with the sparsity-based optimization procedures that have been used earlier for image reconstruction tasks.

[1]: https://arunabh98.github.io/reports/tomography.pdf
