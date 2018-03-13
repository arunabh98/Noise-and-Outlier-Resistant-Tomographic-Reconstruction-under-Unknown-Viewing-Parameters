
ARPord.m does the bulk of the work. However, to run and configure for multiple images, use Driver.m
In Driver.m, all parameters are neatly stacked at the top in the PARAMETERS section

	imagelist: list of image files to be consumed. Is prepended by ‘../images/’ automatically
	noiselist: For each image, each of these noise levels will be tested
	numkeep: number of (distinct) angles to retain out of the possible 180
	numstarts: number of starts to use in the multistart approach for coordinate descent
	Ord: Highest order moments to consider
	randseed: random seed


Other files in the directory:
	assembleA.m: returns the A matrix (as in A*IM = PM) for a given set of angles and moment order
	calculateProjectionMoment.m: calculates the k^th moment of the data P
	denoise.m: Given noisy projections, returns denoised projections, using patch-based PCA denoising
	imageMomentFromImage.m: calculates image moments directly from image (to compare)

	orderbytsp.m, tsp_ga.m: Order the projections as per TSP. Not currently used in ARPord



