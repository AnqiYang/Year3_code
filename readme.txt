### Codebase for Spatially-varying Gain and Binning

There are two main scripts:

- main_binning_theory.m
  This script implements the theory that given light level decides optimal binning size.
	
- main_sv_gain_binning.m
  This script simulate noisy images from an HDR scene with these modes:
  (1), ground-truth noise free
  (2), no binning + const gain
  (3), no binning + spatially-varying (ROI) gain
  (4), no binning + per-pixel gain
  (5), additive binning + const gain
  (6), additive binning + spatially-varying gain
  (7), average binning + const gain
  (8), average binning + spatially-varying gain
  (9), ISP (digital) binning + const gain
  (10), ISP (digital) binning + spatially-varying gain


