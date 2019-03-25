# CSSA
Analysis and Simulation Tools for CRISPR-Cas9 Pooled Screens

CSSA is an R package for simulation and analysis of pooled CRISPR-Cas9 screens. The simulator is immensely helpful in understanding the impact of, well, all kinds of variables on the outcome of a CRISPR-Cas9 screen. Install the package and read the documentation of CRISPRsim to learn more! Hopefully I will have a vignette added soon. The other functions are all related to analyzing a pooled CRISPR-Cas9 screen. With a few commands you can run the analysis in a manner similar to drugZ, which I have found to be a robust analysis tool. But perhaps more interesting, try the getdeg function, which tries to calculate actual effect sizes caused by gene knockout! Finally, there are a few functions included to calculate rate ratios. The simplest of these can add a number of artificial read counts to all features (guides) to prevent divided-by-zeros. But also give radjust a shot, which gives rate ratios based on confidence limits calculated through a binomial distributions! 

If I can make one last pitch for this package, it is that it has very straightforward steps to analyzing both sensitizer and synthetic lethality screens. These have been especially lacking in many of the previously available methods.

For any questions, contact me at j.poell@vumc.nl
