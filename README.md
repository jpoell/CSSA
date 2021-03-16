# CSSA
Analysis and Simulation Tools for CRISPR-Cas9 Pooled Screens

CSSA is an R package for simulation and analysis of pooled CRISPR-Cas9 screens. The simulator is immensely helpful in understanding the impact of, well, all kinds of variables on the outcome of a CRISPR-Cas9 screen. Install the package and read the documentation of CRISPRsim and sortingsim to learn more! Hopefully I will have a vignette added at some point. The other functions are all related to analyzing a pooled CRISPR-Cas9 screen. With a few commands you can run the analysis in a manner similar to drugZ, which I have found to be a robust analysis tool. But perhaps more interesting, try the getdeg function (or degrep to incorporate replicate arms), which tries to calculate actual effect sizes caused by gene knockout! The geteffect function also zooms in on effect sizes, but taking a likelihood approach. There are a few functions included to calculate rate ratios. The simplest of these can add a number of artificial read counts to all features (guides) to prevent divided-by-zeros. But also give radjust a shot, which gives rate ratios based on confidence limits calculated through a binomial distributions! Finally, I have included an analysis method based on odds and nonparametric distribution of guides. Give it a whirl! I expect it is most useful for selection-based screens (I specifically made the functions for analysis of a FACS-based CRISPR screen), but there is no reason why it might not work well in other contexts.

If I can make one last pitch for this package, it is that it has very straightforward steps to analyzing both sensitizer and synthetic lethality screens. These have been especially lacking in many of the previously available methods!

For any questions, contact me at j.poell@amsterdamumc.nl
