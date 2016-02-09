# What is the LETKF? #

  * LETKF stands for Local Ensemble Transform Kalman Filter, invented by [Hunt et al. (2007)](http://dx.doi.org/10.1016/j.physd.2006.11.008) at the [University of Maryland](http://www.umd.edu/).
  * It is an advanced data assimilation method for many possible applications.
  * It has been tested with numerical weather prediction (NWP) models, storm to global scales.

# What is available here? #

  * LETKF source codes and run scripts
  * Some documentation

# Feature #

  * All-in-one package to run SPEEDY-LETKF, a great tool to learn and play with.
  * Written by Fortran 90/95 language and B-shell scripting
  * Parallel processing with the MPI (Message Passing Interface) library and OpenMP directives, making most use of parallel architecture computers including multi-core processors, cluster environment, and even supercomputer (tested with the Japanese [Earth Simulator](http://en.wikipedia.org/wiki/Earth_Simulator)).
  * It would be relatively easy to adapt to other applications.
  * Any kinds of computer simulations of real world that we can observe.
  * Possible collaborators, please contact [Takemasa Miyoshi](http://www.atmos.umd.edu/~miyoshi/)!

# Getting Started #

  * I recommend general users start with [SPEEDY-LETKF](SpeedyLetkf.md).
  * Those who are already familiar with the codes may go directly to [Source](http://code.google.com/p/miyoshi/source/checkout).

# Mailing lists #

  * Announcements of important changes: [miyoshi-code-announce](http://groups.google.com/group/miyoshi-code-announce)
> This is just for announcements. Users cannot post emails.
  * Discussions (comments, questions, etc.): [miyoshi-code-discussions](http://groups.google.com/group/miyoshi-code-discussions)
> This is for exchanging comments. User's comments and questions are welcome. Anyone can answer. I would like to see users help each other.
  * Developers: [miyoshi-code-devel](http://groups.google.com/group/miyoshi-code-devel)
> This is a closed list for developers. Automated announcements of commits and other actions to the project, and discussions among developers. If you would like to receive detailed announcements of each commit and other action to the project, please contact [Takemasa Miyoshi](http://www.atmos.umd.edu/~miyoshi/).

# About this google code project #

  * Maintained by [Takemasa Miyoshi](http://www.atmos.umd.edu/~miyoshi/).
  * All contents are provided "as-is" without warranties of any kind. I welcome feedbacks, but please do not expect support.

# References #

**System upgrades**

  * Miyoshi, T., 2010: The Gaussian Approach to Adaptive Covariance Inflation and Its Implementation with the Local Ensemble Transform Kalman Filter. Mon. Wea. Rev., in press. doi:10.1175/2010MWR3570.1
  * Miyoshi, T., S. Yamane, and T. Enomoto, 2007: Localizing the Error Covariance by Physical Distances within a Local Ensemble Transform Kalman Filter (LETKF). SOLA, 3, 89-92. doi:10.2151/sola.2007-023

**WRF-LETKF**

  * Miyoshi, T. and M. Kunii, 2011: The local ensemble transform Kalman filter with the Weather Research and Forecasting model: experiments with real observations. Pure and Appl. Geophys., in press.

**AFES-LETKF**

  * Miyoshi, T. and S. Yamane, 2007: Local ensemble transform Kalman filtering with an AGCM at a T159/L48 resolution. Mon. Wea. Rev., 135, 3841-3861. doi:10.1175/2007MWR1873.1
  * Miyoshi, T., S. Yamane, and T. Enomoto, 2007: The AFES-LETKF Experimental Ensemble Reanalysis: ALERA. SOLA, 3, 45-48. doi:10.2151/sola.2007-012

**JMA**

  * Miyoshi, T., Y. Sato, and T. Kadowaki, 2010: Ensemble Kalman filter and 4D-Var inter-comparison with the Japanese operational global analysis and prediction system. Mon. Wea. Rev., 138, 2846-2866. doi:10.1175/2010MWR3209.1
  * Miyoshi, T. and Y. Sato, 2007: Assimilating Satellite Radiances with a Local Ensemble Transform Kalman Filter (LETKF) Applied to the JMA Global Model (GSM). SOLA, 3, 37-40. doi:10.2151/sola.2007-010
  * Miyoshi, T. and K. Aranami 2006: Applying a Four-dimensional Local Ensemble Transform Kalman Filter (4D-LETKF) to the JMA Nonhydrostatic Model (NHM). SOLA, 2, 128-131. doi:10.2151/sola.2006-033

**Initial version**

  * Miyoshi, T., 2005: Ensemble Kalman filter experiments with a primitive-equation global model. Ph.D. dissertation, University of Maryland, College Park, 197pp.