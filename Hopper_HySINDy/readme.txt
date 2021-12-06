This code is the implementation of our efficient data-driven sparse
identification algorithm of dynamical systems.
The code aims at reconstructing the different sets of governing equations
and identifying discontinuity surfaces in hybrid systems when the number of
discontinuities is known a priori.

The main code named maincode_Hopping.m is autonomous.
Running maincode_Hopping.m both the training and validation steps are
execute.

The folders "plot_performance" and "Analysis_RecoveryRate_ClusterDimension_Noise"
contain additional codes to perform more in-depth analysis about how the algorithm works.

The other folders ("hopper", "SINDy", "validation") contain the utils
functions called by the main code.
