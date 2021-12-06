This code is the implementation of our efficient data-driven sparse
identification algorithm of dynamical systems.
The code aims at reconstructing the different sets of governing equations
and identifying discontinuity surfaces in hybrid systems when the number of
discontinuities is known a priori.

The main codes named maincode_FIEH_ex1.m and maincode_FIEH_ex2_transient_data.m is autonomous.
Running these main codes both the training and validation steps are
execute.

The folders "plot_performance" contain additional code to perform more in-depth analysis about how the algorithm works.

The folders images contain some intresting plot related to the aforementioned algorithm.

The other folders ("FIEH", "SINDy", "validation") contain the utils
functions called by the main code.

A more in-depth validation and analysis of the code could be find in the
script called "ParetoFrontApproxModel.m"