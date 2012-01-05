Automated Parameter Estimation and Model Selection Toolkit (APEMoST)
=====================================================================

If you are unfamiliar with Bayesian inference using MCMC (Markov Chain Monte Carlo) sampling/updates, you can find some information in `the Links section <others.html#introduction>`_.

The program was first written by Michael Gruberbauer (at the time University of Vienna).
Later it was rewritten for parallelization by Johannes Buchner (at the time University of Vienna, Technical University of Vienna), which was funded by Werner W. Weiss (University of Vienna). 
Now it is maintained by Johannes Buchner.

Specification
~~~~~~~~~~~~~~~

APEMoST enables the user to perform parameter estimation in a simple way: 
Only three things have to be provided:

#. The likelihood function
   
	Specified as a C function (usually ~15 lines of code)
   
#. The parameter space
   
	Specified as a text file, containing the borders and starting points
   
#. A data file
   
	A table of data the likelihood function can work on is read in and 
	available as a matrix.


Fundamentals
~~~~~~~~~~~~~~~

APEMoST is written in ANSI C and uses the GNU Scientific Library (GSL). This makes it a very fast sampler. 
In fact, if run on a million iterations, less than a minute is spent in APEMoST, the rest of the time the CPU evaluates the likelihood function.

APEMoST utilizes all CPUs on a computer, by using OpenMP. 

Features
~~~~~~~~~~

- **Automatic stepwidth calibration** of the proposal distribution. 
    
    The stepwidth of each parameter is optimally calibrated for the desired acceptance rate.

- **Simulated Tempering and Parallel Tempering**: 
    
    A number of chains (usually 20) is run and swap their positions. This can be disabled 
    by setting the number of chains to 1.

- **Automagic**: 

    For example, the frequency of swaps, the betas of the hottest chain and many more 
    values are calculated to have a sane value the user usually doesn't have to mess with.
    
    Furthermore, we discovered a way of predicting the stepwidths of chains once two chains
    are calibrated. This allows skipping time consuming calibration and yet reaching the 
    desired acceptance rates

We also have some `differences to other software <other.html>`_.

Ready to dig in? **Get started** with the `user manual <manual.html>`_, it explains APEMoST step by step.

