.. include:: header.txt

-------------------------------------------------------------
User Manual
-------------------------------------------------------------

A general introduction, frequently asked questions and the license can be found at `the main site <index.html>`_. The website of APEMoST is http://apemost.sourceforge.net/. The latest version of this manual can be found there aswell.

.. contents::

--------------------------
 Prerequisites
--------------------------

You need to have the following software installed:

- gcc (>= 4.1, for OpenMP)
- GNU scientific library (libgsl0-dev or gsl-devel)
- Boehm garbage collector (boehmgc or gc)
	
	The garbage collector is optional, but can be enabled by compiling with
	make WITH_GARBAGE_COLLECTOR=1 

For the version control (keeping track of different versions, experimenting with code and patches)

- git (git-core or git)

For generating the documentation (make doc):

- doxygen 
- docutils (python-docutils)
- rst2pdf

It may be handy to keep this manual open, and the `GSL manual <http://www.gnu.org/software/gsl/manual/>`_ as well as the `api documentation <api/html/index.html>`_.

----------------------------
 Getting APEMoST
----------------------------

.. image:: http://sflogo.sourceforge.net/sflogo.php?group_id=275812&type=5
  :height: 62
  :width: 210
  :alt: Sourceforge project page of APEMoST
  :target: https://sourceforge.net/projects/apemost/
  :align: right

You can get a release of APEMoST from the `Sourceforge project page http://sourceforge.net/projects/apemost/ <http://sourceforge.net/projects/apemost/>`_

or `fetch the latest code using git <https://sourceforge.net/scm/?type=git&group_id=275812>`_.


---------------------------
 The general concept
---------------------------

A quick overview of the components. All of these will be explained in detail in the following sections.

There are 2 pieces of software:

#. The MCMC engine, APEMoST

	Its code can be found in the folder src/

#. The user-provided likelihood function

	The code has to be placed in apps/

Say your project is called "foo".
The idea is that you implement one function, calc_model, in apps/foo.c .
Then, you run 

	$ CCFLAGS="-DMAX_ITERATIONS=1000000" make foo.exe

and this compiles APEMoST and your function together. All code is compiled at once to get the best optimization.
You can specify the algorithms behaviour and many other things with CCFLAGS. This will be explained later. 

Furthermore, when you APEMoST (that has your likelihood function in it), you need two more files that you need to place where you start your program foo.exe (in the working directory). They are called data and params.

.. _params:

#. The definition of the parameter space, "params"

	Each line characterizes a dimension of parameter space. The columns are:
	
	#. Start value
	#. Min value
	#. Max value
	#. Name (a string)
	#. Initial step width (before calibration)
		
	Use a name without whitespaces. You can set the step width to -1, which results in automatically using
	10% of (max - min).
	
	.. _data:
	
#. A data table, "data"
	
	This file will be read into a matrix that is accessible for the likelihood calculation.
	
	If you have your data elsewhere, you can also use a symlink.


General remarks
~~~~~~~~~~~~~~~~~~~~~~~

The program operates on logarithmic probability values (loglikelihood). 
Everywhere and exclusively. Double precision is used.

The values read from files can be of many formats (everything scanf can read). 


------------------------
Phase 0 - Preparation
------------------------

In this phase of your progress, you will want to

- specify your likelihood function from a formula

- compile the program the first time

- get used to CCFLAGS

-----------------------------------------------------------

Equipped with the GSL manual, take a look at the file apps/simplesin.c. 
You will write something very similar.

You have to implement two functions::

	void calc_model(mcmc * m, const gsl_vector * old_values);
	
	void calc_model_for(mcmc * m, const unsigned int i, const double old_value) {
		calc_model(m, NULL);
	}

The second one, calc_model_for, is only used in the calibration: Only the ith parameter has been changed. 
The old values are given as parameters in case you can do some optimizations 
(e.g. if the parameter is just an offset, simply add/subtract something from the probability).
The default (as above) is to just recalculate everything it using calc_model().

The first function, calc_model() has to 

- look at the parameter values

	get_params(m) returns all parameters as a vector.

	get_params_for(m, i) returns the ith parameter value.

- look at the data table

	m->data is the matrix read in from the file "data".

	You can do read operations on this matrix, e.g. gsl_matrix_get(m->data, i, j)

- calculate and set the prior
	
	The program has to keep track of the prior, since we need the probability both with
	and without the prior.
	
	Use set_prior(m, myprior);
	
- calculate and set the probability
	
	The probability has to contain the prior, and has to incorporate beta.
	
	Example: set_prob(m, get_prior(m) + get_beta(m) * myprob);

	As your MCMC papers will tell you, the priors should not be exponentiated by beta.

Keep in mind that everything is logarithmic (loglikelihoods!).

*m*, more precisely the mcmc structure, represents one chain.

First compilation
~~~~~~~~~~~~~~~~~~~~~~~~~

Lets try to compile your program. I'll take simplesin as an example (replace simplesin with your project name).

I run::

	$ make simplesin.exe
	$ make eval_simplesin.exe
	$ make benchmark_simplesin.exe

If everything works out, I get three executables: simplesin.exe, eval_simplesin.exe and benchmark_simplesin.exe (replace simplesin with your project name).

Lets see if our loglikelihood function is correct, and evaluate it at.

We change into a empty directory we want to work from, and put two files there. 

In "data" (for example, also see data_)::

	101	0.67
	102	1.01
	103	7.9e-1
	104	1.34
	and so on

In "params" (see above at params_)::
	
	0	0	2	amplitude	-1
	0	0	0.3	frequency	-1
	0	0	1.0	phase	-1
	0	0	2	offset	-1

You can specify the values in different formats (e.g. 0.13, 1.3e-1) and use tabs or spaces as you like (I would recommend tabs).

.. _eval:

Now we can try out the likelihood function (replace apemost-directory with 
where the apemost code and the Makefile is)::

	$ apemost-directory/eval_simplesin.exe
	(you enter:) 	1 0.2 1 0
	(output:) 	-1.480898044165363e+01	0.000000000000000e+00
	(quit with Ctrl-C or Ctrl-D)

The first is the probability, the second the prior.

You may get this error, which can be a little confusing::

	gsl: ../gsl/gsl_vector_double.h:177: ERROR: index out of range
	Default GSL error handler invoked.
	Aborted

This means you tried to access a element beyond the size of the vector (or matrix). 
In that case, the function expects a different number of parameters than the params file provides.

.. _benchmark:

Although this is less relevant for the first read, you can also benchmark your likelihood function with 
the benchmark_simplesin.exe you produced. It takes the number of evaluations as arguments.

The third way of accessing the MCMC engine is the really interesting one::

	$ apemost-directory/simplesin.exe
	$ apemost-directory/simplesin.exe check

This also outputs some inline help about the phases.

You can find the main() functions of these three programs in apps/generic_main.c, apps/eval.c and apps/benchmark.c.

CCFLAGS
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**It is essential that you understand the CCFLAGS variable**. This will be the main "interface"
how you tinker with the program, change its default values and its behaviour.

For compilation, you can do something like::

	CCFLAGS="-DMAX_ITERATION=100000 -DWITHOUT_GARBAGE_COLLECTOR" make simplesin.exe

This tells the compiler to set preprocessor values. Here, I call e.g. WITHOUT_GARBAGE_COLLECTOR 
a "flag", and you "set the flag" by appending it to your CCFLAGS string with "-DFLAG"
and you "set the flag to a value" using -DFLAG=value.

The check subcommand outputs the values currently set (after compilation)::

	$ apemost-directory/simplesin.exe check

A full list of flags can be found in the `api documentation <api/html/index.html>`_, with their 
meaning and default values. This is a good resource that you should keep open.

If you were to write a new calibration algorithm, or use a different adaptive MCMC algorithm, 
you would use "#ifdef MYFLAG" preprocessor directives and enable/disable the use of the algorithm by
a flag.

The perhaps most important flag is DEBUG, which enables some debug output. 

**Note**: Smart readers will notice that you have to rebuild the program
when you want to change a flag something. 

---------------------------------------------------------------------

In this phase of your progress, you learned how to

- specify your likelihood function from a formula

- compile the program the first time

- get used to CCFLAGS

Very good! You get a cookie.

-----------------------------------------------
Phase 1 - Calibration, first chain (beta = 1)
-----------------------------------------------

A good MCMC sampling should have a good acceptance rate. Different sources state 
different things, something between 30% and 80% should be right.

To reach this acceptance rate, a calibration algorithm tinkers with the stepwidths of
the proposal distribution (lets assume the default, a multivariate normal distribution).

The inline help shows which flags are relevant::

	$ apemost-directory/simplesin.exe help calibrate_first

You will want to enable the DEBUG flag, otherwise you won't see much if stuff goes wrong.
Ideally, you don't have to care about it, practically you will want to see which stepwidths 
scale up, which scale down.

You can run the calibration with::

	$ apemost-directory/simplesin.exe calibrate_first

The result of the calibration will be a file that stores the calibrated stepwidths, "calibration_results".
The rows are defined as::
	
	beta	param1_stepwidth	param2_stepwidth	...	param1_value	param2_value	...

Each chain will get one such row in the next phase. For now, just one row in this file.

Also, the program suggests a new params file ("params_suggested") that contains 
the new stepwidths (last column).
If you use these stepwidths in your params file, this will make your next calibrate_first run go faster.

Possible problems in this phase
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- stepwidth gets too large
	
	You may want to increase the parameter space.

	This can also mean the posterior distribution is independent of this parameter!

- stepwidth gets too small

- calibration fails
	
	You can increase ITER_LIMIT.

- calibration takes too long and doesn't find a good end point.

	Bad. 

	Among many things, you can try altering MAX_AR_DEVIATION.

	Among the less recommended, but possible solutions are: 
	
		- manually setting some stepwidths

	You can also add another calibration algorithm to APEMoST (we'd be happy).

**Note**: You can watch the progress of the calibration by plotting the file "calibration_progress.data".
The columns are defined as:
	
	#. parameter number (starting with 0)
	#. number of iteration done
	#. stepwidth ([0..1], normalized to parameter space)
	#. acceptance rate
	#. accuracy of the acceptance rate estimate (-1 if not available)

For example::

	0	200	0.058824	0.590000	-1.000000
	1	200	0.050000	0.590000	-1.000000
	2	200	0.058824	0.590000	-1.000000
	3	200	0.058824	0.590000	-1.000000
	0	400	0.069204	0.395000	-1.000000
	and so on

------------------------------------------------
Phase 1 - Calibration, other chains (beta <= 1)
------------------------------------------------

Now we just have to do the same with the hot chains.

There are some interesting facts about the hot chains in APEMoST, for example

#. Per default, beta is not distributed equally, but using a chebyshev scheme

	This proved to be quite good so far, I tested out several methods (also see my 
	bachelor thesis).

#. The hottest chain's beta is automatically determined so that the stepwidths will be maximally the size of the parameter space.
	
	If the beta_0 seems suspicously low, you can set the flag BETA_0 to something sensible,
	0.01 is often used.
	
#. There is a mechanism that allows skipping the calibration of all but two chains
	
	This saves you plenty of time, and possible calibration failures with the 
	hottest chain. 
	
	You can enable it with the flag SKIP_CALIBRATE_ALLCHAINS.
	
	Although this technique, developed by us (see my bachelor thesis), is not based
	on a sound mathematical proof (yet?), I have yet to see a scenario where this technique
	is inappropriate.

If you want a different number of chains, set N_BETA.

With our knowledge from the previous chapter, we look up the inline help::

	$ apemost-directory/simplesin.exe help calibrate_rest

and run::

	$ apemost-directory/simplesin.exe calibrate_rest


**Note**: You should also be aware that when you change a flag, the likelihood function or
the parameter space, you may have to do the calibration again, as the stepwidths will not 
be appropriate anymore.

Example output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A output of the calibration phase can look like this (ideal case, DEBUG turned off, SKIP_CALIBRATE_ALLCHAINS turned on)::

	$ ../simplesin.exe calibrate_first
	Initializing 20 chains
	Starting markov chain calibration
	wrote calibration results for 1 chains to calibration_results
	new suggested parameters file has been written
	$ ../simplesin.exe calibrate_rest
	Initializing 20 chains
	Calibrating chains
	Calibrating second chain to infer stepwidth factor
		Chain  1 - beta = 0.993277	steps: Vector4d[0.052377;0.000060;0.036247;0.037842]
	stepwidth factors: Vector4d[0.887411;0.887411;1.044013;0.887411]
	automatic beta_0: 0.013294
		Chain  1 - beta = 0.993271	steps: Vector4d[0.046480;0.000053;0.037842;0.033582]
		Chain  2 - beta = 0.973269	steps: Vector4d[0.046955;0.000054;0.038229;0.033925]
		Chain  3 - beta = 0.940538	steps: Vector4d[0.047765;0.000055;0.038889;0.034510]
		Chain  4 - beta = 0.895972	steps: Vector4d[0.048939;0.000056;0.039844;0.035358]
		Chain  5 - beta = 0.840786	steps: Vector4d[0.050519;0.000058;0.041131;0.036500]
		Chain  6 - beta = 0.776486	steps: Vector4d[0.052569;0.000060;0.042800;0.037981]
		Chain  7 - beta = 0.704825	steps: Vector4d[0.055177;0.000063;0.044923;0.039866]
		Chain  8 - beta = 0.627758	steps: Vector4d[0.058466;0.000067;0.047601;0.042242]
		Chain  9 - beta = 0.547388	steps: Vector4d[0.062611;0.000072;0.050976;0.045237]
		Chain 10 - beta = 0.465906	steps: Vector4d[0.067866;0.000078;0.055254;0.049033]
		Chain 11 - beta = 0.385536	steps: Vector4d[0.074605;0.000085;0.060741;0.053902]
		Chain 12 - beta = 0.308470	steps: Vector4d[0.083405;0.000095;0.067906;0.060260]
		Chain 13 - beta = 0.236809	steps: Vector4d[0.095192;0.000109;0.077502;0.068776]
		Chain 14 - beta = 0.172508	steps: Vector4d[0.111531;0.000128;0.090805;0.080581]
		Chain 15 - beta = 0.117323	steps: Vector4d[0.135241;0.000155;0.110109;0.097712]
		Chain 16 - beta = 0.072756	steps: Vector4d[0.171737;0.000197;0.139823;0.124080]
		Chain 17 - beta = 0.040026	steps: Vector4d[0.231543;0.000265;0.188514;0.167290]
		Chain 18 - beta = 0.020023	steps: Vector4d[0.327367;0.000375;0.266531;0.236523]
		Chain 19 - beta = 0.013294	steps: Vector4d[0.401759;0.000460;0.327099;0.290271]
	all chains calibrated.
		Chain  0 - beta = 1.000000 	steps: Vector4d[0.052201;0.000060;0.036125;0.037715]
		Chain  1 - beta = 0.993271 	steps: Vector4d[0.046480;0.000053;0.037842;0.033582]
		Chain  2 - beta = 0.973269 	steps: Vector4d[0.046955;0.000054;0.038229;0.033925]
		Chain  3 - beta = 0.940538 	steps: Vector4d[0.047765;0.000055;0.038889;0.034510]
		Chain  4 - beta = 0.895972 	steps: Vector4d[0.048939;0.000056;0.039844;0.035358]
		Chain  5 - beta = 0.840786 	steps: Vector4d[0.050519;0.000058;0.041131;0.036500]
		Chain  6 - beta = 0.776486 	steps: Vector4d[0.052569;0.000060;0.042800;0.037981]
		Chain  7 - beta = 0.704825 	steps: Vector4d[0.055177;0.000063;0.044923;0.039866]
		Chain  8 - beta = 0.627758 	steps: Vector4d[0.058466;0.000067;0.047601;0.042242]
		Chain  9 - beta = 0.547388 	steps: Vector4d[0.062611;0.000072;0.050976;0.045237]
		Chain 10 - beta = 0.465906 	steps: Vector4d[0.067866;0.000078;0.055254;0.049033]
		Chain 11 - beta = 0.385536 	steps: Vector4d[0.074605;0.000085;0.060741;0.053902]
		Chain 12 - beta = 0.308470 	steps: Vector4d[0.083405;0.000095;0.067906;0.060260]
		Chain 13 - beta = 0.236809 	steps: Vector4d[0.095192;0.000109;0.077502;0.068776]
		Chain 14 - beta = 0.172508 	steps: Vector4d[0.111531;0.000128;0.090805;0.080581]
		Chain 15 - beta = 0.117323 	steps: Vector4d[0.135241;0.000155;0.110109;0.097712]
		Chain 16 - beta = 0.072756 	steps: Vector4d[0.171737;0.000197;0.139823;0.124080]
		Chain 17 - beta = 0.040026 	steps: Vector4d[0.231543;0.000265;0.188514;0.167290]
		Chain 18 - beta = 0.020023 	steps: Vector4d[0.327367;0.000375;0.266531;0.236523]
		Chain 19 - beta = 0.013294 	steps: Vector4d[0.401759;0.000460;0.327099;0.290271]
	calibration summary has been written
	wrote calibration results for 20 chains to calibration_results
	$ 

A more readable output (especially when you used DEBUG) is available in the file "calibration_summary".

---------------------------------
Phase 2 - Running
---------------------------------

If you made it this far, you have almost won! 
You have calibrated chains (with the burn in already done). 

In this phase the program will do the actual sampling, parallel tempering and write out 

#. The visited parameter values of chain0 (beta = 1)
	
	The files are named by the scheme paramname-chain-0.prob.dump.
	These just consist of the visited values for each iteration (doubles for rejects).
	
	These will be used for parameter estimation.
	
#. The probabilities of all chains
	
	The files are named by the scheme prob-chain<chain number>.dump.
	They consist of two columns:
	
	#. posterior probability including prior (as set by the likelihood function
	#. likelihood (excluding prior) as calculated by the likelihood function, but the prior subtracted.
	
	These will be used for the data probability and model selection.

#. "acceptance_rate.dump" allows you to watch the acceptance rates. 

	Its first column is the iteration count, the succeeding columns are the number of accepts.
	
	For convenience, a gnuplot file, acceptance_rate.dump.gnuplot is written that allows you to 
	make a nice plot and press "refresh" in the gnuplot window to watch the progress
	while the program runs (also try "set key left").
	
	On the one hand it would be nice to have the acceptance rates as percentages, but this way we present
	two pieces of information at once: The acceptance rate can be inferred by subtracting 
	the previous row, or estimated by adding 0.5*x to the plot. But it also allows us to see
	when chains get seriously stuck (the plot goes horizontal).

The first two are called "dump files". They can easily reach hundreds of megabytes.
Unless you specify --append, the existing dump files will be overwritten.

The online help is as always available with::

	$ apemost-directory/simplesin.exe help run

and run::

	$ apemost-directory/simplesin.exe run

It is important to realise that the speed of the calculation is *only* limited by the loglikelihood function,
and not by the output written to stdout or the files.

Stopping the run
~~~~~~~~~~~~~~~~~~~~~~~

**Note Bene**: The last line of output files may be invalid. A analysis tool that looks at
the output in real time should ignore it. Read on for why:

For speed purposes, the output to the files is unbuffered. This means the last two lines could be::

	3.592794839126184e-01(newline)
	3.367089(no newline)

And the rest not yet written. This is done efficiently by the operating system, which operates on 
blocks, not on lines.

To force a flush, you can send the USR1 signal to the program:

	$ killall -SIGUSR1 simplesin.exe 

Which will cause the program to flush all files, and then continue to run.

To stop the program, press Ctrl-C or send the TERM signal with kill.
This will also cause a flush, and the files will be cleanly finished.

Unless you specified MAX_ITERATIONS, the program will happily run forever.

You can also pause and continue the program using normal job control (see the manual 
of your shell on how to send STOP and CONT signals).


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Speeding up the run
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can use the benchmark_ program to evaluate the speed of your loglikelihood function.
For example, `pow(a*b, 2)` is faster than `a*a*b*b`. 

You can also get speed improvements from setting N_PARAMETERS. The program will then 
expect the given number of parameters. This allows the compiler to do loop unrolling.



--------------------------------------------

At this point, you are probably waiting for the program to reach a million iterations. 
You deserve a banana (APEmost, get it?).


------------------------------
Phase 3 - Analyse
------------------------------

Since we not only want to fill our hard disks, at some point we will want to 
analyse our data. 

In this phase, all the dump files are read in again. This is often not limited by the CPU, but the 
hard disk speed. As noted in the FAQ, you can analyse your files independently on a different computer,
or paste several dump files together.

So far, APEMoST can produce the following statistics:

#. Marginal distribution histograms

	This gets you the pretty pictures you are looking for, i.e. the full
	posterior probabilities for each parameter. 

	The files are -- appropriately -- named "paramname.histogram".

	NBINS and HISTOGRAMS_MINMAX are flags you might be interested in.

	For convenience, a gnuplot file is written, "marginal_distributions.gnuplot".
	If you remove the leading '#' and run it with gnuplot, it will give you
	a nice graphic of all histograms. For your publication you probably want to use a 
	eps file or a different plot program.

#. MCMC error estimate

	Essentially, this tells you how much the mean of a histogram changes over time. 
	The sigma should be less than 1% of the histogram sigma. (It will say "** high!" 
	if that is not the case.)

	The formula is from `here <http://www.stat.umn.edu/geyer/mcmc/talk/mcmc.pdf>`_.

	You should include this estimate in your publication. 

	This does not do a clean overlapping batch estimate, just analyses a batch of the
	length sqrt(total number of iterations) after another. since the number of iterations
	is high, this should be sufficient (batch length > 500).

#. Model selection / data probability

	This will output the model probability and will let you compare this model to others.

	example output::
	
		Model probability ln(p(D|M, I)): [about 10^-59] -135.52659
		
		Table to compare support against other models (Jeffrey):
		 other model ln(p(D|M,I)) | supporting evidence for this model
		 --------------------------------- 
			>  -135.5 	negative (supports other model)
		  -135.5 .. -145.5 	Barely worth mentioning
		  -145.5 .. -158.5 	Substantial
		  -158.5 .. -169.5 	Strong
		  -169.5 .. -181.5 	Very strong
			<  -181.5 	Decisive
	
	If you have evaluated another model, look up its logarithmic (ln) model probability in this table.


------------------------------------------------

Pretty neat, eh? No cookie now, you got your histograms.

Possible problems in this phase
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#. Straight peaks in the histograms

	These mean a chain got stuck. Bad. 

	Either run until this peak vanishes, change the calibration, ...
	I think you could increase NBINS, and take a average of the neighbouring bins,
	throwing away extreme outliers.

#. The results may be unexpected, or you are not sure if they are correct

	Thinking about it, or simulating the data with the resulting parameters
	may help.

#. Some possible values in the parameter space may have not been detected

	This is one real mean danger, because you probably will never know. 
	A as high number iteration as possible helps.
	
	If two or more peaks have been detected already, you can try to find out
	after how many iterations the last peak showed up. Maybe you should run
	for another so many iterations.
	
	You can try to increase or decrease BETA_0, the beta value of the hottest 
	chain.
	
	You can also try
	to tinker with the calibration or the proposal distribution (e.g. using a 
	distribution with a wider tail such as logit).

	It is a good idea to run the sampling several times and also with 
	different starting points.

#. The heights of different, independent peaks in the histograms do not correctly represent the probability relations.

	This will almost always be the case. Since the runtime is finite, the 
	frequency of visits will be distorted.

	You should evaluate the likelihood function at the peaks to get their real values.
	The eval_ executable and peaks.exe will help you with this.

	peaks.exe will retrieve the median and quartiles of any independent peak in the marginal distribution.
	(independent means 1% of parameter space is unused in between). 
	Since peaks.exe does not use a histogram, it is exact! Prefer it to measuring out the histogram.
	

----------------------------------------
 Hacking APEMoST
----------------------------------------

Feel free to read all the source, write and change algorithms and everything. 

Feedback, ideas, remarks and problems are welcome and will be added to the `FAQ <faq.html>`_.

As the `license <license.html>`_ states, since we worked so hard on APEMoST and you get it for free,
you are expected to contribute changes back to us, so everyone can profit. 

Ideally, get familiar with git, which is the version control system in use.
Some resources are here:

- http://cworth.org/hgbook-git/tour/
- http://git-scm.com/ http://book.git-scm.com/1_welcome_to_git.html
- http://zrusin.blogspot.com/2007/09/git-cheat-sheet.html

The most important commands are "git pull", "git commit" and "git format-patch". 
The last allows you to send us a patch of your changes, so everyone can profit from it.

You can also set up your own repository (which is very easy, e.g. on github), and just tell 
me that you will contribute there. This will allow me to pull your changes.

**If this is all too much for you** -- before you decide not to contribute back -- a
tarball or zip file is also welcome. The contact address can be found at the 
`contact page <contact.html>`_.

That said, a version control system is really useful to stay on top of things (e.g. trying out 
some code). Consider using it for your other projects. 
If you don't like git, try hg, which has better GUIs. There is also a hg-git bridge. 

-------------------------------------
Other topics
-------------------------------------
Random generator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
(Pseudo-) random number generation is a very important topic and should be
addressed. We use the default random generator from GSL. This can be influenced
with environment variables, for example setting GSL_RANDOM_SEED and GSL_RNG_TYPE.
See the GSL manual.

Only one random generator is used for the whole program, so setting the 
seed will not result in multiple, synchronized random generators.

Set a different seed for different runs, otherwise you will always obtain the
same results! 

For example, you can use a random number as the initial seed. If you use bash::

	export GSL_RANDOM_SEED=$RANDOM

Mention in your publication that you set or varied the seed. Otherwise you may
be victim to systematic errors!

---------------------------------------
Concluding remarks
---------------------------------------

None.








.. include:: footer.txt

