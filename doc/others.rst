Links and other software
=================================


.. contents::
Introduction to MCMC and Bayesian inference
-----------------------------------------------

If you are unfamiliar with the concept, it can be quite overwhelming what it is all about, 
especially since there are many incomplete introductions. Here, we try to list several sources we consider good
that cover the subject, so you can intercompare and reread in one what remains unclear in the other.

In one paragraph? Lets try. This is just to get a grasp on how the terms relate to each other.

	*You have a model that has several parameters, which should predict or describe an observation. 
	An example of a model with three parameters (A, B, F) is "y(t) = A*sin(F*t + B)". It describes the 
	relation between y and t.*

	*Assume now you can specify how likely it is that the real observation data has been observed under this model.
	This is called the likelihood function p(D|M, I). Now we would need to evaluate this function everywhere
	in parameter space, which is simply unfeasible. We would need a method that evaluates there more densely,
	where the value is high. This is where Markov Chain Monte Carlo comes in (Metropolis algorithm, proposal distribution).
	With MCMC and the Bayesian theorem, you can not only tell what the best fitting values are, you also get
	the probability of each possible parameter value as a marginal distribution. And, you can compare the 
	likelihood of one model to another (in the example above, we could add "+ C").*

	*Parallel Tempering, calibrations, different proposal distributions, adaptive MCMC are enhancements to improve convergence.*

Well, not quite one paragraph.

#. A very good book on the subject:
	*Bayesian logical data analysis for the physical sciences* by P.C. Gregory
	You can find it in `Google Books <http://books.google.com/books?id=yJ_5VFo0zGMC>`_ and at the `Cambridge Catalogue <http://www.cambridge.org/catalogue/catalogue.asp?isbn=052184150X>`_.

#. Then there is always place to promote `my bachelor thesis <http://textfeld.ac.at/text/1589>`_ :-)

	I think it provides a pretty good overview of everything necessary. As a bonus, 
	we detect pulsations on a star.
	
	The tool has a different name in there, and model selection is not covered (which 
	works very nicely on the pulsation example).
	
	It is based on the `work of Michael Gruberbauer <http://arxiv.org/abs/0811.3345>`_.

#. Then there is always Wikipedia, although I found it hard to understand for novices.

	http://en.wikipedia.org/wiki/MCMC http://en.wikipedia.org/wiki/Bayesian_inference

#. Geyer has a `good presentation <http://www.stat.umn.edu/geyer/mcmc/talk/mcmc.pdf>`_ on how MCMC works.

#. There are some astronomy-related papers

	`Bayes in the sky: Bayesian inference  and model selection in cosmology <http://arxiv.org/abs/0803.4089>`_
	provides a good introduction and highlights cosmology as a good application for Bayesian inference with MCMC.

	and `Applications of Bayesian model selection to cosmological parameters <http://arxiv.org/abs/astro-ph/0504022>`_

#. ...
	
	If you know or have any good resources on the topic, you are welcome to :doc:`suggest links <contact>`.

Differences and other programs
--------------------------------

Obviously, APEMoST is not the first program to do Bayesian inference or MCMC sampling. 

~~~~~~~~~~~~~
 Differences
~~~~~~~~~~~~~

We think it differs from existing programs by

- specifying and programming the model in a imperative language

	Although a declarative specification of the likelihood function has some
	beauty, we find a imperative specification to be easier to access, debug and
	faster to evaluate.

	The downside is that complex likelihood relations may not be expressible with APEMoST.

- not providing a GUI (so far)

	We provide a lean and mean command based structure. We think it is easier to 
	understand and debug than other tools.

	The output (text files with double values) can be analysed with any plot program 
	or any program written in any programming language. 
	We do provide tools to analyse the output (e.g. creating histograms, etc.). 

- a simple workflow

	The new user is not overwhelmed by dozens of features when [s]he does not have to be.

- Metropolis sampling
   
	Other programs usually use Gibbs Sampling, we use classic Metropolis sampling (usually without 
	needing the Hastings extension, since the proposal distributions are normally symmetric). 

- general purpose for data analysis

	We are no ''pure statistician'', but we also don't just target Astronomy and Physics. All natural
	sciences are invited to use and extend our program.

- speed

	If another program converges twice as fast due to some sophisticated, experimental concept,
	but APEMoST runs 20 times as many iterations in the same time ...
	
	Regarding speed we can say that APEMoST, running 2 million iterations,
	spends less than a minute outside the likelihood function. So, the 
	bottleneck will definitely not be APEMoST, and C allows the most 
	efficient implementation of your likelihood function.

~~~~~~~~~~~~~~~~~
 Other programs
~~~~~~~~~~~~~~~~~

You are welcome to notify us if you know others!

- BUGS_/WinBUGS/OpenBUGS_/LinBUGS

	GUI-based; used by statisticians. Uses Gibbs sampling.

	The BUGS language is in wide use as declarative description of the model.

.. _OpenBUGS: http://mathstat.helsinki.fi/openbugs/
.. _BUGS: http://www.mrc-bsu.cam.ac.uk/bugs/welcome.shtml

- JAGS http://calvin.iarc.fr/~martyn/software/jags/

	Similar to BUGS. Manual says it has poor performance.

- MCMC R package http://www.stat.umn.edu/geyer/mcmc/

	There is a mcmc package for the R project.
	Has very good slides of an introduction to MCMC sampling http://www.stat.umn.edu/geyer/mcmc/talk/mcmc.pdf

- BioBayes http://www.dcs.gla.ac.uk/biobayes/

	A Software Package for Bayesian Inference in Systems Biology

	Has a very nice, user-friendly GUI (Java).
	Has a very nice video presentation

- COSMOMC http://cosmologist.info/cosmomc/ http://cosmologist.info/notes/

	A mcmc program specialized to cosmological problems. Also uses the GSL.
	Has a nice `presentation <http://cosmologist.info/notes/MCMC.ppt>`_ on what it is about.

- FBM Software for Flexible Bayesian Modeling http://www.cs.utoronto.ca/~radford/fbm.software.html

	ANSI-C
	
	“Flexible Bayesian models for regression and classifica- 
	tion based on neural networks and Gaussian processes, 
	and for probability density estimation using mixtures. 
	Neural net training using early stopping is also sup-
	ported.”
	“Markov chain Monte Carlo methods, and their appli-
	cations to Bayesian modeling, including implementations 
	of Metropolis, hybrid Monte Carlo, slice sampling, and      
	tempering methods. “
	
	This looks like the most similar approach (being ANSI-C). Looks very powerful and complete.

- less relevant software follows

- Bassist http://www.cs.helsinki.fi/research/fdk/bassist/

	generates C++ code: Bayesian model
	data -> posterior distribution of model parameters

- mrbayes http://mrbayes.csit.fsu.edu/

	bayesian inference with biology models (several discrete options)

- BEAST http://www.beastsoftware.org/
	
	modeling population models. Java

- YADAS http://www.ccs.lanl.gov/ccs6/yadas/

- There is a MCMC sampler written in Python, `pyMC <http://code.google.com/p/pymc/>`_

- Programs that are not available for download are not listed (e.g. "bayesiananalysis", and the one from Do Kester)

Some notes: several software packages are abandoned since a few years. 

I (Johannes) find it hard to get into the programs and to understand them. A toy example that works both in e.g. BUGS/JAGS, BioBayes, fbm and APEMoST would be great.


