FAQ - Frequently Asked Questions
===========================================

.. contents::

- *Is APEMoST free? How is APEMoST licensed? What are my obligations?*

	See the `license <license.html>`_.

- *What do I have to do in my publication?*

	State the version of APEMoST, and that you used it (see next question). 

	Also double-check that your results make sense.
	
	Include error estimates in your publication.
	
	Mention that you set the PRNG seed.

	If you made modifications to APEMoST, `contribute them back <contact.html>`_ or publish them otherwise. 
	If you don't, others can not check your results and uncover errors.

	You can publish your likelihood function, parameter and data file in the appendix of your paper.
	If you don't, others may have a hard time to reconstruct your results.
	You are also encouraged to send them to us. We will publish them on the sourceforge website 
	so you can just put a link (that always stays the same) in your paper. This also allows us to show how 
	APEMoST is used by others.

- *How do I cite APEMoST correctly?*

	A publication is in the works. Please write us a email, so we can notify you of the correct
	citation. In the meantime, you can use the website as a placeholder.

- *Can I run APEMoST on a Cluster/Grid?*

	Since Parallel Tempering requires the chains to intercompare and swap their status frequently,
	I estimated that a distributed computation would actually be slowed down by the synchronization overhead.
	This of course depends on how long it takes to evaluate your likelihood function, but with 300 data tuples
	and a 6-parametric model, we can have 1000 iterations per second (in each chain) on a single-core machine.
	With comparing chains every 100 iterations or so you can see that waiting for synchronisation over the 
	network is problematic. 

	Ideally, get a fast computer with many CPUs.

	You can also run the program independently on several machines, then paste the output (visited values)
	together, and analyze the combined runs.

	Other software packages use MPI for running on clusters.

- *Can I run APEMoST on Linux/Unix/Mac?*

	Yes, all Unix derivates are the primary target of APEMoST. 

- *Can I run APEMoST on Windows?*

	We don't test APEMoST on Windows, but there is no reason you can't. You need the software
	APEMoST is based on, and the gcc compiler. 

	Alternatively, if you run into trouble, you can try working in a cygwin environment
	or run a Linux in a virtual machine.

- *How can I stay up to date?*

	The easiest to get notified is to `contact us <contact.html>`_.

	You can follow the `news on our project page <http://sourceforge.net/projects/apemost/>`_
	and also the `code changes <http://apemost.git.sourceforge.net/git/gitweb.cgi?p=apemost/apemost;a=summary>`_.
	

- *My question has not been answered.*

	If you have a technical question, maybe the answer is in the manual?

	Otherwise feel free to `contact us <contact.html>`_.


