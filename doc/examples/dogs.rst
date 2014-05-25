.. index:: Examples; Dogs: Loglinear Model for Binary Data

.. _example-Dogs:

Dogs: Loglinear Model for Binary Data
-------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1`, Lindley and Smith :cite:`lindley:1972:BEL`, and Kalbfleisch :cite:`kalbfleisch:1985:PSI` concerning the Solomon-Wynne experiment on dogs.  In the experiment, 30 dogs were subjected to 25 trials.  On each trial, a barrier was raised, and an electric shock was administered 10 seconds later if the dog did not jump the barrier.

Model
^^^^^

Failures to jump the barriers in time are modelled as

.. math::

	y_{i,j} &= \text{Bernoulli}(\pi_{i,j}) \quad\quad i=1,\ldots,30; j=2,\ldots,25 \\
	\log(\pi_{i,j}) &= \alpha x_{i,j-1} + \beta (j - 1 - x_{i,j-1}) \\
	\alpha, \beta &\sim \text{Flat}(-\infty, -0.00001) \\
	
where :math:`y_{i,j} = 1` if dog :math:`i` fails to jump the barrier before the shock on trial :math:`j`, and 0 otherwise; :math:`x_{i,j-1}` is the number of successful jumps prior to trial :math:`j`; and :math:`\pi_{i,j}` is the probability of failure.

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: dogs.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia
