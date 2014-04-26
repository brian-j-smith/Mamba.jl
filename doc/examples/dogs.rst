.. index:: Examples; Dogs: Loglinear Model for Binary Data

Dogs: Loglinear Model for Binary Data
-------------------------------------

An example from Kalbfleisch :cite:`kalbfleisch:1985:PSI` concerning the Solomon-Wynne experiment on dogs.  In the experiment, 30 dogs were subjected to 25 trials.  On each trial, a barrier was raised, and an electric shock was administered 10 seconds later if the dog did not jump the barrier.

Model
^^^^^

Failure to jump the barrier in time is modelled as

.. math::

	y_{i,j} &= \text{Bernoulli}(\pi_{i,j}) \quad\quad i=1,\ldots,30; j=2,\ldots,25 \\
	\log(\pi_{i,j}) &= \alpha x_{i,j} + \beta (j - x_{i,j}) \\
	\alpha &\sim \text{Uniform}(-1000, 0) \\
	\beta &\sim \text{Uniform}(-1000, 0),
	
where :math:`y_{i,j} = 1` if dog :math:`i` does not jump the barrier before the shock (failure) on trial :math:`j`, and 0 otherwise; :math:`x_{i,j}` is the number of successes prior to trial :math:`j`; and :math:`\pi_{i,j}` is the probability of failure.

.. literalinclude:: dogs.jl
	:language: julia
