.. index:: Examples; Line: Approximate Bayesian Computation

.. _example-Line_ABC:

Line: Approximate Bayesian Computation
--------------------------------------

A simple example to demonstrate the Approximate Bayesian Computation (ABC) sampler within the MCMC framework, based on the linear regression model defined in the :ref:`section-Line` section.  ABC sampling is applied separately to the ``:beta`` and ``:s2`` parameter blocks.  Different summary statistics are specified to show a range of functions that could be used.  More common practice is to use the same data summaries for all ABC-sampled parameters.

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: line_abc.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 1001:10000
    Thinning interval = 1
    Chains = 1,2,3
    Samples per chain = 9000

    Empirical Posterior Estimates:
               Mean        SD       Naive SE       MCSE       ESS
         s2 1.30743333 1.99877929 0.0121641834 0.083739029 569.73624
    beta[1] 0.72349922 1.03842764 0.0063196694 0.039413390 694.16848
    beta[2] 0.77469344 0.31702542 0.0019293553 0.011392989 774.30630

    Quantiles:
                2.5%        25.0%      50.0%      75.0%     97.5%
         s2  0.048095084 0.23351203 0.57947788 1.45858829 7.7639321
    beta[1] -1.309713807 0.12616636 0.67263204 1.27579373 3.1735176
    beta[2]  0.107216316 0.59367961 0.77867235 0.95156086 1.4043715
