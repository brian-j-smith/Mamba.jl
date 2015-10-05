.. index:: Examples; Line: Block-Specific Sampling with AMWG and Slice

.. _example-Line_AMWG_Slice:

Line: Block-Specific Sampling with AMWG and Slice
-------------------------------------------------

An example based on the linear regression model defined in the :ref:`section-Line` section.  The program below illustrates use of the stand-alone :func:`amwg!` and :func:`slice!` functions to sample different parameter blocks within the same MCMC algorithm.

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: line_amwg_slice.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 1:10000
    Thinning interval = 1
    Chains = 1
    Samples per chain = 10000

    Empirical Posterior Estimates:
          Mean        SD       Naive SE       MCSE       ESS
    b0 0.64401798 0.99315634 0.0099315634 0.060725564  267.4805
    b1 0.78985612 0.29790444 0.0029790444 0.017888106  277.3481
    s2 1.20785292 2.96033511 0.0296033511 0.062566344 2238.7222

    Quantiles:
           2.5%       25.0%      50.0%     75.0%     97.5%
    b0 -1.33127385 0.075527035 0.6403226 1.1902679 2.7665517
    b1  0.16055439 0.625740352 0.7923275 0.9527029 1.3903861
    s2  0.16673189 0.381185645 0.6538295 1.2373814 5.6065938
