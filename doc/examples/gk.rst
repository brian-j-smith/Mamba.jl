.. index:: Examples; GK: Approximate Bayesian Computation

.. _example-gk:

GK: Approximate Bayesian Computation
--------------------------------------

Approximate Bayesian Computation (ABC) is often useful when the likelihood is costly to compute. The generalized GK distribution :cite:`rayner:2002:GK`  is a distribution defined by the inverse of its cumulative distribution function. It is therefore easy to sample from, but difficult to obtain analytical likeihood functions. These properties make the GK distribution well suited for ABC. The following is a simulation study that mirrors `allingham:2009:ABCGK`.

Model
^^^^^
.. math::

    x_{i} &\sim \text{GK}\left(A, B, g, k\right) \quad\quad i=1,\ldots,10000 \\
    A &\sim \text{Uniform}(0, 10) \\
    B &\sim \text{Uniform}(0, 10) \\
    g &\sim \text{Uniform}(0, 10) \\
    k &\sim \text{Uniform}(0, 10)

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: gk.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia
    Iterations = 4001:15000
    Thinning interval = 1
    Chains = 1,2,3
    Samples per chain = 11000

    Empirical Posterior Estimates:
          Mean         SD        Naive SE       MCSE        ESS   
     k 0.48221758 0.061956076 0.00034105698 0.0021518021  829.0165
     A 2.87330714 0.114382460 0.00062965472 0.0023599697 2349.1243
     B 0.99777598 0.104611155 0.00057586545 0.0035817717  853.0219
     g 5.65630699 2.580219260 0.01420363956 0.0791397180 1062.9779

     Quantiles:
          2.5%      25.0%      50.0%     75.0%     97.5%  
     k 0.35557416 0.4420227 0.48229997 0.5235059 0.6029042
     A 2.66385620 2.7956063 2.86717760 2.9436771 3.1158983
     B 0.80773524 0.9256279 0.99225590 1.0625602 1.2221914
     g 1.53996630 3.3474209 5.62664473 7.9560817 9.8046807
