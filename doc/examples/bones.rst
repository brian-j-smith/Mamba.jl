.. index:: Examples; Bones: Latent Trait Model for Multiple Ordered Categorical Responses

.. _example-Bones:

Bones: Latent Trait Model for Multiple Ordered Categorical Responses
--------------------------------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex`, Roche *et al.* :cite:`roche:1975:SM`, and Thissen :cite:`thissen:1986:MUG` concerning skeletal age in 13 boys predicted from 34 radiograph indicators of skeletal maturity.

Model
^^^^^

Skeletal ages are modelled as

.. math::

    \operatorname{logit}(Q_{i,j,k}) &= \delta_j (\theta_i - \gamma_{j,k}) \quad\quad i=1,\ldots,13; j=1,\ldots,34; k=1,\ldots,4 \\
    \theta_i &\sim \text{Normal}(0, 100),

where :math:`\delta_j` is a discriminability parameter for indicator :math:`j`, :math:`\gamma_{j,k}` is a threshold parameter, and :math:`Q_{i,j,k}` is the cumulative probability that boy :math:`i` with skeletal age :math:`\theta_i` is assigned a more mature grade than :math:`k`.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: bones.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 2502:10000
    Thinning interval = 2
    Chains = 1,2
    Samples per chain = 3750

    Empirical Posterior Estimates:
                  Mean        SD       Naive SE       MCSE         ESS
     theta[1]  0.32603385 0.20640874 0.0023834028 0.0039448110 2737.81276
     theta[2]  1.37861692 0.25824308 0.0029819342 0.0058663965 1937.82503
     theta[3]  2.35227822 0.27998526 0.0032329913 0.0067161153 1737.93707
     theta[4]  2.90165730 0.29713320 0.0034309987 0.0078730921 1424.33353
     theta[5]  5.54427283 0.50242324 0.0058014839 0.0169090038  882.88350
     theta[6]  6.70804782 0.57206890 0.0066056827 0.0221532973  666.83738
     theta[7]  6.49138381 0.60154625 0.0069460578 0.0219158412  753.39330
     theta[8]  8.93701249 0.73636136 0.0085027686 0.0336199950  479.71875
     theta[9]  9.03585289 0.65172497 0.0075254717 0.0233182299  781.15561
    theta[10] 11.93125529 0.69360918 0.0080091090 0.0282955741  600.88678
    theta[11] 11.53686992 0.92271657 0.0106546132 0.0493587234  349.46912
    theta[12] 15.81482824 0.54261736 0.0062656056 0.0210976666  661.48275
    theta[13] 16.93028146 0.72458739 0.0083668145 0.0323302069  502.30161

    Quantiles:
                  2.5%       25.0%       50.0%       75.0%       97.5%
     theta[1] -0.11215555  0.19557824  0.33881555  0.45840506  0.7174563
     theta[2]  0.91705346  1.19969433  1.36116575  1.53751273  1.9466119
     theta[3]  1.78287586  2.17136780  2.35623189  2.53035766  2.9211580
     theta[4]  2.32940825  2.69621746  2.89121336  3.10758151  3.4945343
     theta[5]  4.59142954  5.20543314  5.53392246  5.86525435  6.5867245
     theta[6]  5.56649066  6.30983797  6.70666338  7.09569168  7.8229872
     theta[7]  5.38663728  6.07628064  6.46533033  6.88636840  7.7051374
     theta[8]  7.47304526  8.43125608  8.96072241  9.45344704 10.2856733
     theta[9]  7.80477915  8.60559136  9.01498109  9.46962522 10.3024722
    theta[10] 10.64129157 11.48379528 11.89611699 12.37737647 13.3873043
    theta[11]  9.83558611 10.88717498 11.49029895 12.15757004 13.4263451
    theta[12] 14.79250437 15.45889470 15.79840132 16.15824313 16.9593310
    theta[13] 15.61843069 16.42289725 16.90719268 17.41900248 18.3895761
