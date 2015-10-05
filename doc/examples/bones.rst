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
                  Mean        SD      Naive SE       MCSE        ESS
     theta[1]  0.32603385 0.2064087 0.0023834028 0.005562488 1376.94938
     theta[2]  1.37861692 0.2582431 0.0029819342 0.007697049 1125.66422
     theta[3]  2.35227822 0.2799853 0.0032329913 0.008869417  996.50667
     theta[4]  2.90165730 0.2971332 0.0034309987 0.009730398  932.48357
     theta[5]  5.54427283 0.5024232 0.0058014839 0.022915189  480.72039
     theta[6]  6.70804782 0.5720689 0.0066056827 0.029251931  382.46138
     theta[7]  6.49138381 0.6015462 0.0069460578 0.030356237  392.68306
     theta[8]  8.93701249 0.7363614 0.0085027686 0.042741328  296.81508
     theta[9]  9.03585289 0.6517250 0.0075254717 0.031204552  436.20719
    theta[10] 11.93125529 0.6936092 0.0080091090 0.039048344  315.51820
    theta[11] 11.53686992 0.9227166 0.0106546132 0.065584299  197.94151
    theta[12] 15.81482824 0.5426174 0.0062656056 0.028535483  361.59041
    theta[13] 16.93028146 0.7245874 0.0083668145 0.043348628  279.40285

    Quantiles:
                 2.5%       25.0%      50.0%       75.0%        97.5%
     theta[1] -0.1121555  0.1955782  0.33881555  0.45840506  0.717456283
     theta[2]  0.9170535  1.1996943  1.36116575  1.53751273  1.946611947
     theta[3]  1.7828759  2.1713678  2.35623189  2.53035766  2.921158047
     theta[4]  2.3294082  2.6962175  2.89121336  3.10758151  3.494534287
     theta[5]  4.5914295  5.2054331  5.53392246  5.86525435  6.586724493
     theta[6]  5.5664907  6.3098380  6.70666338  7.09569168  7.822987248
     theta[7]  5.3866373  6.0762806  6.46533033  6.88636840  7.705137414
     theta[8]  7.4730453  8.4312561  8.96072241  9.45344704 10.285673300
     theta[9]  7.8047792  8.6055914  9.01498109  9.46962522 10.302472161
    theta[10] 10.6412916 11.4837953 11.89611699 12.37737647 13.387304288
    theta[11]  9.8355861 10.8871750 11.49029895 12.15757004 13.426345115
    theta[12] 14.7925044 15.4588947 15.79840132 16.15824313 16.959330990
    theta[13] 15.6184307 16.4228972 16.90719268 17.41900248 18.389576101
