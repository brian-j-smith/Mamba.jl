.. index:: Examples; Bones: Latent Trait Model for Multiple Ordered Categorical Responses

.. _example-Bones:

Bones: Latent Trait Model for Multiple Ordered Categorical Responses
--------------------------------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1`, Roche *et al.* :cite:`roche:1975:SM`, and Thissen :cite:`thissen:1986:MUG` concerning skeletal age in 13 boys predicted from 34 radiograph indicators of skeletal maturity. 

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
    14x6 Array{Any,2}:
     ""             "Mean"    "SD"      "Naive SE"   "MCSE"         "ESS"
     "theta[1]"    0.326034  0.206409  0.0023834    0.00556249  1376.95
     "theta[2]"    1.37862   0.258243  0.00298193   0.00769705  1125.66
     "theta[3]"    2.35228   0.279985  0.00323299   0.00886942   996.507
     "theta[4]"    2.90166   0.297133  0.003431     0.0097304    932.484
     "theta[5]"    5.54427   0.502423  0.00580148   0.0229152    480.72
     "theta[6]"    6.70805   0.572069  0.00660568   0.0292519    382.461
     "theta[7]"    6.49138   0.601546  0.00694606   0.0303562    392.683
     "theta[8]"    8.93701   0.736361  0.00850277   0.0427413    296.815
     "theta[9]"    9.03585   0.651725  0.00752547   0.0312046    436.207
     "theta[10]"  11.9313    0.693609  0.00800911   0.0390483    315.518
     "theta[11]"  11.5369    0.922717  0.0106546    0.0655843    197.942
     "theta[12]"  15.8148    0.542617  0.00626561   0.0285355    361.59
     "theta[13]"  16.9303    0.724587  0.00836681   0.0433486    279.403

    Quantiles:
    14x6 Array{Any,2}:
     ""             "2.5%"     "25.0%"    "50.0%"    "75.0%"    "97.5%"
     "theta[1]"   -0.112156   0.195578   0.338816   0.458405   0.717456
     "theta[2]"    0.917053   1.19969    1.36117    1.53751    1.94661
     "theta[3]"    1.78288    2.17137    2.35623    2.53036    2.92116
     "theta[4]"    2.32941    2.69622    2.89121    3.10758    3.49453
     "theta[5]"    4.59143    5.20543    5.53392    5.86525    6.58672
     "theta[6]"    5.56649    6.30984    6.70666    7.09569    7.82299
     "theta[7]"    5.38664    6.07628    6.46533    6.88637    7.70514
     "theta[8]"    7.47305    8.43126    8.96072    9.45345   10.2857
     "theta[9]"    7.80478    8.60559    9.01498    9.46963   10.3025
     "theta[10]"  10.6413    11.4838    11.8961    12.3774    13.3873
     "theta[11]"   9.83559   10.8872    11.4903    12.1576    13.4263
     "theta[12]"  14.7925    15.4589    15.7984    16.1582    16.9593
     "theta[13]"  15.6184    16.4229    16.9072    17.419     18.3896
