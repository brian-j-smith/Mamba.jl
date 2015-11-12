.. index:: Examples; Magnesium: Meta-Analysis Prior Sensitivity

.. _example-Magnesium:

Magnesium: Meta-Analysis Prior Sensitivity
------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex`.

Model
^^^^^

Number of events reported for treatment and control subjects in 8 studies is modelled as

.. math::

    r^c_j &\sim \text{Binomial}(n^c_j, p^c_{i,j}) \quad\quad i=1,\ldots,6; j=1,\ldots,8 \\
    p^c_{i,j} &\sim \text{Uniform}(0, 1) \\
    r^t_j &\sim \text{Binomial}(n^t_j, p^t_{i,j}) \\
    \operatorname{logit}(p^t_{i,j}) &= \theta_{i,j} + \operatorname{logit}(p^c_{i,j}) \\
    \theta_{i,j} &\sim \text{Normal}(\mu_i, \tau_i) \\
    \mu_i &\sim \text{Uniform}(-10, 10) \\
    \tau_i &\sim \text{Different Priors},

where :math:`r^c_j` is the number of control group events, out of :math:`n^c_j`, in study :math:`j`; :math:`r^t_j` is the number of treatment group events; and :math:`i` indexes differ prior specifications.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: magnesium.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 2502:12500
    Thinning interval = 2
    Chains = 1,2
    Samples per chain = 5000

    Empirical Posterior Estimates:
              Mean        SD       Naive SE       MCSE         ESS
    tau[1] 0.55516794 0.39604206 0.0039604206 0.0249542854  251.87922
    tau[2] 1.10064876 0.59493159 0.0059493159 0.0259756957  524.56544
    tau[3] 0.81391916 0.50206143 0.0050206143 0.0222820890  507.69387
    tau[4] 0.46761803 0.27315180 0.0027315180 0.0152718027  319.90979
    tau[5] 0.50273023 0.35215203 0.0035215203 0.0199145134  312.69505
    tau[6] 0.56535848 0.19240620 0.0019240620 0.0068362546  792.13952
     OR[1] 0.45972229 0.14882861 0.0014882861 0.0068020036  478.73933
     OR[2] 0.42192295 0.24621559 0.0024621559 0.0066694219 1362.87094
     OR[3] 0.44135115 0.21012901 0.0021012901 0.0068759175  933.92227
     OR[4] 0.47081945 0.14219883 0.0014219883 0.0071945264  390.64984
     OR[5] 0.47900077 0.14791826 0.0014791826 0.0074993502  389.04186
     OR[6] 0.44708748 0.14007210 0.0014007210 0.0054715042  655.37495

    Quantiles:
               2.5%       25.0%      50.0%      75.0%      97.5%
    tau[1] 0.053249224 0.26393250 0.47726923 0.74426901 1.58474129
    tau[2] 0.323712679 0.68508835 0.97547583 1.37889339 2.63318913
    tau[3] 0.140369113 0.48170109 0.72504452 1.03577966 2.03904609
    tau[4] 0.093345642 0.27431059 0.41201604 0.59915029 1.13215956
    tau[5] 0.036043879 0.24783727 0.43547026 0.68157942 1.36172745
    tau[6] 0.184141434 0.43413529 0.56424751 0.69657737 0.94066205
     OR[1] 0.184900548 0.36280207 0.45956086 0.55354433 0.76626505
     OR[2] 0.119614981 0.27183699 0.38455909 0.51625616 0.95025887
     OR[3] 0.152008568 0.31438014 0.42065000 0.53999576 0.86382432
     OR[4] 0.196690338 0.37280299 0.46795381 0.56729578 0.74612951
     OR[5] 0.202086157 0.38165881 0.47934537 0.57673386 0.76497223
     OR[6] 0.214535135 0.34754950 0.43349452 0.52942101 0.76370767
