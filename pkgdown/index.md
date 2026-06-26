<div style="display: flex; align-items: center; gap: 10px; margin-bottom: 0.75rem;">
  <img src="reference/figures/logo.png"
       alt="bfbin2arm logo"
       style="height: 84px; width: auto;">
  <h2 style="margin: 0;">bfbin2arm</h2>
</div>

## Sequential Bayesian trial design for clinical phase II trials

Sequential hypothesis tests are an important tool to improve the efficiency of clinical trials.
In contrast to study designs with a fixed sample size, interim analyses are carried out which allow to stop a trial early for futility when a novel drug or treatment is ineffective.
Often, such designs are applied in phase II proof of concept trials, where the primary endpoint measures the binary response (success or failure) of each enrolled patient to the novel drug or treatment.
The R software package `bfbin2arm` provides methodology and software for planning, design and 
calibration of such trials, and allows to analyze the operating-characteristic (such as the 
power and type-I-error) of sequential Bayesian phase II clinical trial with binary endpoint. 
The package focuses on single-arm and two-arm Bayesian designs, including one-stage and 
two-stage settings, with tools for design calibration and visualisation of trial operating 
behavior. Importantly, the package allows for different calibration modes for a trial design, 
including Bayesian, frequentist and hybrid calibration.

The package website collects methodological articles, worked examples, and reference documentation intended to support both applied use and methodological development.

## Funding

Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - Project number 549296018.

For any inquiries, contact Dr. Riko Kelter, [riko.kelter@uni-koeln.de](mailto:riko.kelter@uni-koeln.de).


<div style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 24px; align-items: center; justify-items: center; margin-top: 1em; margin-bottom: 1.5em;">

<div style="grid-column: 1 / span 3; text-align: center;">
<img src="man/figures/imsb_logo.png" alt="Institute of Medical Statistics and Computational Biology logo" style="max-height: 110px; width: auto;">
</div>

<div style="text-align: center;">
<img src="man/figures/ukoeln_logo.png" alt="University of Cologne logo" style="max-height: 90px; width: auto;">
</div>

<div style="text-align: center;">
<img src="man/figures/ukkoeln_logo.png" alt="University Hospital Cologne logo" style="max-height: 90px; width: auto;">
</div>

<div style="text-align: center;">
<img src="man/figures/dfg_logo.jpg" alt="Deutsche Forschungsgemeinschaft logo" style="max-height: 90px; width: auto;">
</div>

</div>

## Main articles

The software package `bfbin2arm` emerged from clinical trial designs for two-arm clinical phase II trials with binary endpoints, using Bayes factors as the primary measure of evidence.
Over time, more designs were implemented, including one- and two-stage designs, single-arm and two-arm designs, designs for equivalence, non-inferiority and superiority tests, optimal
designs which minimize the expected sample size under the null hypothesis, as well as different calibration modes for the available designs such as Bayesian, frequentist or hybrid calibration.

<div style="text-align: center; margin: 1.5em 0;">
  <img src="man/figures/bfbin2arm_overview.png" 
       alt="Diagram showing that bfbin2arm includes single-arm and two-arm, one-stage and two-stage designs, multiple calibration modes, and directional, two-sided, and ROPE equivalence tests." 
       style="max-width: 100%; width: 700px;">
</div>

An overview about the implemented designs is given in:

- [Clinical trial designs in bfbin2arm](articles/bfbin2arm-overview.html)

Equivalence testing in single-arm one-stage designs:

- [ROPE-based trial design for single-arm one-stage phase II trials with binary endpoints](articles/bfbin2arm-rope-singlearm-onestage-design.html)
- [Frequentist and hybrid calibration of one-stage ROPE-based designs for single-arm phase II trials](articles/bfbin2arm-rope-singlearm-onestage-calibration.html)

One-stage single-arm designs:

- [Calibration of Bayesian one-stage designs for single-arm phase II trials with binary endpoints](articles/bfbin2arm-singlearm-onestage.html)

Two-stage (optimal) single-arm designs:

- [Optimal Bayesian calibration for single-arm two-stage Bayes factor designs](articles/bfbin2arm-singlearm-twostage_bayesian.html)
- [Optimal frequentist calibration for single-arm two-stage Bayes factor designs with binary endpoints](articles/bfbin2arm-singlearm_twostage_frequentist.html)
- [Optimal hybrid calibration for single-arm two-stage Bayes factor designs with binary endpoints](articles/bfbin2arm-singlearm_twostage_hybrid.html)
- [Optimal full calibration for single-arm two-stage Bayes factor designs with binary endpoints](articles/bfbin2arm-singlearm_twostage_full.html)

One-stage two-arm designs:

- [Bayesian calibration of two-arm one-stage Bayes factor designs with binary endpoints](articles/bfbin2arm-twoarm_onestage_Bayesian.html)

Two-stage (optimal) two-arm designs:

- [Optimal Bayesian calibration of two-arm two-stage Bayes factor designs with binary endpoints](articles/bfbin2arm-twoarm-twostage_Bayesian.html)

## Reference documentation

- [Function reference](reference/index.html)

## Source code

- [GitHub repository](https://github.com/rikokelter/bfbin2arm)