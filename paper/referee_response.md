# Referee response

We thank the referee for their precious work and for the comment they made.
We read carefully throught the proposed changes and we amended the text accordingly. The referee can find below thes answers to their questions/concerns

## PnP comments

- Please add citations to the LIGO and Virgo instrument papers the first time these detectors are mentioned (1411.4547, 1408.3978).
	+ Done: page 9

- Please fix the broken citation at the beginning of section IV to "the analytical PSD computed for LIGO Handford interferometer, released together with the GWCT-1 catalog”
	+ Done, added also the link to DCC entry where the PSD are taken from

- Please link to the correct citation at the start of Section VB, as [27] links to the climate data. "We focus on the public data released by the LIGO/Virgo collaboration [27]. We reconstruct the PSD both assuming the CAT and FPE loss functions on 1000s of data [27] from the Livingston observatory starting from GPS time 1164603392. The data are sampled at 4096 Hz.”
	+ Done

- "Moreover, if a glitch is detected, its shape can be estimated (as well as the confidence level) by subtracting the expected signal with the actual signal.” Here the Bayeswave papers should be cited (1410.3835, 1410.3852) as well as other key LIGO papers describing glitches, their properties, and subtraction (1508.07316, 1602.03844, 1808.03619), and I believe [28] should point to the Zevin et. al paper in [29].
	+ Done: Bayeswave (as well as a study of its accuracy) is cited in sec.IV, page 10, together with 1907.06540 that assess its accuracy; papers on glitches are cited on sec. 5.2

- "for short time series, comparable with the length of binary black hole systems as observed by LIGO, Virgo and KAGRA,” Please again add the citations for the LIGO and Virgo detector papers, and also add a citation to the KAGRA detector papers (1306.6747, 1111.7185).
	+ Done: pag 14

- "hence MESA could be employed for a simultaneous inference of the PSD and of the gravitational wave signal during a Bayesian parameter estimation exercise, effectively marginalizing over the detector noise PSD”: Please add citations to other papers from LVC members discussing marginalization over PSD detector noise (0804.3853, 1109.0442, 1909.01934, 1409.7215, 1307.8195, 2004.05149, 2006.05292, 2101.01200, 2004.07515, 1506.00185 and the two Bayeswave papers mentioned above).
	+ Perhaps here there was a misundenstarding, as we did not make clear our idea. We do not aim to include in the PE the uncertainties of the PSD (as some works indicated by the referee do). Instead, our idea is to estimate the PSD at _every_ LL evaluation during a PE run. We tried to clarify this in the conclusion in the text and in a dedicated footnote.

- Please add a citation to https://arxiv.org/abs/gr-qc/0011041, which also demonstrates the applicability of autoregressive processes to gravitational-wave power spectral density estimation.
	+ Done: on sec. 5.2

- I’m not sure if this qualifies as a PnP comment, but a significant omission is the lack of discussion of the Bayeswave method of PSD estimation, which is the official method for calculating PSDs for parameter estimation results by the LVC. In addition to citing the papers described above, I think it would be appropriate to add a discussion, either in the Introduction or the Final Remarks comparing MESA to this “on-source” PSD estimation and citing the relvant papers disucssing the Bayeswave PSD method in comparison to the Welch method, some repeated from the above requests (1907.06540, 2006.05292).
	+ We found that the comparison with Bayesline goes Beyond the scope of our work. Bayesline is a parametric model (unlike the non parametric Welch and MESA) and it is taylored to GW data analysis. Our discussion here aimed to be inspired from but broader than GW data analysis. Thus, a comparison with BayesLine would narrow too much the applicability and the interest of our work.
	We mentioned this point in a newly added paragraph in the introduction (where we cited the two papers suggested).

## Personal comments

We fixed all the typos that the referee pointed out. We are thankful for that.

### Validations of the model:
- I wasn't able to get a good grasp of the practical applications of the estimation of the autoregressive process order while reading this section. Is gravitational-wave data expected to be stationary enough to reliably model as an autoregressive process? If not, what are some other examples? It would be very helpful to elaborate on this, perhaps when the connection to AR is introduced.
	+ We added at the beginning of the relevant section (3.4) a note on this, trying to explain the importance of such study.


### Comparison with Welch method:
- Maybe it’s because of the broken citation, but I don’t understand which PSD is being used in the simulation. Is it the same “analytical PSD computed for the LIGO Hanford interferometer” used in the previous section, or is it a PSD calculated using the real LIGO data from around GW150914? If so, how was that PSD calculated?
	+ The two PSDs are different. In the validation part, we used the design sensitivity of aLIGO (we quoted the relevant paper as well as the DCC link for the file). In the Welch comparison part, we employed the PSD released together with GW150914 (again we quoted the GWTC-1 paper and the link).

### Applications:

- It’s not really clear to me when you perform the forecasting for the temperature data why CAT is the obvious choice for the loss function, since in the previous section you described how CAT is unable to recover the true order of the autoregressive process. Does that mean that it would be worse at forecasting based on the AR(p)? There seems to be a contradiction where CAT does a better job for the temperature data but worse for the simulation presented in section II. It would be helpful to elaborate on this.
	+ We added a claryfing footnote on this (note 10, page 10). We hope it is clearer now.

- "Indeed, a precise prediction of the strain time series can be beneficial in the detection of anomalies in the data and, eventually, their removal.”: How does this compare to the method presented in https://arxiv.org/abs/1908.05644, where the data containing a glitch are “inpainted” with values chosen such that the inverse-PSD-filtered data is 0 at these times?
	+ To answer this matter is not easy and this should be addressed in future works. Thanks for pointing out the reference: we added it

### Final Remarks and Future Prospects

- In addition to adding the citations to other LVC works presenting methods and results for marginalizing over the PSD uncertainty (see PnP comments above), a more detailed discussion of how MESA compares to these other methods (particularly the error in the PSD measurement) would be appropriate here.
	+ We believe we discussed in details the matter in section 4, which we refer the reader to. At the beginning of the section, we addeed two lines explaining why the comparison shall be qualitative (being the results very problem dependent, we believe that every measure of error hardly generalises to other problems).















			
