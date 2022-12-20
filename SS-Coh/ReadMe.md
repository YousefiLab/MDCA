# SS-Coh 

In State Space Coherence (SS-Coh) project, we provide mor robust algorithm for estimating Global Coherence (GCoh). 
To understand some aspects of SS-Coh modelling, you can read [EMBC-2019](EEG%20data%20from%20human%20patients%20under%20general%20anesthesi) and [bioRxiv-2020](https://www.biorxiv.org/content/10.1101/2020.07.13.199034v1.abstract). In the RO1 proposal, we provide more complete modelling of SS-Coh. 


## Data
To assess SS-Coh model, we use EEG data from human patients under general anesthesia. The data set was collected in Emery Brown’s laboratory. 

The complete description of the experimental protocol can be found in [Purdon et al .](https://www.pnas.org/doi/10.1073/pnas.1221180110) . Briefly, ten consenting human volunteers of ages 18-36 years were impaneled for the study approved by the MGH Human Research Committee. For each subject, the induction and emergence from propofol anesthesia were studied by administering a computer-controlled (StanPump) infusion of propofol using the target control protocol based on the Schnider pharmacokinetic-pharmacodynamic model, while the subject executed a behavioral task to identify the points of loss and recovery of consciousness. Neural activity was recorded from 64 channels of EEG at a 250 Hz sampling rate.

Using this dataset, we will study the causal relationship between propofol blood concentration, level of consciousness, and spatio-temporal patterns of functional connectivity. Spectrograms, sliding window GCoh, and LDCM will be applied to identify network modes and their changes associated with loss of consciousness.

## Results

![ **Figure 1. SS-GCoh analysis in Alpha band for anesthesia EEG data and its correspondence with empirical results. A)** Inferred unconsciousness state estimation over 2 hours of anesthesia. The latent state represents the unconsciousness level; the state transition is modeled by a random-walk model. **B)** Scalp heat-map of the dominant eigenmodes for 3 different time points during the experiment. The result using SS-GCoh and empirical measures are similar to each other. **C)** Empirical GCoh and inferred GCoh using SS-GCoh. Not only the inferred coherence matches the empirical one; it attains it at a finer temporal resolution. With SS-GCoh, we also derive higher-order statistics of the coherence such as confidence interval.](imgs/fig1.jpg)