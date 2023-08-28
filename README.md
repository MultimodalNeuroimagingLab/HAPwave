# HAPwave
Repository for code related to the publication of: 
- Ojeda Valencia G, Gregg N, Huang H, Lundstrom B, Brinkmann B, Pal Attia T, Van Gompel J, Bernstein M, In MH, Huston J, Worrell G, Miller KJ, and Hermes D. 2023. Signatures of electrical stimulation driven network interactions in the human limbic system. _Journal of Neuroscience_ (in press).

# Task Description
Patients were resting in the hospital bed, while single pulse stimulation was performed. The stimulation had a duration of 200 microseconds, was biphasic and had an amplitude of 6mA. For subject 7 stimulation amplitude was sometimes reduced to 4mA to minimize interictal responses.

# Dataset
This data is organized according to the Brain Imaging Data Structure specification ([BIDS](https://bids-specification.readthedocs.io/en/stable/) version 1.12.0) and available on [OpenNeuro](https://openneuro.org/datasets/ds004696)

# Dependencies
This code is dependent on the following packages:
- [MatMef](https://github.com/MaxvandenBoom/matmef)
- [MNL_ieegbasics](https://github.com/MultimodalNeuroimagingLab/mnl_ieegBasics)
- [Vistasoft](https://github.com/vistalab/vistasoft)

# How to run
In order to analyse the data, and reproduce figures from the manuscript, use the following steps:
- Analyze the data and calculate statistics: pcc_Allstats_CRP.m
- Make figure 2: pcc_fig2_allLimbicN1.m 
- Make figure 3: pcc_fig3b_singleTrials.m
- Make figure 4: pcc_fig4_sigCCEPs.m
- Make figure 5: pcc_fig5_signatureShift.m
- Make figure 6: pcc_fig6_temporalDelay.m
- Make figure 7: pcc_fig7c_xCorrelation.m



# Acknowledgements
This project was funded by the National Institute Of Mental Health of the National Institutes of Health Brain Initiative under Award Number R01 MH122258, “CRCNS: Processing speed in the human connectome across the lifespan". 

# Contact
Please contact [Dora Hermes](https://github.com/dorahermes) (hermes.dora@mayo.edu) or [Gabriela Ojeda Valencia](https://github.com/GabOjVa) (OjedaValencia.Alma@mayo.edu) for questions.
