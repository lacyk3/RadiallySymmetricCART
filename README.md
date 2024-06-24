# RadiallySymmetricCART
This repository houses code used to generate results discussed in the manuscript *Spatiotemporal dynamics of tumor - CAR T cell interaction following local administration in solid cancers*. Our reaction-diffusion model for tumor and CAR T cell interaction was solved via finite differences on a radially symmetric domain. The repository is organized as follows:

## TumorGrowth: 
contains code for simulating tumor growth in the absence of treatment
- code used to create Figure 2
- code used to generate tumors of different sizes which were then used to simulate treatment
- code used to generate Figure S1

## IntratumoralTreatment:
contains code for simulating intratumoral injection of CAR T cells

- code used to create Figure 3
- code used to create Figure 6 row 1
- code used to generate Figure S2A
- code used to generate Figure S3A
- code used to generate Figure S4
- code used to generate Figure S5-6

## IntracavitaryTreatment: 
contains code for simulating intracavitary administration of CAR T cells to treat solid tumors
- code used to generate Figure 4
- code used to generate Figure 6 row 2
- code used to generate Figure S2B
- code used to generate Figure S3B

## Fit to Data:
contains code used to fit model to murine study data
  - Zhao et al. -- Figure 5A
  - Skovgard et al. -- Figure 5B
  
