*********************************************************************
** REPLICATION FILES FOR CALONICO, CATTANEO AND FARRELL (2018, JASA)
** https://doi.org/10.1080/01621459.2017.1285776
** Software: https://nppackages.github.io/nprobust/
*********************************************************************

This repository includes R scripts to replicate the simulation study reported in the main paper and in the online supplemental appendix. 


Kernel Density Estimation:

1) Replication file is "simuls_kd.R" (the other scripts are called by this file). Choices are: model (model = 1,2,3,4), evaluation point (eval = 1,2,3,4,5) and three scripts for different tables:

	(a) simuls_kd_grid.R: to run simulations over a bandwidth grid.

	(b) simuls_kd_bws.R: to run simulations using population and data-driven bandwidths.

	(c) simuls_kd_3D.R: to run simulations over grids of bandwidth and rho (ratio h/b).

2) Each script generates a .csv file in the output folder.



Local Polynomial Estimation:

1) Replication file is "simuls_lp.R" (the other scripts are called by this file). Choices are: and model (model 1,2,3,4,5,6), evaluation point (eval = 1,2,3,4,5) and three scripts for different tables:

	(a) simuls_lp_grid.R: to run simulations over a bandwidth grid.

	(b) simuls_lp_bws.R: to run simulations using population and data-driven bandwidths.

	(c) simuls_lp_comp.R: to run simulations using MSE data-driven bandwidths with lambda = (0.5, 0.7, 1).

2) Each script generates a .csv file in the output folder.

3) To conduct local polynomial estimation on a particular dataset, we provide the nprobust package:

	(a) lprobust: local polynomial estimation. The main especification implemented in the simulations is:

		lprobust(y, x, c=0, p=1, q=2, deriv=0, kernel="epa", vce="hc3")

	(b) lpbwselect: bandwidth selection for local polynomial estimation using MSE or RBC criterion. The main especification implemented in the simulations is:

		lpbwselect(y, x, c=0, p=1, q=2, deriv=0, kernel="epa", vce="hc3")

	for MSE bandwidth selection. Alternatively, bandwidth selection for local polynomial estimation using CER criterion with a DPI implementation is implemented with:

		lpbwselect(y, x, c=0, p=1, q=2, deriv=0, kernel="epa", vce="hc3", bwselect = "cer")

Additional information is included in the package documentation.



