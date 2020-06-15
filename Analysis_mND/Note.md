# Programs
X = number of normal distributions
- RandomDistribution_#{X}ND.o
- RXMC_opt_#{X}ND.o
- RX_multiple_#{X}ND.o
- RXMC_#{X}ND.o
- Average_analysis.o
- PSRF.o
- Free_energy_1_#{X}ND.o
- Free_energy_2_#{X}ND.o
- Posterior_#{X}ND.o (from RX_multiple_posterior_mND.o)
- WAIC_#{X}ND.o (from RX_multiple_WAIC_mND.o)

We have to change the script name "RX_multiple_mND.o" in RX_multipls.sh to "RX_multiple_#{X}ND.o" before executing.