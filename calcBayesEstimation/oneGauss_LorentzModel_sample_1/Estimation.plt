binWidth=0.1
numData=1000
set boxwidth binWidth
sigma=1.0
mu=0.0
plot "Bayes_oneGauss_LorentzModel_sample_1.dat" w l title "Estimation", "ND_sample_1_histo.dat" using 1:($2/(binWidth*numData)) with boxes title "Histogram", 1/(sqrt(2*pi)*sigma)*exp(-(x-mu)**2/(2*sigma**2)) title "ND(0,1)"
