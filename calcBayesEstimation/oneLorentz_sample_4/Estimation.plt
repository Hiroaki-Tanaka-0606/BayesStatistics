binWidth=0.1
numData=100000
set boxwidth binWidth
sigma=1.0
mu=0.0
set xr [mu-5*sigma:mu+5*sigma]
plot "Esti_oneLorentz_sample_4.dat" w l title "Estimation", "LD_sample_4_histo.dat" using 1:($2/(binWidth*numData)) with boxes title "Histogram", sigma/(pi*((x-mu)**2+sigma**2)) title "LD(0,1)"
