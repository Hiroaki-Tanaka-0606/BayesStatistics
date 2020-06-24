ND(s, m, x)=1/(sqrt(2*pi)*s)*exp(-((x-m)/s)**2/2)
binWidth=0.1
set boxwidth binWidth
num_data=10000
plot "mND_sample_3_histo.dat" using 1:($2/(num_data*binWidth)) with boxes title "Histogram"
rep 0.5*ND(1, 0, x)+0.3*ND(1, 3, x)+0.2*ND(1, -4, x) title "mND function"
rep 0.5*ND(1, 0, x) title "0.5*ND(1, 0)"
rep 0.3*ND(1, 3, x) title "0.3*ND(1, 3)"
rep 0.2*ND(1, -4, x) title "0.2*ND(1, -4)"