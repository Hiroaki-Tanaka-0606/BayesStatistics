ND(s, m, x)=1/(sqrt(2*pi)*s)*exp(-((x-m)/s)**2/2)
binWidth=0.1
set boxwidth binWidth
num_data=10000
plot "mND_sample_2_histo.dat" using 1:($2/(num_data*binWidth)) with boxes title "Histogram"
rep 0.7*ND(1, 0, x)+0.2*ND(1, 0.8, x)+0.1*ND(1, -2, x) title "mND function"
rep 0.7*ND(1, 0, x) title "0.7*ND(1, 0)"
rep 0.2*ND(1, 0.8, x) title "0.2*ND(1, 0.8)"
rep 0.1*ND(1, -2, x) title "0.1*ND(1, -2)"