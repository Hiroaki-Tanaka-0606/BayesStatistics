bin(x,width)=width*floor(x/width)+width/2.0
binwidth = 0.1
sigma=1.0
mu=0.0
stats "LD_sample_1.dat"
set xrange [mu-5*sigma:mu+5*sigma]
set boxwidth binwidth
plot "LD_sample_1.dat" using (bin($1,binwidth)):(1.0/(binwidth*STATS_records)) smooth freq with boxes title "Program Output", sigma/(pi*(sigma**2+(x-mu)**2)) w l lw 3 title "Lorentz Distribution"