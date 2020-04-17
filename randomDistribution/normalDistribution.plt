bin(x,width)=width*floor(x/width)+width/2.0
binwidth = 0.1
sigma=1.0
mu=0.0
stats "ND_sample_1.dat"
set boxwidth binwidth
plot "" using (bin($1,binwidth)):(1.0/(binwidth*STATS_records)) smooth freq with boxes title "Program Output", 1.0/(sqrt(2*pi)*sigma)*exp(-1.0/2*((x-mu)/sigma)**2) w l lw 3 title "Normal Distribution"