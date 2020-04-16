bin(x,width)=width*floor(x/width)+width/2.0
binwidth = 0.1
sigma=2.0
mu=1.0
stats "normalDistribution.dat"
plot "" using (bin($1,binwidth)):(1.0/(binwidth*STATS_records)) smooth freq with boxes title "Program Output", 1.0/(sqrt(2*pi)*sigma)*exp(-1.0/2*((x-mu)/sigma)**2) w l lw 3 title "Normal Distribution"