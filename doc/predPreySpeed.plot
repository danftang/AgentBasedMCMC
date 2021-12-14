plot "predPreySpeed.dat" using ($1*$1*$2*6):4 notitle
f(x) = m*x
fit f(x) "predPreySpeed.dat" using ($1*$1*$2*6):4 via m
replot f(x) notitle
