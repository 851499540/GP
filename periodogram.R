#periodogram

hours <- 1:120
fdata <- as.ts(cos(pi*(hours/12)))
plot(fdata)
fft(fdata) #fast fourier transform
f1 <- spec.pgram(fdata, log='no', plot = FALSE)
f1$freq <- f1$freq*24
plot(f1$freq, f1$spec, type = 'l')


fdata <- as.ts(sin(pi*(hours/6)) + cos(pi*(hours/12)))
plot(fdata)
f2 <- spec.pgram(fdata, log='no', plot = FALSE)
f2$freq <- f2$freq*24
plot(f2$freq, f2$spec, type = 'l')

fdata <- as.ts(0.5*sin(pi*(hours/6)) + cos(pi*(hours/12)))
plot(fdata)
f3 <- spec.pgram(fdata, log='no', plot = FALSE)
f3$freq <- f3$freq*24
plot(f3$freq, f3$spec, type = 'l')

fdata <- as.ts(0.5*sin(pi*(hours/6)) + cos(pi*(hours/12)) + 0.1*cos(pi*(hours/2)))
plot(fdata)
plot(0.1*cos(pi*hours/2), type = 'l')
f4 <- spec.pgram(fdata, log='no', plot = FALSE)
f4$freq <- f4$freq*24
plot(f4$freq, f4$spec, type = 'l')

fdata <- as.ts(0.5*sin(pi*(hours/6)) + cos(pi*(hours/12)) + 0.1*cos(pi*(hours/2)) + rnorm(120))
plot(fdata)
f5 <- spec.pgram(fdata, log='no', plot = FALSE)
f5$freq <- f5$freq*24
plot(f5$freq, f5$spec, type = 'l')

# not consist when estimate, always the same even has more data
#smooth, use tapering


f5 <- spec.pgram(fdata, log='no',spans = c(5,5), plot = FALSE)
f5$freq <- f5$freq*24
plot(f5$freq, f5$spec, type = 'l')

#over smooth, lose signal, then use taper to less smooth

f5 <- spec.pgram(fdata, log='no',spans = c(5,5), plot = FALSE, taper = 0.5)
f5$freq <- f5$freq*24
plot(f5$freq, f5$spec, type = 'l')

#leakage cause the problem : the spiky peak smooth out its area. freq1 leak to freq 2
#solution1: filter
#solution2: welch's method
#     try to imporve periodogram estimate
#     idea: chop the time series to several blocks, calculate periodogram separately, average them
#          remove some noise

#SOI
require(astsa)
data(soi)
plot(soi)

soi.per <- spec.pgram(soi, taper=0, log='no')
soi.persmooth <- spec.pgram(soi, spans=5, taper=0, log='no')
soi.persmoothln <- spec.pgram(soi, spans=9, taper=0, log='yes')


s0 = spec.pgram(soi, spans=c(5,5), log='no', plot=FALSE)
s50 =spec.pgram(soi,taper=0.5, spans=c(5,5), log='no', plot=FALSE)
plot(s0$freq, s0$spec, log='y', type='l', lty=2, ylab='spectrum',
     xlab='frequency')
lines(s50$freq,s50$spec)

#taper make it less smooth




#welch's meathod
soi1 = soi[1:227]
soi2 = soi[228:453]
s1 = spec.pgram(soi1, taper=0, log='no')
s2 = spec.pgram(soi2, taper=0, log='no')

soi.ave = (s1$spec+s2$spec)/2
plot(s1$freq, soi.ave, type='l')





