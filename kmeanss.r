## voice
kk <- c(seq(10, 50, 10), seq(100, 500, 100))

numc <- 5

for(k in 1:length(kk)) {
  print(kk[k])
  fpath <- 'E:/ESCUELA/CIMAT/4 Semestre/ST2/prog/mfcc/'
  # fpath <- 'C:\Users\Estudiante\Dropbox\~ESCUELA\CIMAT\4 Semestre\ST2\prog'
  ffile <- 'calderon3'
  fext <- '.csv'
  
  finput <- paste(fpath, ffile, '_ceps', fext, sep='')
  foutput <- paste(fpath, ffile, '_', kk[k], fext, sep='')
  
  aa <- as.matrix(read.csv(finput, header=F))
  
  pc <- prcomp(t(aa))
  
  cc <- pc$rotation[1:numc, ] %*% aa
  
  bb <- kmeans(t(cc), kk[k], iter.max = 200) #, algorithm=c("Hartigan-Wong"))
  ll <- bb$cluster
  write.table(t(ll), foutput, row.names=F, quote=F, col.names=F, sep=',')
}