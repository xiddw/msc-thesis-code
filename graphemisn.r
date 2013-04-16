library(ggplot2)

setwd("E:/ESCUELA/CIMAT/4 Semestre/ST2/prog/pruebas/prb3")
ems <- read.csv("6to7emsn.mat", header=F)
ems <- data.matrix(ems)

n <- dim(ems)[1];
w <- dim(ems)[2];

ww <- seq(w);
  
dfems <- data.frame();

for(i in seq(n)) {
  b <- data.frame(word=ww,
                  speaker=rep(paste('Speaker', i), w),
                  prc=ems[i, ])
  
  rownames(b) <- w*(i-1) + seq(w)
  
  dfems <- rbind(dfems, b)
}

ggplot(dfems) + 
  geom_area(aes(word, prc, group=speaker, colour=speaker, fill=speaker), 
            position="identity", alpha=0.4)
