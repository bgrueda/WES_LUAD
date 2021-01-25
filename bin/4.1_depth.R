setwd("~/Desktop")
tab <- read.table ("depth.txt", 
                   header= F)
attach (tab)
names (tab)

pdf("depth1.pdf", width=6, height=3)

x <- c(V2)
y <- c(V3)

plot (x,y, 
      main="Depth",
      col= 'navy', 
      type= "l")

dev.off()

