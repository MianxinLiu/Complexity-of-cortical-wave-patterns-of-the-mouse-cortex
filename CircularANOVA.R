library(R.matlab)
library(circular)

radiandata<-readMat("/media/user/4TB/facial/download/Kiki2/Thetas.mat")
radiandata<-radiandata$theta

data <- list(
  group1 = circular(c(radiandata[,1]), 
                 units="radians", template="geographics"),
  group2 = circular(c(radiandata[,2]),
                     units="radians", template="geographics"),
  group3 = circular(c(radiandata[,3]),
                    units="radians", template="geographics")
)

rm(radiandata)
watson.williams.test(data)

watson.wheeler.test(data)
