#setwd(path where everything is saved)

source('fitr.R') #load the model fitting functions

#mTor model fitting

#load in the data
mtor=read.csv('mtor.csv')

#prettify the data
mdoses <- unique(mtor$DoseNum)
cols <- paste0("gray",round(seq(from=30,to=80,length=length(mdoses))))
mtor.agg <- aggregate(logCellNum ~ Day+DoseNum+Cells,mtor,mean)

#intialize vectors to save fitted b values
mtorN = vector(length=length(mdoses))
mtorA = vector(length=length(mdoses))

#fix a carrying capacity
K_Fixed=206822


par(mfrow=c(1,2)) #to plot 2 figs at once

plot(logCellNum ~ Day,mtor,pch="",subset=Cells=="N")
title(main="mTOR inhibitor: Non-adapted cells and Model Fit")
for (idose in 1:length(mdoses)) {
  points(logCellNum ~ Day,mtor,subset=Cells=="N" & DoseNum==mdoses[idose],
         pch=19,col=cols[idose])
  #lines(logCellNum ~ Day,mtor.agg,subset=Cells=="N" & DoseNum==mdoses[idose],
  #       lwd=2,col=cols[idose])
  N0_b = get_N0_and_b('mtor',mdoses[idose],'N',K_Fixed)
  mtorN[idose] = N0_b[2]
  logMod = function(t) log(logisticModel(t,0,N0_b[1],N0_b[2],K_Fixed))
  curve(logMod, add=T,lwd=2,col=cols[idose])
}

plot(logCellNum ~ Day,mtor,pch="",subset=Cells=="A")
title(main="mTOR inhibitor: Adapted cells and Model Fit")
for (idose in 1:length(mdoses)) {
  points(logCellNum ~ Day,mtor,subset=Cells=="A" & DoseNum==mdoses[idose],
         pch=19,col=cols[idose])
  #lines(logCellNum ~ Day,mtor.agg,subset=Cells=="A" & DoseNum==mdoses[idose],
  #       lwd=2,col=cols[idose])
  N0_b = get_N0_and_b('mtor',mdoses[idose],'A',K_Fixed)
  mtorA[idose] = N0_b[2]
  logMod = function(t) log(logisticModel(t,0,N0_b[1],N0_b[2],K_Fixed))
  curve(logMod, add=T,lwd=2,col=cols[idose])
}
par(mfrow=c(1,1))

mtorRs = data.frame(cbind(mdoses,mtorN,mtorA))

# ribo model fitting

ribo=read.csv('ribo.csv')
rdoses <- unique(ribo$DoseNum)
cols <- paste0("gray",round(seq(from=30,to=80,length=length(rdoses))))
ribo.agg <- aggregate(logCellNum ~ Day+DoseNum+Cells,ribo,mean)

riboN = vector(length=length(rdoses))
riboA = vector(length=length(rdoses))

windows()
par(mfrow=c(1,2))

plot(logCellNum ~ Day,ribo,pch="",subset=Cells=="N")
title(main="Ribociclib: Non-adapted cells & Model Fit")
for (idose in 1:length(rdoses)) {
  points(logCellNum ~ Day,ribo,subset=Cells=="N" & DoseNum==rdoses[idose],
         pch=19,col=cols[idose])
  #lines(logCellNum ~ Day,ribo.agg,subset=Cells=="N" & DoseNum==rdoses[idose],
  #       lwd=2,col=cols[idose])
  N0_b = get_N0_and_b('ribo',rdoses[idose],'N',K_Fixed)
  riboN[idose] = N0_b[2]
  logMod = function(t) log(logisticModel(t,0,N0_b[1],N0_b[2],K_Fixed))
  curve(logMod, add=T,lwd=2,col=cols[idose])
}

plot(logCellNum ~ Day,ribo,pch="",subset=Cells=="A")
title(main="Ribociclib: Adapted cells & Model Fit")
for (idose in 1:length(rdoses)) {
  points(logCellNum ~ Day,ribo,subset=Cells=="A" & DoseNum==rdoses[idose],
         pch=19,col=cols[idose])
  #lines(logCellNum ~ Day,ribo.agg,subset=Cells=="A" & DoseNum==rdoses[idose],
  #       lwd=2,col=cols[idose])
  N0_b = get_N0_and_b('ribo',rdoses[idose],'A',K_Fixed)
  riboA[idose] = N0_b[2]
  logMod = function(t) log(logisticModel(t,0,N0_b[1],N0_b[2],K_Fixed))
  curve(logMod, add=T,lwd=2,col=cols[idose])
}
par(mfrow=c(1,1))

riboRs = data.frame(cbind(rdoses,riboN,riboA))
