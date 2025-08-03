
# totalLikelihood_Hongkong.R

totalLikelihood <- function(obsIncidence,obsLabData,genTimeDistr,relTransmiss,rho1){
  rho_d = rep(0,length(obsIncidence));
  for(ii in 1:length(obsIncidence)){
    if(ii == 1){
      rho_d[ii] = rho1;
    }
    else{
      XAVR = relTransmiss*(rev(genTimeDistr[1:(ii-1)])%*% as.matrix(obsIncidence[1:(ii-1)]*rho_d[1:(ii-1)]));
      XAVS = rev(genTimeDistr[1:(ii-1)])%*% as.matrix(obsIncidence[1:(ii-1)]*(1-rho_d[1:(ii-1)]));
      rho_d[ii] = XAVR/(XAVR + XAVS);
    }
  }
  # monthly likelihood
  numMonth = length(obsLabData[,1]);
  numDaysEachMonth = c(0,31,30,31,31,29,31,30,31);
  sumTest = rep(0,numMonth);
  sumRes = rep(0,numMonth);
  for(jj in 1:numMonth){
    sumTest[jj] = sum(obsIncidence[(sum(numDaysEachMonth[1:jj]+1):sum(numDaysEachMonth[1:(jj+1)]))]);
    sumRes[jj] = sum(obsIncidence[(sum(numDaysEachMonth[1:jj]+1):sum(numDaysEachMonth[1:(jj+1)]))]*
                       rho_d[(sum(numDaysEachMonth[1:jj]+1):sum(numDaysEachMonth[1:(jj+1)]))]);
  }
  monthlyAvgRho = sumRes/sumTest;
  # likelihood function
  iiLikelihood = obsLabData[obsLabData[,1]>0,2]*log(monthlyAvgRho[obsLabData[,1]>0])+
    (obsLabData[obsLabData[,1]>0,1]-obsLabData[obsLabData[,1]>0,2])*log(1-monthlyAvgRho[obsLabData[,1]>0]);
  totalLogL = sum(iiLikelihood);
  return(totalLogL);
}