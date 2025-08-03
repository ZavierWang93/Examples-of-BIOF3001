
# for "optim" hongkong
negTotalLikelihood <- function(x,countryText,obsIncidence,obsLabData,genTimeDistr){
  relTransmiss = x[1];
  rho1 = x[2];
    return(-totalLikelihood(obsIncidence,obsLabData,genTimeDistr,relTransmiss,rho1));
}