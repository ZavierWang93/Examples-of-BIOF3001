
logFlatPrior <- function(parms,vecAlpha,vecBeta)
{
  # Flat priors
  logP = -log(vecBeta-vecAlpha);
  return(sum(logP));
}