# MCMC

MCMC_U <- function(mcIter,burnIn,obsIncidence,obsLabData,genTimeDistr,
                 parms,parmStepSize,parmLB,parmUB,minProbAccept,maxProbAccept,adjStepSize,priorA,priorB)
{
  goodchain = FALSE;
  numParms = length(parms);
  pAccept = rep((minProbAccept+maxProbAccept)/2,numParms);
  while(goodchain == FALSE)
  {
    iiRec = matrix(NA,nrow=mcIter,ncol=numParms);
    nAccept = rep(0,numParms);
    # Adjust step size
    parmStepSize[pAccept>maxProbAccept] = parmStepSize[pAccept>maxProbAccept]*(1+adjStepSize);
    parmStepSize[pAccept<minProbAccept] = parmStepSize[pAccept<minProbAccept]*(1-adjStepSize);
    # Current target function
    currParm = parms;
    currLogLikelihood =totalLikelihood(obsIncidence,obsLabData,genTimeDistr,currParm[1],currParm[2]);
    currTarget = currLogLikelihood + logFlatPrior(currParm,priorA,priorB);
    for(iiIter in 1:mcIter)
    {
      # Update one element in the parameter vector every time
      for (iiParm in 1:numParms)
      {
        tempNewParm = mcmcProposal(currParm[iiParm],parmStepSize[iiParm],parmLB[iiParm],parmUB[iiParm]);
        newParm = currParm;
        newParm[iiParm] = tempNewParm;
        # Compute new target function
        newLogLikelihood = totalLikelihood(obsIncidence,obsLabData,genTimeDistr,newParm[1],newParm[2]);
        newTarget = newLogLikelihood + logFlatPrior(newParm,priorA,priorB);
        # Metropolis-Hastings algorithm
        alphaMH = min(1,exp(newTarget-currTarget));
        uuMH = runif(1,min=0,max=1);
        if (uuMH <= alphaMH)
        {
          # Accept the new step
          currParm = newParm;
          currTarget = newTarget;
          nAccept[iiParm] = nAccept[iiParm]+1;
        }
      }
      iiRec[iiIter,] = currParm;
      # Check MCMC, every 10%
      if (iiIter%%(mcIter/10)==0)
      {
       
        pAccept = nAccept/iiIter;
        print(c('pAccept: ',as.character(round(pAccept,digits=2))));
        print(c('MCMC: ',as.character(iiIter/mcIter)));
        if (sum(pAccept>rep(maxProbAccept,numParms))+sum(pAccept<rep(minProbAccept,numParms))>0)
        {
          goodchain = FALSE;
          if (iiIter>=burnIn*mcIter)
          {
            print('Restart MCMC');
            break;
          }
        }
        else
        {
          goodchain = TRUE;
        }
      }
    }
  }
  return(iiRec);
}
