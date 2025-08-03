# source
set.seed(123)
install.packages("gridExtra")
library(gridExtra)
library('ggplot2')
library('reshape2')
source("totalLikelihood_HongKong.R");
source("negTotalLikelihood_HongKong.R");
source("logFlatPrior.R");
source("mcmcProposal.R");
source("MCMC_U.R");
# input parameters
genTime=2;  # the mean of generation time
sdGenTime=1.5; # the sd of generation time 
relTransmiss=1.1; # initial guess of AVR fitness
rho1=0.34; # initial guess of p(t)
mcIter=10000; # number of MCMC iters
# define days and weeks
firstDay = '24/09/2007';
firstDayEveryMonth = c('01/10/2007','01/11/2007','01/12/2007','01/01/2008','01/02/2008','01/03/2008','01/04/2008','01/05/2008');
weekIndex = as.numeric((as.Date(firstDayEveryMonth,'%d/%m/%Y')-as.Date(firstDay,'%d/%m/%Y')+14)/7);
weekLabel = c('Oct 2007','Nov 2007','Dec 2007','Jan 2008','Feb 2008','Mar 2008','Apr 2008','May 2008');

# read data files

# observed daily overall incidence data i(t) 
weeklyFluActivity = read.csv('weekly_country_flu_H1_activity.csv',header = TRUE,stringsAsFactors = FALSE); 

# observed laboratory resistance data Z 
monthlyResData = read.csv('monthly_country_resistance_data.csv',header = TRUE,stringsAsFactors = FALSE);

# part 1. deconvoluting the time series of cases by dates, weekly to daily
fluIsolate = as.numeric(weeklyFluActivity[1,3:length(weeklyFluActivity[1,])]);
# linear adjustment
numWeek = length(fluIsolate);
obsIncidence = rep(0,7*numWeek);
# 4th day of the week
obsIncidence[7*((1:numWeek)-1)+4] = fluIsolate[1:numWeek]/3.5;
# 7th day and 1st of the week
obsIncidence[7*((2:numWeek)-2)+7] = 0.5*(as.numeric(obsIncidence[7*((2:numWeek)-2)+4])+as.numeric(obsIncidence[7*((2:numWeek)-1)+4]));
obsIncidence[7*((2:numWeek)-1)+1] = 0.5*(as.numeric(obsIncidence[7*((2:numWeek)-2)+4])+as.numeric(obsIncidence[7*((2:numWeek)-1)+4]));
# 1st-3rd day of the week
# week 1
obsIncidence[1:3] = as.numeric(obsIncidence[4])/3*(1:3);
# other weeks
obsIncidence[(7*((2:numWeek)-1)+2)] = as.numeric(obsIncidence[7*((2:numWeek)-1)+1])+((as.numeric(obsIncidence[7*((2:numWeek)-1)+4])-as.numeric(obsIncidence[7*((2:numWeek)-1)+1]))/3);
obsIncidence[(7*((2:numWeek)-1)+3)] = as.numeric(obsIncidence[7*((2:numWeek)-1)+1])+((as.numeric(obsIncidence[7*((2:numWeek)-1)+4])-as.numeric(obsIncidence[7*((2:numWeek)-1)+1]))/3)*2;
obsIncidence[(7*((2:numWeek)-1)+5)] = as.numeric(obsIncidence[7*((2:numWeek)-1)+4])+((as.numeric(obsIncidence[7*((2:numWeek)-1)+7])-as.numeric(obsIncidence[7*((2:numWeek)-1)+4]))/3);
obsIncidence[(7*((2:numWeek)-1)+6)] = as.numeric(obsIncidence[7*((2:numWeek)-1)+4])+((as.numeric(obsIncidence[7*((2:numWeek)-1)+7])-as.numeric(obsIncidence[7*((2:numWeek)-1)+4]))/3)*2;
obsIncidence = as.numeric(obsIncidence);
obsIncidence = pmax(obsIncidence,1e-7);


# part 2. antiviral resistance data
temp = monthlyResData[1,3:length(monthlyResData[1,])];
obsLabData = t(rbind(as.numeric(temp[seq(1,length(temp),2)]),as.numeric(temp[seq(2,length(temp),2)])));

# part 3. combine all data for next step
allData = list(fluIsolate,obsLabData);
# flu isolates
fluIsolate = as.data.frame(allData[1]);
colnames(fluIsolate) = 'Isolate';
fluIsolate = cbind(Week=1:length(fluIsolate$Isolate),fluIsolate);
# AVR data
AVRTest = as.data.frame(allData[2]);
colnames(AVRTest) = c('Total','Resistant');
AVRTest[,'Sensitive'] = AVRTest$Total-AVRTest$Resistant;


# barplot
AVRTest = cbind(Month=1:length(AVRTest$Total),AVRTest);
barplotData = melt(AVRTest[,c('Month','Sensitive','Resistant')],id.vars = 'Month');
colnames(barplotData) = c('Month','Antiviral','Number');
plot1  = ggplot(barplotData,aes(x = Month,y = Number,fill = Antiviral)) + 
  geom_bar(stat = "identity") + 
  ylab('Number of AVS/AVR isolates')+
  xlab('Month')+
  expand_limits(x=c(1,8))+
  scale_x_continuous(breaks = 1:8,labels = weekLabel)+
  theme(legend.position = c(0.8, 0.8));

plot1

# lineplot
plot2 = ggplot(fluIsolate[,c('Week','Isolate')])+ 
  geom_line(aes(x = Week, y = Isolate), size = 0.5, alpha = 0.75) + 
  scale_x_continuous(breaks = weekIndex,labels = weekLabel)+
  expand_limits(x=c(1,weekIndex[8]))+
  ylab('Number of A(H1N1) isolates')+
  xlab('Week');
plot2

grid.arrange(plot1, plot2, ncol = 1)


# part 4. generation time distribution
numDay = length(obsIncidence);
genTimeDistr = rep(0,numDay);
genTimeDistr[1] = plnorm(0.5,meanlog = log(genTime),sdlog = log(sdGenTime)); # using log-normal distribution
genTimeDistr[2:length(genTimeDistr)] = 
  plnorm((2:length(genTimeDistr))+0.5,meanlog = log(genTime),sdlog = log(sdGenTime))-
  plnorm((2:length(genTimeDistr))-0.5,meanlog = log(genTime),sdlog = log(sdGenTime));

# part4. non-linear optimization for the initial values
parmLB = c(0.7,1e-7);
parmUB = c(1.3,1-1e-7);
resOptim = optim(par = c(relTransmiss,rho1),
                 fn = negTotalLikelihood,
                 lower = parmLB,
                 upper = parmUB,
                 countryText = countryText,
                 obsIncidence = obsIncidence,
                 obsLabData = obsLabData,
                 genTimeDistr = genTimeDistr,
                 method='L-BFGS-B');

# part5. MCMC for final results
relTransmiss = resOptim$par[1];
rho1 = resOptim$par[2];
relTransmissStep = 0.5/numWeek;
rho1Step =  max(rho1*20/numWeek,0.01);
parms = c(relTransmiss,rho1);
parmStepSize = c(relTransmissStep,rho1Step);
burnIn = 0.2;
minProbAccept = 0.3;
maxProbAccept = 0.7;
adjStepSize = 0.1;
priorA = c(0.9,1e-7);
priorB = c(1.1,1);
mcmcRes = MCMC_U(mcIter,
               burnIn,
               obsIncidence,obsLabData,genTimeDistr,
               parms,
               parmStepSize,
               parmLB,parmUB,minProbAccept,maxProbAccept,adjStepSize,priorA,priorB);

allData = list(fluIsolate,obsLabData,mcmcRes[(burnIn*length(mcmcRes[,1])):length(mcmcRes[,1]),]);
# flu isolates
fluIsolate = as.data.frame(allData[1]);
colnames(fluIsolate) = 'Isolate';
fluIsolate = cbind(Week=1:length(fluIsolate$Isolate),fluIsolate);
# AVR data
AVRTest = as.data.frame(allData[2]);
colnames(AVRTest) = c('Total','Resistant');
AVRTest[,'Sensitive'] = AVRTest$Total-AVRTest$Resistant;
# MCMC data
resMCMC = as.data.frame(allData[3]);
colnames(resMCMC) = c('fitness_AVR','t0_p_AVR');
# MCMC plot
library(ggplot2)
plot3 = ggplot(resMCMC, aes(resMCMC$fitness_AVR)) + 
  geom_histogram()+
  xlab("Fitness estimate")+
  ylab("Count");
# Combine
plot.new()
plot3;

# show all results
grid.arrange(plot1, plot2,plot3, ncol = 1)



