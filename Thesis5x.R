library(ggplot2) # Used to graph efficient frontier, optional
library(quadprog) #Needed for solve.QP
library(plyr)
eff.frontier <- function (returns, short="no", max.allocation=NULL,
                          risk.premium.up=.5, risk.increment=.001, dimx = NULL, 
                          dimlength = NULL)
{  
  covariance <- cov(returns)
  #print(covariance)
  n <- ncol(covariance)
  Amat <- matrix (1, nrow=n)
  bvec <- 1
  meq <- 1
  if(short=="no"){
    Amat <- cbind(1, diag(n))
    bvec <- c(bvec, rep(0, n))
  }
  if(!is.null(max.allocation)){
    if(max.allocation > 1 | max.allocation <0){
      stop("max.allocation must be greater than 0 and less than 1")
    }
    if(max.allocation * n < 1){
      stop("Need to set max.allocation higher; not enough assets to add to 1")
    }
    Amat <- cbind(Amat, -diag(n))
    bvec <- c(bvec, rep(-max.allocation, n))
  }
  loops <- risk.premium.up / risk.increment + 1
  loop <- 1
  eff <- matrix(nrow=loops, ncol=n+3)
  colnames(eff) <- c(dimx[2:dimlength], "Std.Dev", "Exp.Return", "sharpe")
  # Loop through the quadratic program solver
  for (i in seq(from=0, to=risk.premium.up, by=risk.increment)){
    dvec <- colMeans(returns) * i # This moves the solution along the EF
    sol <- solve.QP(covariance, dvec=dvec, Amat=Amat, bvec=bvec, meq=meq)
    eff[loop,"Std.Dev"] <- sqrt(sum(sol$solution*colSums((covariance*sol$solution))))
    eff[loop,"Exp.Return"] <- as.numeric(sol$solution %*% colMeans(returns))
    eff[loop,"sharpe"] <- eff[loop,"Exp.Return"] / eff[loop,"Std.Dev"]
    eff[loop,1:n] <- sol$solution
    loop <- loop+1
  }
  return(as.data.frame(eff))
}
regression.model = function(dimx = NULL,dimlength = NULL, newX = NULL)
{
  regmodelpred = dimx[2]
  for(i in 3:dimlength) 
  {
    regmodelpred = paste(regmodelpred,"+",dimx[i])
    next
  }
  regmodelpred = paste(regmodelpred,"-1")
  RegModel = lm(as.formula(paste(dimx[1],"~",regmodelpred)), data = newX)
  Regression.Weighting = 
    as.matrix(coefficients(RegModel))/(sum(coefficients(RegModel)))
  RegWeights = as.data.frame(t(Regression.Weighting), 
                             row.names = "Regression.Weighting")
  return(RegWeights)
}
unit.model = function(dimlength = NULL, dimx = NULL)
{
  Unit.Weighting = as.matrix(rep((1/(dimlength-1)),(dimlength-1)))
  UnitWeights = t(data.frame(Unit.Weighting, row.names = dimx[2:dimlength]))
  return(UnitWeights)
}
applicants = function(PerfMin = NULL, PerfMax = NULL, dimx = NULL, dimlength = NULL,                       
                      R = NULL, numapps = NULL, SEED = NULL)
{
  U = t(chol(R)) #Creates the Cholesky decomposition of the correlation matrix
  #This next section creates a a random data set
  nvars = dim(U)[1]  #number of variables
  set.seed(SEED+1)  #Choose a seed to generate random numbers
  random.normal = matrix(rnorm(nvars*numapps*100,0,1), nrow=nvars, ncol=numapps*100);
  j = sample(1:100*numapps, numapps)
  random.normal = random.normal[,j]
  X = U %*% random.normal 
  newX = as.data.frame(t(X)) 
  Xranked = newX[order(newX$Performance),] 
  #rescale the predictors
  Xrescaled = (apply(Xranked, MARGIN = 2, 
                     FUN = function(X) (X - min(X))/diff(range(X))))*(PerfMax-PerfMin)+
                    (matrix(rep(9,numapps*dimlength),ncol = dimlength))
  apps.mat = c(Xrescaled, newX)
  return(apps.mat)
  ### Plot the correlations
  # raw = as.data.frame(newX)
  # names(raw) = dimx
  # cor(raw)
  # plot(head(raw, 100)) #plot the corelations 
}
performance = function(PerfMin = NULL, PerfMax = NULL, dimx = NULL, 
                       dimlength = NULL, R = NULL, numobs = NULL, SEED = NULL)
{
  U = t(chol(R)) #Creates the Cholesky decomposition of the correlation matrix
  #This next section creates a a random data set
  nvars = dim(U)[1]  #number of variables
  set.seed(SEED)  #Choose a seed to generate random numbers
  random.normal = matrix(rnorm(nvars*numobs,0,1), nrow=nvars, ncol=numobs); 
  X = U %*% random.normal 
  newX = as.data.frame(t(X)) #transpose the result
  Xranked = newX[order(newX$Performance),] #rank orders the results by performance
  #rescale the predictors
  Xrescaled = (apply(Xranked, MARGIN = 2, 
                     FUN = function(X) (X - min(X))/diff(range(X))))*(PerfMax-PerfMin)+
                    (matrix(rep(9,numobs*dimlength),ncol = dimlength))
  performance.mat = c(Xrescaled, newX)
  return(performance.mat)
  ### Plot the correlations
  # raw = as.data.frame(newX)
  # names(raw) = dimx
  # cor(raw)
  # plot(head(raw, 100)) #plot the corelations 
}
get.returns = function(Xrescaled = NULL, dimlength = NULL)
{
  Xpred = Xrescaled[,2:dimlength] #Removes the performance variable
  returns = (tail(Xpred,-1) - head(Xpred,-1))
  return(returns)
  #hist(returns)
}
##Diagnostic/Graphic
# graph.eff = function(eff = NULL)
# {
#   # Find the optimal portfolio
#   eff.optimal.point <- eff[eff$sharpe==max(eff$sharpe),]
#   # graph efficient frontier
#   # Start with color scheme
#   ealred <- "#7D110C"
#   ealtan <- "#CDC4B6"
#   eallighttan <- "#F7F6F0"
#   ealdark <- "#423C30" 
#   eff.plot = ggplot(eff, aes(x=Std.Dev, y=Exp.Return)) + geom_point(alpha=.1, 
#     color=ealdark) +
#     geom_point(data=eff.optimal.point, aes(x=Std.Dev, y=Exp.Return, label=sharpe),
#                color=ealred, size=5) +
#     annotate(geom="text", x=eff.optimal.point$Std.Dev,
#              y=eff.optimal.point$Exp.Return,
#              label=paste("Risk: ",
#                          round(eff.optimal.point$Std.Dev*100, digits=3),"\nReturn: ",
#                          round(eff.optimal.point$Exp.Return*100, digits=4),"%\nSharpe: ",
#                          round(eff.optimal.point$sharpe*100, digits=2), "%", sep=""), 
#              hjust=1.5, vjust=2) +
#     ggtitle("Efficient Frontier\nand Optimal Portfolio") +
#     labs(x="Risk (standard deviation of portfolio)", y="Return") +
#     theme(panel.background=element_rect(fill=eallighttan),
#           text=element_text(color=ealdark),
#           plot.title=element_text(size=24, color=ealred))
#   ggsave("Efficient Frontier.png")
#   return(eff.plot)
# } 
MPT.model = function(eff = NULL, dimx = NULL, dimlength = NULL)
{
  EFFSolutionsHead = head(eff)
  MPTWeights = EFFSolutionsHead[1,1:(dimlength -1)]
  rownames(MPTWeights) = "MPT.Weighting"
  #MPTRetAndSD = (EFFSolutionsHead[1,dimlength:(dimlength+1)])
  # print(MPTRetAndSD)
  return(MPTWeights)
}
SignalDetection = function(dimx = NULL, Perf = NULL, PerfMin = NULL, PerfMax = NULL, 
                           Signal = NULL, Cutoff = NULL, numapps = NULL)
{
  Serf = Perf
  for(i in 1:numapps)
  {
    for(x in 2:4)
      {
        if((Perf[i,x] > Cutoff) & (Perf[i,1] > Cutoff)) Serf[i,x] = 1
        if((Perf[i,x] < Cutoff) & (Perf[i,1] < Cutoff)) Serf[i,x] = 2
        if((Perf[i,x] > Cutoff) & (Perf[i,1] < Cutoff)) Serf[i,x] = 3
        if((Perf[i,x] < Cutoff) & (Perf[i,1] > Cutoff)) Serf[i,x] = 4
      }
  }
  Signal = (Serf[,2:4])
  return(Signal)
}
#########Simulation 1################################################################
Main = function(SEED = NULL, PerfMin = NULL, PerfMax = NULL, Cutoff = NULL, 
                Weighting = NULL, dimx = NULL, dimlength = NULL,
                R = NULL, numobs = NULL, numapps = NULL, Trials = NULL, 
                SigBox = NULL, Results = NULL, varylab = NULL)
{
for(j in (PerfMin+dimlength):PerfMax)
{  
assign(varylab,j)
for(i in 1:Trials)
{
SEED = SEED + 1
performance.mat = performance(PerfMin = PerfMin,PerfMax = PerfMax, dimx = dimx, 
                    dimlength = dimlength, R = R, numobs = numobs, SEED = SEED)
eff <- eff.frontier(returns = get.returns(Xrescaled = 
                    matrix(as.numeric(performance.mat[1:(numobs*dimlength)]),
                    byrow = FALSE, nrow = numobs, ncol = dimlength, dimnames = 
                    list(c(1:numobs),dimx=dimx)),
                    dimlength = dimlength),dimx = dimx, dimlength = dimlength)
# print(graph.eff(eff = eff)) #diagnostic/graphic
apps.mat = applicants(PerfMin= PerfMin, PerfMax =PerfMax, dimx = dimx, dimlength = 
                      dimlength, R = R, numapps = numapps, SEED = SEED)
if(is.null(Weighting) == TRUE) Weighting = rbind(UnitWeights = 
                  as.data.frame(unit.model(dimlength = dimlength, dimx = dimx)), 
                  RegWeights = regression.model(dimx = dimx, dimlength = dimlength, 
                  newX = as.data.frame(do.call(cbind, performance.mat[
                  (dimlength*numobs + 1):(numobs*dimlength +dimlength)]))), 
                  MPT.model(eff = eff, dimx = dimx, dimlength = dimlength))
if(is.null(Weighting) == FALSE) Weighting = rbind(Weighting,(rbind(UnitWeights = 
                  as.data.frame(unit.model(dimlength = dimlength, dimx = dimx)), 
                  RegWeights = regression.model(dimx = dimx, dimlength = dimlength, 
                  newX = as.data.frame(do.call(cbind, performance.mat[
                  (dimlength*numobs + 1):(numobs*dimlength +dimlength)]))), 
                  MPT.model(eff = eff, dimlength = dimlength, dimx = dimx))))
#print(Weighting) #diagnostic
Signal = SignalDetection(Perf = cbind(matrix(as.numeric(apps.mat[1:(numapps*dimlength)]),
                  byrow = FALSE, nrow = numapps, ncol = dimlength, dimnames = 
                  list(c(1:numapps), dimx))[,1],t(as.matrix(Weighting[
                  (3*(i-1)+1):(3*(i-1)+3),])%*%t(as.matrix(matrix(as.numeric(apps.mat[
                  1:(numapps*dimlength)]), byrow = FALSE, nrow = numapps, ncol = 
                  dimlength, dimnames = list(c(1:numapps),dimx))[,2:dimlength])))),
                  PerfMin = PerfMin, PerfMax = PerfMax, Signal = NULL, Cutoff = 
                  Cutoff, numapps = numapps)
# print(Signal)
SigTot = cbind(table(factor(Signal[,1],lev = 1:4)),table(factor(Signal[,2],lev = 1:4)),
                  table(factor(Signal[,3],lev = 1:4)))
colnames(SigTot) = c("Unit Weighting","Multiple Regression Weighting", "MPT Weighting")
if(is.null(SigBox) == FALSE) SigBox = rbind(SigBox, SigTot)
if(is.null(SigBox) == TRUE) SigBox = SigTot
options("scipen"= 100, "digits"= 2)
}
CorrectSelections = colSums(SigBox[seq(1, nrow(SigBox), 4),])
CorrectRejections = colSums(SigBox[seq(2, nrow(SigBox), 4),])
FalsePositives = colSums(SigBox[seq(3, nrow(SigBox), 4),])
FalseNegatives = colSums(SigBox[seq(4, nrow(SigBox), 4),])
SigBox = rbind(CorrectSelections, CorrectRejections, FalsePositives, FalseNegatives)
rownames(SigBox) = c("Correct Selections","Correct Rejections","False Positives",
                     "False Negatives")
Sensitivity = SigBox[1,]/(SigBox[1,]+SigBox[4,])
Specificity = SigBox[2,]/(SigBox[2,]+SigBox[3,])
FalsePositiveRate = SigBox[3,]/(SigBox[1,]+SigBox[3,])
FalseNegativeRate = SigBox[4,]/(SigBox[4,]+SigBox[2,])
WeightingAVG = rbind(colMeans(Weighting[seq(1,ncol(Weighting),3),]),colMeans(Weighting[
                    seq(2,ncol(Weighting),3),]), colMeans(Weighting[
                    seq(3,ncol(Weighting),3),]))
SigBox = rbind(SigBox,Sensitivity, Specificity, FalsePositiveRate, 
                    FalseNegativeRate,t(WeightingAVG))
if(is.null(Results) == FALSE) Results = rbind(Results, SigBox)
if(is.null(Results) == TRUE) Results = SigBox
SigBox = NULL
Weighting = NULL
}
return(Results)
}
Results = Main(SEED = 8675309, PerfMin = 1, PerfMax = 100, Cutoff = 60, 
               Weighting = NULL, dimx = c("Performance", "Conscientiousness", 
              "Cognitive.Ability","Work.Sample"), dimlength = 4,
               R = matrix(cbind(  1,  .20,  .51, .57,
                                  .20,  1,  .01, .09,
                                  .51,  .01,  1, .34,
                                  .57,  .09,  .34, 1), nrow = 4,
              dimnames = list(c("Performance", "Conscientiousness", 
              "Cognitive.Ability","Work.Sample"), c("Performance", "Conscientiousness", 
              "Cognitive.Ability","Work.Sample"))),
               numobs = 100,   # number of observations
               numapps = 100,  #number of applicants
               Trials = 10,  # number of times to test each level
               SigBox = NULL, Results = NULL, varylab = "Cutoff")
varylab = "Cutoff"
matplot(Results[seq(1,nrow(Results),11),], type ="l", xlab = varylab, ylab = 
          "Correct Selections", main = "Correct Selections")
matplot(Results[seq(2,nrow(Results),11),], type ="l", xlab = varylab, ylab = 
          "Correct Rejections", main = "Correct Rejections")
matplot(Results[seq(3,nrow(Results),11),], type ="l", xlab = varylab, ylab = 
          "False Positives", main = "False Positives")
matplot(Results[seq(4,nrow(Results),11),], type ="l", xlab = varylab, ylab = 
          "False Negatives", main = "False Negatives")
matplot(Results[seq(5,nrow(Results),11),], type ="l", xlab = varylab, ylab = 
          "Sensitivity", main = "Sensitivity")
matplot(Results[seq(6,nrow(Results),11),], type ="l", xlab = varylab, ylab = 
          "Specificity", main = "Specificity")
matplot(Results[seq(7,nrow(Results),11),], type ="l", xlab = varylab, ylab = 
          "False Positive Rate", main = "False Positive Rate")
matplot(Results[seq(8,nrow(Results),11),], type ="l", xlab = varylab, ylab = 
          "False Negative Rate", main = "False Negative Rate")
matplot(Results[seq(9,nrow(Results),11),], type ="l", xlab = varylab, ylab = 
          "Weights", main = "Conscientiousness")
matplot(Results[seq(10,nrow(Results),11),], type ="l", xlab = varylab, ylab = 
          "Weights", main = "Cognitive Ability")
matplot(Results[seq(11,nrow(Results),11),], type ="l", xlab = varylab, ylab = 
          "Weights", main = "Work Sample")
##############Simulation 2###########################################################
Main = function(SEED = NULL, PerfMin = NULL, PerfMax = NULL, Cutoff = NULL, 
                Weighting = NULL, dimx = NULL, dimlength = NULL,
                R = NULL,numobsmax = NULL,numobs = NULL, numapps = NULL, 
                Trials = NULL, SigBox = NULL, Results = NULL, varylab = NULL)
{
  for(j in 10:numobsmax)
  {  
    numobs = j
    for(i in 1:Trials)
    {
      SEED = SEED + 1
      performance.mat = performance(PerfMin = PerfMin,PerfMax = PerfMax, dimx = dimx, 
                                    dimlength = dimlength, R = R, numobs = numobs, 
                                    SEED = SEED)
      eff <- eff.frontier(returns = get.returns(Xrescaled = 
                                    matrix(as.numeric(performance.mat[1:(numobs*dimlength)]),
                                    byrow = FALSE, nrow = numobs, ncol = dimlength, 
                                    dimnames = list(c(1:numobs),dimx=dimx)), dimlength = 
                                    dimlength),dimx = dimx, dimlength = dimlength)
      # print(graph.eff(eff = eff)) #diagnostic/graphic
      apps.mat = applicants(PerfMin= PerfMin, PerfMax =PerfMax, dimx = dimx, 
                            dimlength = dimlength, R, numapps, SEED)
      if(is.null(Weighting) == TRUE) Weighting = 
        rbind(UnitWeights = as.data.frame(unit.model(dimlength = dimlength, dimx = dimx)), 
        RegWeights = regression.model(dimx = dimx, dimlength = dimlength, 
        newX = as.data.frame(do.call(cbind, 
        performance.mat[(dimlength*numobs + 1):(numobs*dimlength +dimlength)]))), 
        MPT.model(eff = eff, dimx = dimx, dimlength = dimlength))
      if(is.null(Weighting) == FALSE) Weighting = rbind(Weighting,(rbind(UnitWeights = 
        as.data.frame(unit.model(dimlength = dimlength, dimx = dimx)), 
        RegWeights = regression.model(dimx = dimx, dimlength = dimlength, 
        newX = as.data.frame(do.call(cbind, 
        performance.mat[(dimlength*numobs + 1):(numobs*dimlength +dimlength)]))), 
        MPT.model(eff = eff, dimlength = dimlength, dimx = dimx))))
      #print(Weighting) #diagnostic
      Signal = SignalDetection(Perf = 
        cbind(matrix(as.numeric(apps.mat[1:(numapps*dimlength)]),
        byrow = FALSE, nrow = numapps, ncol = dimlength, dimnames = list(c(1:numapps),
        dimx))[,1],t(as.matrix(Weighting[(3*(i-1)+1):(3*(i-1)+3),])%*%
        t(as.matrix(matrix(as.numeric(apps.mat[1:(numapps*dimlength)]),                                                                                                                                                                                                 
        byrow = FALSE, nrow = numapps, ncol = dimlength, 
        dimnames = list(c(1:numapps),dimx))[,2:dimlength])))),
        PerfMin = PerfMin, PerfMax = PerfMax, Signal = NULL, Cutoff = Cutoff, 
        numapps = numapps)
      # print(Signal)
      SigTot = cbind(table(factor(Signal[,1],lev = 1:4)),table(factor(Signal[,2],
                           lev = 1:4)),table(factor(Signal[,3],lev = 1:4)))
      colnames(SigTot) = c("Unit Weighting","Multiple Regression Weighting", 
                           "MPT Weighting")
      if(is.null(SigBox) == FALSE) SigBox = rbind(SigBox, SigTot)
      if(is.null(SigBox) == TRUE) SigBox = SigTot
      options("scipen"= 100, "digits"= 2)
    }
    CorrectSelections = colSums(SigBox[seq(1, nrow(SigBox), 4),])
    CorrectRejections = colSums(SigBox[seq(2, nrow(SigBox), 4),])
    FalsePositives = colSums(SigBox[seq(3, nrow(SigBox), 4),])
    FalseNegatives = colSums(SigBox[seq(4, nrow(SigBox), 4),])
    SigBox = rbind(CorrectSelections, CorrectRejections, FalsePositives, 
                   FalseNegatives)
    rownames(SigBox) = c("Correct Selections","Correct Rejections","False Positives",
                         "False Negatives")
    Sensitivity = SigBox[1,]/(SigBox[1,]+SigBox[4,])
    Specificity = SigBox[2,]/(SigBox[2,]+SigBox[3,])
    FalsePositiveRate = SigBox[3,]/(SigBox[1,]+SigBox[3,])
    FalseNegativeRate = SigBox[4,]/(SigBox[4,]+SigBox[2,])
    WeightingAVG = rbind(colMeans(Weighting[seq(1,ncol(Weighting),3),]),
                         colMeans(Weighting[seq(2,ncol(Weighting),3),]),
                         colMeans(Weighting[seq(3,ncol(Weighting),3),]))
    SigBox = rbind(SigBox,Sensitivity, Specificity, FalsePositiveRate, 
                   FalseNegativeRate,t(WeightingAVG))
    if(is.null(Results) == FALSE) Results = rbind(Results, SigBox)
    if(is.null(Results) == TRUE) Results = SigBox
    SigBox = NULL
    Weighting = NULL
  }
  return(Results)
}
Results3 = Main(SEED = 8675309, PerfMin = 1, PerfMax = 100, Cutoff = 60, 
               Weighting = NULL, dimx = c("Performance", "Conscientiousness", 
               "Cognitive.Ability","Work.Sample"), dimlength = 4,
               R = matrix(cbind(  1,  .20,  .51, .57,
                                  .20,  1,  .01, .09,
                                  .51,  .01,  1, .34,
                                  .57,  .09,  .34, 1), nrow = 4,
               dimnames = list(c("Performance", "Conscientiousness", 
               "Cognitive.Ability","Work.Sample"), c("Performance", "Conscientiousness", 
               "Cognitive.Ability","Work.Sample"))),
               numobs = 300,   # number of observations
               numapps = 100,  #number of applicants
               Trials = 10000,  # number of times to test each level
               SigBox = NULL, Results = NULL, varylab = "numobs", numobsmax = 300)
# print(SigBox)
# print(Weighting)

varylab = "Number of Cases Used to Generate Weights"
matplot(Results3[seq(1,nrow(Results3),11),], type ="l", xlab = varylab, ylab = 
          "Correct Selections", main = "Correct Selections")
matplot(Results3[seq(2,nrow(Results3),11),], type ="l", xlab = varylab, ylab = 
          "Correct Rejections", main = "Correct Rejections")
matplot(Results3[seq(3,nrow(Results3),11),], type ="l", xlab = varylab, ylab = 
          "False Positives", main = "False Positives")
matplot(Results3[seq(4,nrow(Results3),11),], type ="l", xlab = varylab, ylab = 
          "False Negatives", main = "False Negatives")
matplot(Results3[seq(5,nrow(Results3),11),], type ="l", xlab = varylab, ylab = 
          "Sensitivity", main = "Sensitivity")
matplot(Results3[seq(6,nrow(Results3),11),], type ="l", xlab = varylab, ylab = 
          "Specificity", main = "Specificity")
matplot(Results3[seq(7,nrow(Results3),11),], type ="l", xlab = varylab, ylab = 
          "False Positive Rate", main = "False Positive Rate")
matplot(Results3[seq(8,nrow(Results3),11),], type ="l", xlab = varylab, ylab = 
          "False Negative Rate", main = "False Negative Rate")
matplot(Results3[seq(9,nrow(Results3),11),], type ="l", xlab = varylab, ylab = 
          "Weights", main = "Conscientiousness")
matplot(Results3[seq(10,nrow(Results3),11),], type ="l", xlab = varylab, ylab = 
          "Weights", main = "Cognitive Ability")
matplot(Results3[seq(11,nrow(Results3),11),], type ="l", xlab = varylab, ylab = 
          "Weights", main = "Work Sample")
############Simulation 3###########################################################
Main = function(SEED = NULL, PerfMin = NULL, PerfMax = NULL, Cutoff = NULL, 
                Weighting = NULL, dimx = NULL, dimlength = NULL,
                R = NULL,numobsmax = NULL,numobs = NULL, numapps = NULL, Trials = NULL, 
                SigBox = NULL, Results = NULL, varylab = NULL)
{
    for(i in 1:Trials)
    {
      SEED = SEED + 1
      performance.mat = performance(PerfMin = PerfMin,PerfMax = PerfMax, dimx = dimx, 
                        dimlength = dimlength, 
                        R = R, numobs = numobs, SEED = SEED)
      eff <- eff.frontier(returns = get.returns(Xrescaled = 
                        matrix(as.numeric(performance.mat[1:(numobs*dimlength)]),
                        byrow = FALSE, nrow = numobs, ncol = dimlength, dimnames = 
                        list(c(1:numobs),dimx=dimx)), dimlength = dimlength),dimx = 
                        dimx, dimlength = dimlength)
      # print(graph.eff(eff = eff)) #diagnostic/graphic
      apps.mat = applicants(PerfMin= PerfMin, PerfMax =PerfMax, dimx = dimx, dimlength = 
                        dimlength, R, numapps, SEED)
      if(is.null(Weighting) == TRUE) Weighting = rbind(UnitWeights = 
                  as.data.frame(unit.model(dimlength = dimlength, dimx = dimx)), 
                  RegWeights = regression.model(dimx = dimx, dimlength = dimlength, 
                  newX = as.data.frame(do.call(cbind, 
                  performance.mat[(dimlength*numobs + 1):(numobs*dimlength +dimlength)]))), 
                  MPT.model(eff = eff, dimx = dimx, dimlength = dimlength))
      if(is.null(Weighting) == FALSE) Weighting = rbind(Weighting,(rbind(UnitWeights = 
                  as.data.frame(unit.model(dimlength = dimlength, dimx = dimx)), 
                  RegWeights = regression.model(dimx = dimx, dimlength = dimlength, 
                  newX = as.data.frame(do.call(cbind, performance.mat[
                  (dimlength*numobs + 1):(numobs*dimlength +dimlength)]))), 
                  MPT.model(eff = eff, dimlength = dimlength, dimx = dimx))))
      #print(Weighting) #diagnostic
      Signal = SignalDetection(Perf = cbind(matrix(as.numeric(apps.mat[
                  1:(numapps*dimlength)]), byrow = FALSE, nrow = numapps, 
                  ncol = dimlength, dimnames = list(c(1:numapps), dimx))[,1],
                  t(as.matrix(Weighting[(3*(i-1)+1):(3*(i-1)+3),])%*%
                  t(as.matrix(matrix(as.numeric(apps.mat[1:(numapps*dimlength)]),                                                                                                                                                                                                 
                  byrow = FALSE, nrow = numapps, ncol = dimlength, dimnames = 
                  list(c(1:numapps),dimx))[,2:dimlength])))),
                  PerfMin = PerfMin, PerfMax = PerfMax, Signal = NULL, Cutoff = 
                  Cutoff, numapps = numapps)
      # print(Signal)
      SigTot = cbind(table(factor(Signal[,1],lev = 1:4)),table(factor(Signal[,2],
                  lev = 1:4)), table(factor(Signal[,3],lev = 1:4)))
      colnames(SigTot) = c("Unit Weighting","Multiple Regression Weighting", 
                  "MPT Weighting")
      if(is.null(SigBox) == FALSE) SigBox = rbind(SigBox, SigTot)
      if(is.null(SigBox) == TRUE) SigBox = SigTot
      options("scipen"= 100, "digits"= 2)
    }
    CorrectSelections = colSums(SigBox[seq(1, nrow(SigBox), 4),])
    CorrectRejections = colSums(SigBox[seq(2, nrow(SigBox), 4),])
    FalsePositives = colSums(SigBox[seq(3, nrow(SigBox), 4),])
    FalseNegatives = colSums(SigBox[seq(4, nrow(SigBox), 4),])
    SigBox = rbind(CorrectSelections, CorrectRejections, FalsePositives, 
                   FalseNegatives)
    rownames(SigBox) = c("Correct Selections","Correct Rejections","False Positives",
                         "False Negatives")
    Sensitivity = SigBox[1,]/(SigBox[1,]+SigBox[4,])
    Specificity = SigBox[2,]/(SigBox[2,]+SigBox[3,])
    FalsePositiveRate = SigBox[3,]/(SigBox[1,]+SigBox[3,])
    FalseNegativeRate = SigBox[4,]/(SigBox[4,]+SigBox[2,])
    WeightingAVG = rbind(colMeans(Weighting[seq(1,ncol(Weighting),3),]),
                         colMeans(Weighting[seq(2,ncol(Weighting),3),]),
                         colMeans(Weighting[seq(3,ncol(Weighting),3),]))
    SigBox = rbind(SigBox,Sensitivity, Specificity, FalsePositiveRate, 
                   FalseNegativeRate,t(WeightingAVG))
    if(is.null(Results) == FALSE) Results = rbind(Results, SigBox)
    if(is.null(Results) == TRUE) Results = SigBox
    SigBox = NULL
    Weighting = NULL
  return(Results)
}
Results3 = Main(SEED = 8675309, PerfMin = 1, PerfMax = 100, Cutoff = 60, 
                Weighting = NULL, dimx = c("Performance", "Conscientiousness", 
                "Cognitive.Ability","Work.Sample"), dimlength = 4,
                R = matrix(cbind(  1,  .20,  .51, .57,
                                   .20,  1,  .01, .09,
                                   .51,  .01,  1, .34,
                                   .57,  .09,  .34, 1), nrow = 4,
                dimnames = list(c("Performance", "Conscientiousness", "Cognitive.Ability",
                "Work.Sample"), c("Performance", "Conscientiousness", "Cognitive.Ability",
                "Work.Sample"))), numobs = 100,   # number of observations
                numapps = 100,  #number of applicants
                Trials = 10000,  # number of times to test each level
                SigBox = NULL, Results = NULL, varylab = "numobs", numobsmax = 300)
print(Results3)
################Simulation 4#######################################################
Main = function(SEED = NULL, PerfMin = NULL, PerfMax = NULL, Cutoff = NULL, 
                Weighting = NULL, dimx = NULL, dimlength = NULL, R = NULL,
                numobsmax = NULL, numobs = NULL, numapps = NULL, Trials = NULL, 
                SigBox = NULL, Results = NULL, varylab = NULL)
{
  for(i in 1:Trials)
  {
    SEED = SEED + 1
    performance.mat = performance(PerfMin = PerfMin,PerfMax = PerfMax, dimx = dimx, 
                                  dimlength = dimlength, 
                                  R = R, numobs = numobs, SEED = SEED)
    eff <- eff.frontier(returns = get.returns(Xrescaled = 
                                  matrix(as.numeric(performance.mat[1:(numobs*dimlength)]),
                                  byrow = FALSE, nrow = numobs, ncol = dimlength, 
                                  dimnames = list(c(1:numobs),dimx=dimx)),
                                  dimlength = dimlength),dimx = dimx, dimlength = dimlength)
    # print(graph.eff(eff = eff)) #diagnostic/graphic
    apps.mat = applicants(PerfMin= PerfMin, PerfMax =PerfMax, dimx = dimx, 
                          dimlength = dimlength, R, numapps, SEED)
    if(is.null(Weighting) == TRUE) Weighting = rbind(UnitWeights = 
                          as.data.frame(unit.model(dimlength = dimlength, dimx = dimx)), 
                          RegWeights = regression.model(dimx = dimx, dimlength = dimlength, 
                          newX = as.data.frame(do.call(cbind, performance.mat[
                          (dimlength*numobs + 1):(numobs*dimlength +dimlength)]))), 
                          MPT.model(eff = eff, dimx = dimx, dimlength = dimlength))
    if(is.null(Weighting) == FALSE) Weighting = rbind(Weighting,(rbind(UnitWeights = 
                          as.data.frame(unit.model(dimlength = dimlength, dimx = dimx)), 
                          RegWeights = regression.model(dimx = dimx, dimlength = dimlength, 
                          newX = as.data.frame(do.call(cbind, performance.mat[
                          (dimlength*numobs + 1):(numobs*dimlength +dimlength)]))), 
                          MPT.model(eff = eff, dimlength = dimlength, dimx = dimx))))
    #print(Weighting) #diagnostic
    Signal = SignalDetection(Perf = cbind(matrix(as.numeric(apps.mat[1:(numapps*dimlength)]),
                          byrow = FALSE, nrow = numapps, ncol = dimlength, 
                          dimnames = list(c(1:numapps),dimx))[,1],
                          t(as.matrix(Weighting[(3*(i-1)+1):(3*(i-1)+3),])%*%
                          t(as.matrix(matrix(as.numeric(apps.mat[1:(numapps*dimlength)]),                                                                                                                                                                                                 
                          byrow = FALSE, nrow = numapps, ncol = dimlength, dimnames = 
                          list(c(1:numapps),dimx))[,2:dimlength])))),
                          PerfMin = PerfMin, PerfMax = PerfMax, 
                          Signal = NULL, Cutoff = Cutoff, numapps = numapps)
    # print(Signal)
    SigTot = cbind(table(factor(Signal[,1],lev = 1:4)),table(factor(Signal[,2],lev = 1:4)),
                          table(factor(Signal[,3],lev = 1:4)))
    colnames(SigTot) = c("Unit Weighting","Multiple Regression Weighting", "MPT Weighting")
    if(is.null(SigBox) == FALSE) SigBox = rbind(SigBox, SigTot)
    if(is.null(SigBox) == TRUE) SigBox = SigTot
    options("scipen"= 100, "digits"= 2)
  }
  CorrectSelections = colSums(SigBox[seq(1, nrow(SigBox), 4),])
  CorrectRejections = colSums(SigBox[seq(2, nrow(SigBox), 4),])
  FalsePositives = colSums(SigBox[seq(3, nrow(SigBox), 4),])
  FalseNegatives = colSums(SigBox[seq(4, nrow(SigBox), 4),])
  SigBox = rbind(CorrectSelections, CorrectRejections, FalsePositives, FalseNegatives)
  rownames(SigBox) = c("Correct Selections","Correct Rejections","False Positives",
                       "False Negatives")
  Sensitivity = SigBox[1,]/(SigBox[1,]+SigBox[4,])
  Specificity = SigBox[2,]/(SigBox[2,]+SigBox[3,])
  FalsePositiveRate = SigBox[3,]/(SigBox[1,]+SigBox[3,])
  FalseNegativeRate = SigBox[4,]/(SigBox[4,]+SigBox[2,])
  WeightingAVG = rbind(colMeans(Weighting[seq(1,ncol(Weighting),3),]),
                       colMeans(Weighting[seq(2,ncol(Weighting),3),]),
                       colMeans(Weighting[seq(3,ncol(Weighting),3),]))
  SigBox = rbind(SigBox,Sensitivity, Specificity, FalsePositiveRate, 
                 FalseNegativeRate,t(WeightingAVG))
  if(is.null(Results) == FALSE) Results = rbind(Results, SigBox)
  if(is.null(Results) == TRUE) Results = SigBox
  SigBox = NULL
  Weighting = NULL
  return(Results)
}
Results4 = Main(SEED = 8675309, PerfMin = 1, PerfMax = 100, Cutoff = 60, 
                Weighting = NULL, dimx = c("Performance", "Conscientiousness", 
                "Cognitive.Ability","Work.Sample"), dimlength = 4,
                R = matrix(cbind(  1,  .20,  .51, .34, 
                                 .20,    1,  .01,-.25, 
                                 .51,  .01,    1, .22,
                                 .34, -.25,  .22,   1), nrow=4,
                dimnames = list(c("Performance", "Conscientiousness", 
                "Cognitive.Ability","Risk.Taking"), c("Performance", "Conscientiousness", 
                "Cognitive.Ability","Risk.Taking"))),
                numobs = 100,   # number of observations
                numapps = 100,  #number of applicants
                Trials = 10000,  # number of times to test each level
                SigBox = NULL, Results = NULL, varylab = "numobs", numobsmax = 300)
print(Results4)