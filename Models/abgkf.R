library(ggplot2)
library(plyr)
library(reshape2)
library(magrittr)
library(foreach)
library(readxl)
library(doParallel)
registerDoParallel(cores=19)
library(pomp)
library(panelPomp)
stopifnot(packageVersion("pomp") >= "1.4.7")

Mesocosm_data = read_excel("Mesocosm data.xls")
dentNoPara <- Mesocosm_data[1:100, ]
dentNoPara <- subset(dentNoPara, select = c(rep, day, dentadult))
dentNoPara <- dentNoPara[100: 1, ]
dentNoPara$day = (dentNoPara$day - 1) * 5 + 7
head(dentNoPara)

data = list()
trails = c("A","B","C","D","E","F","G","H","I","J")
for (i in 1: length(trails)){
  data[[i]]=subset(dentNoPara, select = c("day", "dentadult"),
                   dentNoPara$rep == trails[i])
}

dent_rpro <- Csnippet("
                      double Fbirth, Fpred, Sgrow, Sdeath, Sspl;
                      double noiB, noiG ;
                      double delta = 0.013;
                      
                      noiB = rnorm(0, sigB * sqrt(dt));
                      noiG = rnorm(0, sigG * sqrt(dt));
                      
                      
                      Fbirth = F * (alpha + noiB) * (1 - F/kf) * dt;
                      Fpred =  S * F * (Beta * dt);
                      Sgrow = S * F * (Theta * dt + noiG);
                      Sdeath = S * (gamma * dt);
                      Sspl = S * (delta * dt);
                      
                      F += Fbirth - Fpred;
                      S += Sgrow - Sdeath  - Sspl;
                      
                      if (S <= 1.0) {
                      S = 1.0;
                      error_count += 1;
                      }
                      if (F <= 0) {
                      F = 1000.0;
                      error_count += 1;
                      }
                      ")



# Assume t0 = day 4
dent_init = Csnippet("
                     S = 3; //3= 45/15L
                     F =  1.667e3;
                     error_count = 0;
                     ")



dmeas = Csnippet("
                 double tol = 1.0e-18;
                 double f = 0.0;
                 double delta = 0.013;
                 
                 if (error_count > 0.0) {
                 lik = tol;
                 } else {
                 f += dbinom(dentadult, nearbyint((1/delta)*S), delta, 1) + tol;
                 lik = (give_log) ? f : exp(f);
                 }
                 ")

rmeas = Csnippet("
                 double delta = 0.013;
                 dentadult = rbinom(nearbyint((1/delta)*S),delta);
                 ")

dent_fromEstimationScale <- Csnippet("
                                     Talpha = exp(alpha);
                                     TBeta = exp(Beta);
                                     TTheta = exp(Theta);
                                     Tgamma = exp(gamma);
                                     TsigB = exp(sigB);
                                     TsigG = exp(sigG);
                                     Tkf = exp(kf);
                                     ")



dent_ToEstimationScale <- Csnippet("
                                   Talpha =log(alpha);
                                   TBeta = log(Beta);
                                   TTheta = log(Theta);
                                   Tgamma = log(gamma);
                                   TsigB = log(sigB); 
                                   TsigG = log(sigG);
                                   Tkf = log(kf);
                                   ")

pomplist=list()
parameters = c(0.875100208, 0.0041770577,0.004621913,7.991700e-04, 303.5295,0.0007095320,
               8.143122e-06)
names(parameters) =  c("alpha","Beta", "gamma", "Theta","kf","sigB","sigG")

for (i in 1:10){
  pomp(data = data[[i]], 
       time="day",t0=1,
       rprocess=euler.sim(dent_rpro,delta.t=1/4),
       initializer=dent_init,
       dmeasure = dmeas, 
       rmeasure = rmeas,
       fromEstimationScale=dent_fromEstimationScale,
       toEstimationScale=dent_ToEstimationScale,
       obsnames = "dentadult",
       zeronames = "error_count",
       paramnames = c("alpha","Beta", "Theta", "gamma",  "sigG", "sigB", "kf"), 
       statenames = c("S", "error_count", "F")
  ) -> pomplist[[i]]
  coef(pomplist[[i]])=parameters
}
names(pomplist)=paste("u", 1:10,sep = "")

#constuct panelPOMP
shared_parameter =c(7.991700e-04,0.0007095320,8.143122e-06)
names(shared_parameter) = c("Theta", "sigB","sigG")

specific_mat = matrix(data = c(rep(0.875100208, 10), rep( 0.0041770577,10), rep(0.004621913,10),
                               rep(303.5295,10)), nrow = 4, nc = 10, byrow = T, 
dim = list(rownames = c("alpha","Beta", "gamma", "kf"),
           colnames=c("u1","u2","u3","u4","u5", "u6","u7","u8","u9","u10")))
panelfood = panelPomp(pomplist, shared=shared_parameter,specific = specific_mat)

run_level <- 3
algorithmic.params <- list(
  Np =     c(50, 1e4, 1e5),
  Np_rep = c( 2,  10,  10),
  Mp =     c(50, 1e4,30e3),
  Nmif =   c( 2,  200, 250)
)

names.of.estimated.shared.params <-  c("Theta", "sigB","sigG")
names.of.estimated.specific.params <- c("alpha","Beta", "gamma","kf")
dent_rw.sd= 0.00125
parameter_candidates = list(shared_parameter, specific_mat)
names(parameter_candidates)=c("shared", "specific")
U = length(panelfood)

bake(file = "abgkf.rds",
     seed = 334388458L,
     kind = "L'Ecuyer" ,
     {
       foreach(
         i = 1:getDoParWorkers(),
         .packages = c("pomp", "panelPomp"),
         .inorder = FALSE,
         .options.multicore = list(set.seed = TRUE)
       ) %dopar%
       {
         guessed.parameter.values <- parameter_candidates
         if (!is.null(names.of.estimated.shared.params)) {
           guessed.parameter.values$shared[names.of.estimated.shared.params] <-
             rlnorm(
               n = length(names.of.estimated.shared.params),
               meanlog = log(guessed.parameter.values$shared[names.of.estimated.shared.params]),
               sdlog = 1
             )
         }
         if (!is.null(names.of.estimated.specific.params)) {
           for (i.u in 1:U) {
             guessed.parameter.values$specific[names.of.estimated.specific.params, i.u] <-
               rlnorm(
                 n = length(names.of.estimated.specific.params),
                 meanlog = log(guessed.parameter.values$specific[names.of.estimated.specific.params, i.u]),
                 sdlog = 1
               )
           }
         }
         mif2(
           panelfood,
           Nmif = algorithmic.params$Nmif[run_level],
           shared.start = guessed.parameter.values$shared,
           specific.start = coef(panelfood)$specific,
           transform = TRUE,
           rw.sd = rw.sd(Beta=dent_rw.sd,
                         alpha=dent_rw.sd,
                         sigG=dent_rw.sd,
                         sigB=dent_rw.sd,
                         gamma=dent_rw.sd,
                         Theta=dent_rw.sd,
                         kf = dent_rw.sd),
           cooling.type = "geometric",
           cooling.fraction.50 = 0.5,
           Np = algorithmic.params$Mp[run_level]
         ) -> m1
         
         ll <- replicate(n = algorithmic.params$Np_rep[run_level],
                         unitlogLik(pfilter(m1,
                                            Np = algorithmic.params$Np[run_level])))
         list(mif = m1,
              ll = panel_logmeanexp(x = ll,
                                    MARGIN = 1,
                                    se = TRUE))
       }
     }) -> mf

lls <- sapply(mf, getElement, "ll")
best <- which.max(sapply(mf, getElement, "ll")[1,])
mif.estimate <- coef(mf[[best]]$mif)
pf.loglik.of.mif.estimate <- unname(mf[[best]]$ll[1])
s.e.of.pf.loglik.of.mif.estimate <- unname(mf[[best]]$ll[2])

save(mf,lls,best,mif.estimate,pf.loglik.of.mif.estimate, s.e.of.pf.loglik.of.mif.estimate, 
     file = "abgkf.RData")


