# main algorithm to fit the LTCA model

pointest <- function(dat,num_class,id,time,num_obs,features,Y_dist,covx,ipw,stop,tol=0.005,max=50,varest,balanced=T,verbose=T){
  # let tau0 be a matrix
  #require(VGAM)
  #require(geepack)
  requireNamespace("VGAM")
  requireNamespace("geepack")

  n=length(unique(dat[,id]))

  # sort data by id and time
  dat <- dat[order(dat[,id],dat[,time]),]

  num_feature <- length(Y_dist)
  num_covx <- length(covx)

  covx_lb <- which(colnames(dat) %in% covx)

  colnames(dat)[which(colnames(dat) %in% id)] <- "id"
  colnames(dat)[which(colnames(dat) %in% time)] <- "time"
  colnames(dat)[which(colnames(dat) %in% num_obs)] <- "num_obs"
  dat$id <- as.numeric(as.factor(dat$id))
  colnames(dat)[which(colnames(dat) %in% features)] <- paste("y.",1:num_feature,sep="")
  colnames(dat)[which(colnames(dat) %in% covx)] = paste("baselinecov.",1:num_covx,sep="")

  # extract baseline data
  baseline <- dat[match(unique(dat$id), dat$id),]

  # replicate baseline data to create an observation for each latent class
  # this is for the pseudo observations used in multinomial regression
  lab_class = rep(1,nrow(baseline)*num_class)
  for (i in 2:num_class){
    lab_class[((i-1)*nrow(baseline)+1):((i)*nrow(baseline))] = i
  }
  # initialization scheme: randomly assign latent class memberships to
  # 25*num_class/n subjects

  tau0 <- matrix(0,nrow=n,ncol=num_class)
  prop <- 100*num_class/n
  u <- runif(n,0,1)
  for (i in 1:n){
    if (u[i] < prop){
      tau0[i,ceiling(num_class*runif(1))] = 1
    }
  }
  tau0[tau0 < 1e-8] = 1e-8

  ##### Fit the GEE models #####

  # initialize paramter estimates
  obj_beta0 <- matrix(0,ncol=num_class,nrow=num_feature)
  obj_beta1 <- matrix(0,ncol=num_class,nrow=num_feature)
  obj_phi <- matrix(0,ncol=num_class,nrow=num_feature)
  obj_gamma <- matrix(0,ncol=num_class,nrow=num_feature)

  # intiailize qic
  eqic = 0

  for (c in 1:num_class){
    tau = rep(0,nrow(dat))

    # assign posterior membership probability w.r.t class c to the data
    for (i in 1:n){
      tau[dat$id==i] = tau0[i,c]
    }

    for (j in 1:num_feature){

      # yy is the jth longitudinal marker
      yy <- as.numeric(dat[,paste('y.',j,sep='')])
      nalabel <- is.na(yy)
      tau <- tau[!nalabel]

      # fit corresponding GEE model with weights tau and AR1 correlation structure
      if(Y_dist[j] == 'normal'){
        geefit <- geepack::geeglm(yy[!nalabel]~time,family=gaussian,id=id,waves=num_obs,data=dat[!nalabel,],weight=tau*ipw,corstr='ar1')#,scale.fix=T,scale.value = 1)
      }else if (Y_dist[j] == 'poi'){
        geefit <- geepack::geeglm(yy[!nalabel]~time,family=poisson,id=id,waves=num_obs,data=dat[!nalabel,],weight=tau*ipw,corstr='ar1')#,scale.fix=T,scale.value = 1)
      }else if (Y_dist[j] == 'bin'){
        geefit <- geepack::geeglm(((yy[!nalabel]))~time,family=binomial,id=id,waves=num_obs,data=dat[!nalabel,],weight=tau*ipw,corstr='ar1')#,scale.fix=T,scale.value = 1)
      }

      # obtain point estimates
      obj_beta0[j,c] <- coef(geefit)[1]
      obj_beta1[j,c] <- coef(geefit)[2]
      obj_phi[j,c] <- geefit$geese$gamma
      obj_gamma[j,c] <- geefit$geese$alpha

      # use ic function to obtain eqic
      eqic = eqic + ic(geefit,Y_dist[j],yy[!nalabel])
    }

    # for each class c, compute the class-specific entropy
    # and count it to ceeqic in the iteration w.r.t c

    entropy = sum(weights(geefit) * log(weights(geefit)))
    ceeqic = eqic - 2*entropy
  }

  # finalize ICs by adding number of parameters
  eqica = eqic + 2*(4*num_class*num_feature+(num_class-1)*2) # here (num_class-1)*2 is the number of parameters in the latent class model.
  eqicb = eqic + (4*num_class*num_feature+(num_class-1)*2)*log(n)
  ceeqic = ceeqic + (4*num_class*num_feature+(num_class-1)*2)*log(n)

  ######################################################################

  # prepare a matrix y to store markers

  y <- dat$y.1
  for (j in 2:num_feature){
    y <- cbind(y,dat[,paste('y.',j,sep='')])
  }

  # obtain approximated likelihood ratios exp(w)

  ew <- LinProj(obj_beta0,obj_beta1,obj_phi,obj_gamma,dat,y,Y_dist,balanced)

  # restrict the upper and lower bound of ew

  ew[ew>1e2] = 1e2
  ew[ew<1e-2] = 1e-2

  ###### Fit the multinomial logistic regression model - estimate alpha ######

  # use vglm function

  # note that the data used here is aggregated by pseudo observations
  # with corresponding pseudo classes and posterior weights

  vars <- paste(colnames(baseline)[covx_lb],collapse="+")
  regression <- paste0("as.factor(class)", " ~ ", vars)
  obj_alpha_vec <- coef(VGAM::vglm(as.formula(regression),family = VGAM::multinomial(refLevel = 1),
                             weight=tau0,
                             data=data.frame(do.call("rbind", rep(list(baseline), num_class)),
                                             class = lab_class, tau0 = as.vector(tau0)
                             )
  )
  )

  # avoid obtaining 0

  obj_alpha_vec[abs(obj_alpha_vec) < 1e-8] = 1e-8

  # assign the results to a list called obj_alpha
  # obj_alpha[[i]] saves parameters associated with latent class (i+1)

  obj_alpha = list()
  for (i in 1:(num_class-1)){
    obj_alpha[[i]] = obj_alpha_vec[seq(from=i,to=(num_class-1)*num_covx+i,by=num_class-1)]
  }

  ###### Obtain posterior membership probability tau0 ######


  # p is the fitted membership probability by the multinomial model

  p <- matrix(0,ncol=num_class,nrow=n)
  p[,1] = 1
  for (i in 2:num_class){
    p[,i] = exp(cbind(1,as.matrix(baseline[,covx_lb])) %*% obj_alpha[[i-1]])
  }
  psum = rowSums(p)
  p = apply(p,2,function(x) x/psum)

  # tau0 is the posterior membership probability

  pew <- p*ew
  pewsum = rowSums(pew)
  tau0 = apply(pew,2,function(x) x/pewsum)
  tau0[tau0<1e-8] = 1e-8

  ###### set diff=10 and count=0 and start iterative algorithm ######

  diff = 10
  count=0

  while(diff > tol & count < max){

    count=count+1
    alpha_past <- obj_alpha
    beta0_past <- obj_beta0
    beta1_past <- obj_beta1
    tau_past <- tau0

    ###### repeating the previous GEE fitting steps ######

    # intiailize qic
    eqic = 0

    for (c in 1:num_class){
      tau = rep(0,nrow(dat))

      # assign posterior membership probability w.r.t class c to the data
      for (i in 1:n){
        tau[dat$id==i] = tau0[i,c]
      }

      for (j in 1:num_feature){

        # yy is the jth longitudinal marker
        yy <- as.numeric(dat[,paste('y.',j,sep='')])
        nalabel <- is.na(yy)
        tau <- tau[!nalabel]

        # fit corresponding GEE model with weights tau and AR1 correlation structure
        if(Y_dist[j] == 'normal'){
          geefit <- geepack::geeglm(yy[!nalabel]~time,family=gaussian,id=id,waves=num_obs,data=dat[!nalabel,],weight=tau*ipw,corstr='ar1')#,scale.fix=T,scale.value = 1)
        }else if (Y_dist[j] == 'poi'){
          geefit <- geepack::geeglm(yy[!nalabel]~time,family=poisson,id=id,waves=num_obs,data=dat[!nalabel,],weight=tau*ipw,corstr='ar1')#,scale.fix=T,scale.value = 1)
        }else if (Y_dist[j] == 'bin'){
          geefit <- geepack::geeglm(((yy[!nalabel]))~time,family=binomial,id=id,waves=num_obs,data=dat[!nalabel,],weight=tau*ipw,corstr='ar1')#,scale.fix=T,scale.value = 1)
        }

        # obtain point estimates
        obj_beta0[j,c] <- coef(geefit)[1]
        obj_beta1[j,c] <- coef(geefit)[2]
        obj_phi[j,c] <- geefit$geese$gamma
        obj_gamma[j,c] <- geefit$geese$alpha

        # use ic function to obtain eqic
        eqic = eqic + ic(geefit,Y_dist[j],yy[!nalabel])
      }

      # for each class c, compute the class-specific entropy
      # and count it to ceeqic in the iteration w.r.t c

      entropy = sum(weights(geefit) * log(weights(geefit)))
      ceeqic = eqic - 2*entropy
    }

    # finalize ICs by adding number of parameters
    eqica = eqic + 2*(4*num_class*num_feature+(num_class-1)*2) # here (num_class-1)*2 is the number of parameters in the latent class model.
    eqicb = eqic + (4*num_class*num_feature+(num_class-1)*2)*log(n)
    ceeqic = ceeqic + (4*num_class*num_feature+(num_class-1)*2)*log(n)

    ######################################################################

    # prepare a matrix y to store markers

    y <- dat$y.1
    for (j in 2:num_feature){
      y <- cbind(y,dat[,paste('y.',j,sep='')])
    }

    # obtain approximated likelihood ratios exp(w)

    #ew <- LinProj(obj_beta0,obj_beta1,obj_phi,obj_gamma,dat,y,Y_dist)
    ew <- LinProj(obj_beta0,obj_beta1,obj_phi,obj_gamma,dat,y,Y_dist,balanced)

    # restrict the upper and lower bound of ew

    ew[ew>1e2] = 1e2
    ew[ew<1e-2] = 1e-2

    ###### Fit the multinomial logistic regression model - estimate alpha ######

    # use vglm function

    # note that the data used here is aggregated by pseudo observations
    # with corresponding pseudo classes and posterior weights

    vars <- paste(colnames(baseline)[covx_lb],collapse="+")
    regression <- paste0("as.factor(class)", " ~ ", vars)
    obj_alpha_vec <- coef(VGAM::vglm(as.formula(regression),family = VGAM::multinomial(refLevel = 1),
                               weight=tau0,
                               data=data.frame(do.call("rbind", rep(list(baseline), num_class)),
                                               class = lab_class, tau0 = as.vector(tau0)
                               )
    )
    )

    # avoid obtaining 0

    obj_alpha_vec[abs(obj_alpha_vec) < 1e-8] = 1e-8

    # assign the results to a list called obj_alpha
    # obj_alpha[[i]] saves parameters associated with latent class (i+1)

    obj_alpha = list()
    for (i in 1:(num_class-1)){
      obj_alpha[[i]] = obj_alpha_vec[seq(from=i,to=(num_class-1)*num_covx+i,by=num_class-1)]
    }

    ###### Obtain posterior membership probability tau0 ######


    # p is the fitted membership probability by the multinomial model

    p <- matrix(0,ncol=num_class,nrow=n)
    p[,1] = 1
    for (i in 2:num_class){
      p[,i] = exp(cbind(1,as.matrix(baseline[,covx_lb])) %*% obj_alpha[[i-1]])
    }
    psum = rowSums(p)
    p = apply(p,2,function(x) x/psum)

    # tau0 is the posterior membership probability

    pew <- p*ew
    pewsum = rowSums(pew)
    tau0 = apply(pew,2,function(x) x/pewsum)
    tau0[tau0<1e-8] = 1e-8

    ##### compute covergence criteria,
    # which is the L2 norm of difference between tau and updated tau

    if(stop == 'tau'){
      diff = sum((tau0-tau_past)^2)
    }else if (stop == 'par'){
      diff = max(abs(unlist(obj_alpha)-unlist(alpha_past)),abs(obj_beta0-beta0_past),abs(obj_beta1-beta1_past))
    }else{
      print("Error: stop undefined.")
      break
    }

    if (verbose==T){
      cat(paste('iteration: ',count,'\n','convergence criteria: ',diff,'\n','CEEQIC: ',ceeqic,'\n',sep=''))
    }


  }


  ###### repeating the previous GEE fitting steps ######
  # intiailize qic
  eqic = 0

  for (c in 1:num_class){
    tau = rep(0,nrow(dat))

    # assign posterior membership probability w.r.t class c to the data
    for (i in 1:n){
      tau[dat$id==i] = tau0[i,c]
    }

    for (j in 1:num_feature){

      # yy is the jth longitudinal marker
      yy <- as.numeric(dat[,paste('y.',j,sep='')])
      nalabel <- is.na(yy)
      tau <- tau[!nalabel]

      # fit corresponding GEE model with weights tau and AR1 correlation structure
      if(Y_dist[j] == 'normal'){
        geefit <- geepack::geeglm(yy[!nalabel]~time,family=gaussian,id=id,waves=num_obs,data=dat[!nalabel,],weight=tau*ipw,corstr='ar1')#,scale.fix=T,scale.value = 1)
      }else if (Y_dist[j] == 'poi'){
        geefit <- geepack::geeglm(yy[!nalabel]~time,family=poisson,id=id,waves=num_obs,data=dat[!nalabel,],weight=tau*ipw,corstr='ar1')#,scale.fix=T,scale.value = 1)
      }else if (Y_dist[j] == 'bin'){
        geefit <- geepack::geeglm(((yy[!nalabel]))~time,family=binomial,id=id,waves=num_obs,data=dat[!nalabel,],weight=tau*ipw,corstr='ar1')#,scale.fix=T,scale.value = 1)
      }

      # obtain point estimates
      obj_beta0[j,c] <- coef(geefit)[1]
      obj_beta1[j,c] <- coef(geefit)[2]
      obj_phi[j,c] <- geefit$geese$gamma
      obj_gamma[j,c] <- geefit$geese$alpha

      # use ic function to obtain eqic
      eqic = eqic + ic(geefit,Y_dist[j],yy[!nalabel])
    }

    # for each class c, compute the class-specific entropy
    # and count it to ceeqic in the iteration w.r.t c

    entropy = sum(weights(geefit) * log(weights(geefit)))
    ceeqic = eqic - 2*entropy
  }

  ######################################################################

  # prepare a matrix y to store markers

  y <- dat$y.1
  for (j in 2:num_feature){
    y <- cbind(y,dat[,paste('y.',j,sep='')])
  }

  # obtain approximated likelihood ratios exp(w)

  #ew <- LinProj(obj_beta0,obj_beta1,obj_phi,obj_gamma,dat,y,Y_dist)
  ew <- LinProj(obj_beta0,obj_beta1,obj_phi,obj_gamma,dat,y,Y_dist,balanced)

  # restrict the upper and lower bound of ew

  ew[ew>1e3] = 1e3
  ew[ew<1e-3] = 1e-3

  ###### Fit the multinomial logistic regression model - estimate alpha ######

  # use vglm function

  # note that the data used here is aggregated by pseudo observations
  # with corresponding pseudo classes and posterior weights

  vars <- paste(colnames(baseline)[covx_lb],collapse="+")
  regression <- paste0("as.factor(class)", " ~ ", vars)
  obj_alpha_vec <- coef(VGAM::vglm(as.formula(regression),family = VGAM::multinomial(refLevel = 1),
                             weight=tau0,
                             data=data.frame(do.call("rbind", rep(list(baseline), num_class)),
                                             class = lab_class, tau0 = as.vector(tau0)
                             )
  )
  )

  # avoid obtaining 0

  obj_alpha_vec[abs(obj_alpha_vec) < 1e-8] = 1e-8

  # assign the results to a list called obj_alpha
  # obj_alpha[[i]] saves parameters associated with latent class (i+1)

  obj_alpha = list()
  for (i in 1:(num_class-1)){
    obj_alpha[[i]] =  obj_alpha_vec[seq(from=i,to=(num_class-1)*num_covx+i,by=num_class-1)]
  }

  ###### Obtain posterior membership probability tau0 ######


  # p is the fitted membership probability by the multinomial model

  p <- matrix(0,ncol=num_class,nrow=n)
  p[,1] = 1
  for (i in 2:num_class){
    p[,i] = exp(cbind(1,as.matrix(baseline[,covx_lb])) %*% obj_alpha[[i-1]])
  }
  psum = rowSums(p)
  p = apply(p,2,function(x) x/psum)

  # tau0 is the posterior membership probability

  pew <- p*ew
  pewsum = rowSums(pew)
  tau0 = apply(pew,2,function(x) x/pewsum)
  tau0[tau0<1e-8] = 1e-8


  ###### repeating the previous GEE fitting steps ######

  # intiailize qic
  eqic = 0

  for (c in 1:num_class){
    tau = rep(0,nrow(dat))

    # assign posterior membership probability w.r.t class c to the data
    for (i in 1:n){
      tau[dat$id==i] = tau0[i,c]
    }

    for (j in 1:num_feature){

      # yy is the jth longitudinal marker
      yy <- as.numeric(dat[,paste('y.',j,sep='')])
      nalabel <- is.na(yy)
      tau <- tau[!nalabel]

      # fit corresponding GEE model with weights tau and AR1 correlation structure
      if(Y_dist[j] == 'normal'){
        geefit <- geepack::geeglm(yy[!nalabel]~time,family=gaussian,id=id,waves=num_obs,data=dat[!nalabel,],weight=tau*ipw,corstr='ar1')#,scale.fix=T,scale.value = 1)
      }else if (Y_dist[j] == 'poi'){
        geefit <- geepack::geeglm(yy[!nalabel]~time,family=poisson,id=id,waves=num_obs,data=dat[!nalabel,],weight=tau*ipw,corstr='ar1')#,scale.fix=T,scale.value = 1)
      }else if (Y_dist[j] == 'bin'){
        geefit <- geepack::geeglm(((yy[!nalabel]))~time,family=binomial,id=id,waves=num_obs,data=dat[!nalabel,],weight=tau*ipw,corstr='ar1')#,scale.fix=T,scale.value = 1)
      }

      # obtain point estimates
      obj_beta0[j,c] <- coef(geefit)[1]
      obj_beta1[j,c] <- coef(geefit)[2]
      obj_phi[j,c] <- geefit$geese$gamma
      obj_gamma[j,c] <- geefit$geese$alpha

      # use ic function to obtain eqic
      eqic = eqic + ic(geefit,Y_dist[j],yy[!nalabel])
    }

    # for each class c, compute the class-specific entropy
    # and count it to ceeqic in the iteration w.r.t c

    entropy = sum(weights(geefit) * log(weights(geefit)))
    ceeqic = eqic - 2*entropy
  }

  ######################################################################

  # prepare a matrix y to store markers

  y <- dat$y.1
  for (j in 2:num_feature){
    y <- cbind(y,dat[,paste('y.',j,sep='')])
  }

  # obtain approximated likelihood ratios exp(w)

  #ew <- LinProj(obj_beta0,obj_beta1,obj_phi,obj_gamma,dat,y,Y_dist)
  ew <- LinProj(obj_beta0,obj_beta1,obj_phi,obj_gamma,dat,y,Y_dist,balanced)

  # restrict the upper and lower bound of ew

  ew[ew>1e3] = 1e3
  ew[ew<1e-3] = 1e-3

  ###### Fit the multinomial logistic regression model - estimate alpha ######

  # use vglm function

  # note that the data used here is aggregated by pseudo observations
  # with corresponding pseudo classes and posterior weights

  vars <- paste(colnames(baseline)[covx_lb],collapse="+")
  regression <- paste0("as.factor(class)", " ~ ", vars)
  obj_alpha_vec <- coef(VGAM::vglm(as.formula(regression),family = VGAM::multinomial(refLevel = 1),
                             weight=tau0,
                             data=data.frame(do.call("rbind", rep(list(baseline), num_class)),
                                             class = lab_class, tau0 = as.vector(tau0)
                             )
  )
  )

  # avoid obtaining 0

  obj_alpha_vec[abs(obj_alpha_vec) < 1e-8] = 1e-8

  # assign the results to a list called obj_alpha
  # obj_alpha[[i]] saves parameters associated with latent class (i+1)

  obj_alpha = list()
  for (i in 1:(num_class-1)){
    obj_alpha[[i]] =  obj_alpha_vec[seq(from=i,to=(num_class-1)*num_covx+i,by=num_class-1)]
  }

  ###### Obtain posterior membership probability tau0 ######


  # p is the fitted membership probability by the multinomial model

  p <- matrix(0,ncol=num_class,nrow=n)
  p[,1] = 1
  for (i in 2:num_class){
    p[,i] = exp(cbind(1,as.matrix(baseline[,covx_lb])) %*% obj_alpha[[i-1]])
  }
  psum = rowSums(p)
  p = apply(p,2,function(x) x/psum)

  # tau0 is the posterior membership probability

  pew <- p*ew
  pewsum = rowSums(pew)
  tau0 = apply(pew,2,function(x) x/pewsum)
  tau0[tau0<1e-8] = 1e-8

  ##### compute variance estimation

  if (varest == T){
    cat('Variance estimation started \n')
    x <- as.matrix(baseline[,covx_lb])
    ASE <- VarEst(obj_beta0,obj_beta1,obj_phi,obj_gamma,tau0,p,dat,x,y,Y_dist,balanced)
    rownames(ASE$alpha) = c('intercept',covx)
    rownames(ASE$beta0) = features
    rownames(ASE$beta1) = features
  }

  ######################################################################

  obj_alpha = matrix(unlist(obj_alpha),ncol=num_class-1,nrow=(1+ncol(as.matrix(baseline[,covx_lb]))))
  rownames(obj_alpha) = c('intercept',covx)
  rownames(obj_beta0) = features
  rownames(obj_beta1) = features
  rownames(obj_phi) = features
  rownames(obj_gamma) = features


  ##### output results as a list

  output = list()
  output$alpha = obj_alpha
  output$beta0 = obj_beta0
  output$beta1 = obj_beta1
  output$phi = obj_phi
  output$gamma = obj_gamma
  if(varest == T){
    output$ASE = ASE
  }
  output$tau = tau0
  output$qic = list(EQICA=eqica,EQICB=eqicb,CEEQIC=ceeqic)
  output$diff = diff

  # list(alpha = obj_alpha,  # parameters for latent class model
  #      beta0 = obj_beta0,  # intercepts for GEE model (features as rows, class as columns)
  #      beta1 = obj_beta1,  # slope for GEE model
  #      phi = obj_phi,      # scale parameter
  #      gamma = obj_gamma,  # AR1 structure parameter
  #      ASE = ASE,          # variance estimates
  #      tau = tau0,         # posterior membership prob
  #      qic = list(EQICA=eqica,EQICB=eqicb,CEEQIC=ceeqic),  # information criteria
  #      diff=diff#,          # convergence criteria
  #      #ew=ew             # exp(w)
  # )

  return(output)
}
