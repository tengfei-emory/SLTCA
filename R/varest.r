VarEst <- function(beta0,beta1,phi,gamma,tau1,p,dat,x,y,Y_dist,balanced=T){
  #require(Matrix)
  requireNamespace("Matrix")
  num_class = ncol(phi)
  num_feature = nrow(phi)

  nmax <- max(dat$num_obs)
  n <- length(unique(dat$id))

  bigA <- list()

  if (balanced){
    mu <- list()
    v <- list()
    dmu <- list()
    for (c in 1:num_class){
      mu[[c]] <- list()
      v[[c]] <- list()
      dmu[[c]] <-list()
      bigA[[c]] <- 0
      for (j in 1:num_feature){
        if (Y_dist[j] == 'normal'){
          mulink <- function(x) x
          varlink <- function(x) rep(1,length(x))
          dmulink <- function(x) 1
        }else if (Y_dist[j] == 'poi'){
          mulink <- function(x) exp(x)
          varlink <- function(lambda) lambda
          dmulink <- function(x) exp(x)
        }else if (Y_dist[j] == 'bin'){
          mulink <- function(x) exp(x)/(1+exp(x))
          varlink <- function(p) p*(1-p)
          dmulink <- function(x) exp(x)/(1+exp(x))^2 #(exp(x)*(1+exp(x)) - exp(x)*exp(x))/(1+exp(x))^2
        }
        mu[[c]][[j]] <- mulink(t(beta0[j,c] + beta1[j,c]%*%t(dat$time[dat$id == dat$id[dat$num_obs == nmax][1]])))
        v[[c]][[j]] <- phi[j,c]*diag(as.vector(sqrt(varlink(mu[[c]][[j]])))) %*% (gamma[j,c]^as.matrix(dist(1:nmax))) %*% diag(as.vector(sqrt(varlink(mu[[c]][[j]]))))
        dmu[[c]][[j]] <- cbind(dmulink(t(beta0[j,c] + beta1[j,c]%*%t(dat$time[dat$id == dat$id[dat$num_obs == nmax][1]]))),
                               dmulink(t(beta0[j,c] + beta1[j,c]%*%t(dat$time[dat$id == dat$id[dat$num_obs == nmax][1]])))*dat$time[dat$id == dat$id[dat$num_obs == nmax][1]])

      }

    }

    #print(dmu)

    bigB <- 0
    bigBB <- 0
    bigBalpha <- 0
    bigAalpha <- 0

    for(i in 1:n){
      ii = unique(dat$id)[i]
      ilen <- nrow(dat[dat$id==ii,])
      bigBi <- list()
      bigBa <- 0
      xi <- c(1,x[i,])
      #qalpha =qalpha - kronecker(p[i,],xi)
      for (c in 1:num_class){

        qalpha <- rep(0,(num_class-1)*(ncol(x)+1))
        muc <- as.vector(matrix(unlist(mu[[c]]),ncol=num_feature)[1:ilen,])
        vc <- lapply(v[[c]],function(x) x[1:ilen,1:ilen])
        dmuc <- lapply(dmu[[c]],function(x) x[1:ilen,])
        q <- as.vector(t(as.matrix(Matrix::bdiag(dmuc))) %*% solve(as.matrix(Matrix::bdiag(vc))) %*% (as.vector(y[dat$id==ii,])-unlist(muc)))
        bigBi[[c]] <- (tau1[i,c] * q)#%o%(tau1[i,c] * q)
        bigA[[c]] <- bigA[[c]] + p[i,c]*t(as.matrix(Matrix::bdiag(dmuc))) %*% as.matrix(solve(Matrix::bdiag(vc))) %*% as.matrix(Matrix::bdiag(dmuc)) - tau1[i,c]*q%o%q

        Dalpha <- -outer(p[i,-1],p[i,-1])
        diag(Dalpha) <- p[i,-1]*(1-p[i,-1])
        for (j in 2:num_class){
          qalpha[((j-1)*(ncol(x)+1)-ncol(x)):((j-1)*(ncol(x)+1))] = qalpha[((j-1)*(ncol(x)+1)-ncol(x)):((j-1)*(ncol(x)+1))] - p[i,j]*c(1,x[i,])
        }
        if (c > 1){
          qalpha[((c-1)*(ncol(x)+1)-ncol(x)):((c-1)*(ncol(x)+1))] = qalpha[((c-1)*(ncol(x)+1)-ncol(x)):((c-1)*(ncol(x)+1))] + c(1,x[i,])
        }
        bigBa = bigBa + tau1[i,c]*qalpha
        bigAalpha = bigAalpha + p[i,c]*kronecker(Dalpha,outer(xi,xi)) - tau1[i,c]*qalpha%o%qalpha

      }
      bigBalpha = bigBalpha + outer(bigBa,bigBa)
      bigBB = bigBB + outer(unlist(bigBi),unlist(bigBi))
      bigB = bigB + outer(c(bigBa,unlist(bigBi)),c(bigBa,unlist(bigBi)))
    }

    bigAAA <- Matrix::bdiag(c(list(bigAalpha),bigA)) + bigB
    Sigma <- solve(bigAAA) %*% bigB %*% solve(bigAAA)
    ASE0 <- sqrt(diag(Sigma))

    # bigAA <- bdiag(bigA) + bigBB
    # Sigma <- solve(bigAA)%*%bigBB%*%solve(bigAA)
    ASE <- ASE0[((ncol(x)+1)*(num_class-1)+1):length(ASE0)]

    # bigAalpha = bigAalpha + bigBalpha
    # Sigmaalpha <- solve(bigAalpha)%*%bigBalpha%*%solve(bigAalpha)
    ASEalpha <- ASE0[1:((ncol(x)+1)*(num_class-1))]

  }else{

    mu <- list()
    v <- list()
    dmu <- list()

    bigB <- 0
    bigBB <- 0
    bigBalpha <- 0
    bigAalpha <- 0

    for(i in 1:n){
      ii = unique(dat$id)[i]
      #ilen <- nrow(dat[dat$id==ii,])
      bigBi <- list()
      bigBa <- 0
      xi <- c(1,x[i,])
      #qalpha =qalpha - kronecker(p[i,],xi)
      for (c in 1:num_class){
        mu[[c]] <- list()
        v[[c]] <- list()
        dmu[[c]] <-list()
        bigA[[c]] <- 0

        qalpha <- rep(0,(num_class-1)*(ncol(x)+1))

        adj = 0
        nalabel2 <- NULL

        for (j in 1:num_feature){
          if (Y_dist[j] == 'normal'){
            mulink <- function(x) x
            varlink <- function(x) rep(1,length(x))
            dmulink <- function(x) 1
          }else if (Y_dist[j] == 'poi'){
            mulink <- function(x) exp(x)
            varlink <- function(lambda) lambda
            dmulink <- function(x) exp(x)
          }else if (Y_dist[j] == 'bin'){
            mulink <- function(x) exp(x)/(1+exp(x))
            varlink <- function(p) p*(1-p)
            dmulink <- function(x) exp(x)/(1+exp(x))^2 #(exp(x)*(1+exp(x)) - exp(x)*exp(x))/(1+exp(x))^2
          }

          nalabel <- is.na(y[dat$id==ii,j])

          if (sum(!nalabel)>=1){
            mu[[c]][[j-adj]] <- mulink(t(beta0[j,c] + beta1[j,c]%*%t(unique(dat$time[dat$id == ii][!nalabel]))))
            if (sum(!nalabel)==1 | length(unique(dat$time[dat$id == ii][!nalabel])) == 1){
              v[[c]][[j-adj]] <- phi[j,c]*diag((sqrt(varlink(mu[[c]][[j-adj]])))) %*% (gamma[j,c]^as.matrix(dist(unique(dat$time[dat$id==ii][!nalabel])))) %*% diag((sqrt(varlink(mu[[c]][[j-adj]]))))

            }else{
              v[[c]][[j-adj]] <- phi[j,c]*diag(as.vector(sqrt(varlink(mu[[c]][[j-adj]])))) %*% (gamma[j,c]^as.matrix(dist(unique(dat$time[dat$id==ii][!nalabel])))) %*% diag(as.vector(sqrt(varlink(mu[[c]][[j-adj]]))))
            }

            dmu[[c]][[j-adj]] <- cbind(dmulink(t(beta0[j,c] + beta1[j,c]%*%t(unique(dat$time[dat$id==ii][!nalabel])))),
                                   dmulink(t(beta0[j,c] + beta1[j,c]%*%t(unique(dat$time[dat$id==ii][!nalabel]))))*unique(dat$time[dat$id==ii][!nalabel]))

          }else{
            adj = adj + 1
            nalabel2 <- c(nalabel2,j)
          }

        }

        nalabel <- is.na(y[dat$id==ii,])
        muc <- unlist(mu[[c]])
        vc <- v[[c]]
        dmuc <- dmu[[c]]

        if (length(unique(dat$time[dat$id == ii])) < length(dat$time[dat$id == ii])){
          q <- as.vector(t(as.matrix(Matrix::bdiag(dmuc))) %*% as.matrix(solve(Matrix::bdiag(vc))) %*% (as.vector(y[min(which(dat$id==ii)),][!nalabel[1,]])-muc))
        }else{
          q <- as.vector(t(as.matrix(Matrix::bdiag(dmuc))) %*% as.matrix(solve(Matrix::bdiag(vc))) %*% (as.vector(y[dat$id==ii,][!nalabel])-muc))
        }

        if (is.null(nalabel2)){
          bigBi[[c]] <- (tau1[i,c] * q)#%o%(tau1[i,c] * q)
        }else{
          bigBi[[c]] <- rep(0,2*num_feature)
          bigBi[[c]][-c(nalabel2,num_feature+nalabel2)] = (tau1[i,c] * q)
          bigBi[[c]] <- as.matrix(bigBi[[c]])
        }
        bigA[[c]] <- bigA[[c]] + p[i,c]*t(as.matrix(Matrix::bdiag(dmuc))) %*% solve(as.matrix(Matrix::bdiag(vc))) %*% as.matrix(Matrix::bdiag(dmuc)) - tau1[i,c]*q%o%q

        Dalpha <- -outer(p[i,-1],p[i,-1])
        diag(Dalpha) <- p[i,-1]*(1-p[i,-1])
        for (l in 2:num_class){
          qalpha[((l-1)*(ncol(x)+1)-(ncol(x)+1)+1):((l-1)*(ncol(x)+1))] = qalpha[((l-1)*(ncol(x)+1)-(ncol(x)+1)+1):((l-1)*(ncol(x)+1))] - p[i,l]*c(1,x[i,])
        }
        if (c > 1){
          qalpha[((c-1)*(ncol(x)+1)-(ncol(x)+1)+1):((c-1)*(ncol(x)+1))] = qalpha[((c-1)*(ncol(x)+1)-(ncol(x)+1)+1):((c-1)*(ncol(x)+1))] + c(1,x[i,])
        }
        bigBa = bigBa + tau1[i,c]*qalpha
        bigAalpha = bigAalpha + p[i,c]*kronecker(Dalpha,outer(xi,xi)) - tau1[i,c]*qalpha%o%qalpha
        #cat(sum(diag(bigAalpha)),' ',i,' ',c,'\n')

      }

      bigB = bigB + outer(c(bigBa,as.vector(unlist(bigBi))),c(bigBa,as.vector(unlist(bigBi))))
      bigBalpha = bigBalpha + outer(bigBa,bigBa)
      bigBB = bigBB + outer(as.vector(unlist(bigBi)),as.vector(unlist(bigBi)))
    }


    bigAAA <- Matrix::bdiag(c(list(bigAalpha),bigA)) + bigB
    Sigma <- solve(bigAAA) %*% bigB %*% solve(bigAAA)
    ASE0 <- sqrt(diag(Sigma))

    # bigAA <- bdiag(bigA) + bigBB
    # Sigma <- solve(bigAA)%*%bigBB%*%solve(bigAA)
    ASE <- ASE0[((ncol(x)+1)*(num_class-1)+1):length(ASE0)]

    # bigAalpha = bigAalpha + bigBalpha
    # Sigmaalpha <- solve(bigAalpha)%*%bigBalpha%*%solve(bigAalpha)
    ASEalpha <- ASE0[1:((ncol(x)+1)*(num_class-1))]

    # bigAA <- bdiag(bigA) + bigBB
    # Sigma <- solve(bigAA)%*%bigBB%*%solve(bigAA)
    # ASE <- sqrt(diag(Sigma))
    #
    # bigAalpha = bigAalpha + bigBalpha
    # Sigmaalpha <- solve(bigAalpha)%*%bigBalpha%*%solve(bigAalpha)
    # ASEalpha <- sqrt(diag(Sigmaalpha))

  }

  ASE_beta0 <- ASE[seq(1,num_feature*2*num_class,2)]
  ASE_beta1 <- ASE[seq(2,num_feature*2*num_class+1,2)]

  output = list()
  output$alpha = matrix(ASEalpha,ncol=num_class-1,nrow=(1+ncol(x)))
  output$beta0 = matrix(ASE_beta0,ncol=num_class,nrow=num_feature)
  output$beta1 = matrix(ASE_beta1,ncol=num_class,nrow=num_feature)

  return(output)
}
