LinProj <- function(beta0,beta1,phi,gamma,dat,y,Y_dist,balanced=F){

  # compute the approximated likelihood exp(w)
  #require(Matrix)
  requireNamespace("Matrix")
  num_class = ncol(phi)
  num_feature = nrow(phi)

  nmax <- max(dat$num_obs)
  n <- length(unique(dat$id))

  # LP initialized below saves approximated log likelihood w

  LP <- matrix(0,nrow=n,ncol=num_class)

  # I created two versions, one for equally spaced observation times, one for irregularly spaced observations.
  # The following shows the version used on equally spaced observations.

  if (balanced){
    mu <- list()
    v <- list()
    for (c in 1:num_class){
      mu[[c]] <- list()
      v[[c]] <- list()
      for (j in 1:num_feature){

        # define link functions

        if (Y_dist[j] == 'normal'){
          mulink <- function(x) x
          varlink <- function(x) rep(1,length(x))
        }else if (Y_dist[j] == 'poi'){
          mulink <- function(x) exp(x)
          varlink <- function(x) x
        }else if (Y_dist[j] == 'bin'){
          mulink <- function(x) exp(x)/(1+exp(x))
          varlink <- function(x) x*(1-x)
        }

        # find the mean vector and covariance matrix for jth marker in class c
        # from time 0 to the max time observed among all observations

        mu[[c]][[j]] <- mulink(t(beta0[j,c] + beta1[j,c]%*%t(dat$time[dat$id == dat$id[dat$num_obs == nmax][1]])))
        v[[c]][[j]] <- phi[j,c]*diag(as.vector(sqrt(varlink(mu[[c]][[j]])))) %*% (gamma[j,c]^as.matrix(dist(1:nmax))) %*% diag(as.vector(sqrt(varlink(mu[[c]][[j]]))))
      }

    }

    for (i in 1:n){

      # ii is the scalar observation id
      ii = unique(dat$id)[i]

      for (j in 1:ncol(LP)){
        if (j == 1){
          LP[i,j] = 0
        }else{

          # ilen is the label of observations for the jth feature for observation ii.
          ilen <- nrow(dat[dat$id==ii,]) - sum(is.na(y[dat$id==ii]))

          # When the observation times are equally spaced, can directly subset
          # already craeted mean vector and covariance matrix from time 0 to max

          muj <- as.vector(matrix(unlist(mu[[j]]),ncol=num_feature)[1:ilen,])
          mu1 <- as.vector(matrix(unlist(mu[[1]]),ncol=num_feature)[1:ilen,])
          vj <- lapply(v[[j]],function(x) x[1:ilen,1:ilen])
          v1 <- lapply(v[[1]],function(x) x[1:ilen,1:ilen])

          # compute w

          LP[i,j] = as.numeric(0.5*(unlist(muj)-unlist(mu1))%*%(solve(Matrix::bdiag(vj))%*%(as.vector(y[dat$id==ii,])-unlist(muj)) + solve(Matrix::bdiag(v1))%*%(as.vector(y[dat$id==ii,])-unlist(mu1))))
        }
      }
    }

  }else{

    # the following applies when the observations aer irregularly spaced

    mu <- list()
    v <- list()
    for (c in 1:num_class){
      mu[[c]] <- list()
      v[[c]] <- list()
      for (i in 1:n){
        mu[[c]][[i]] <- list()
        v[[c]][[i]] <- list()
        for (j in 1:num_feature){
          if (Y_dist[j] == 'normal'){
            mulink <- function(x) x
            varlink <- function(x) rep(1,length(x))
          }else if (Y_dist[j] == 'poi'){
            mulink <- function(x) exp(x)
            varlink <- function(x) x
          }else if (Y_dist[j] == 'bin'){
            mulink <- function(x) exp(x)/(1+exp(x))
            varlink <- function(x) x*(1-x)
          }
          mu[[c]][[i]][[j]] <- mulink(t(beta0[j,c] + beta1[j,c]%*%t(dat$time[dat$id==i])))

          if(length(mu[[c]][[i]][[j]])==1){
            v[[c]][[i]][[j]] <- phi[j,c]*varlink(mu[[c]][[i]][[j]])
          }else{
            v[[c]][[i]][[j]] <- phi[j,c]*diag(as.vector(sqrt(varlink(mu[[c]][[i]][[j]])))) %*% (gamma[j,c]^as.matrix(dist(dat$num_obs[dat$id==i]))) %*% diag(as.vector(sqrt(varlink(mu[[c]][[i]][[j]]))))
          }
        }
      }
    }

    for (i in 1:n){
      #ii = unique(dat$id)[i]
      for (c in 1:ncol(LP)){
        if (c == 1){
          LP[i,c] = 0
        }else{
          for (j in 1:num_feature){
            yj <- y[dat$id==i,j]
            nalab <- !is.na(yj)
            if(sum(nalab) > 0){
              if(length(mu[[c]][[i]][[j]])==1){
                muc = mu[[c]][[i]][[j]][nalab]
                mu1 = mu[[1]][[i]][[j]][nalab]
                vc = v[[c]][[i]][[j]][nalab]
                v1 = v[[1]][[i]][[j]][nalab]
              }else{
                muc = mu[[c]][[i]][[j]][nalab]
                mu1 = mu[[1]][[i]][[j]][nalab]
                vc = v[[c]][[i]][[j]][nalab,nalab]
                v1 = v[[1]][[i]][[j]][nalab,nalab]
              }
              yj = yj[nalab]
              LP[i,c] = LP[i,c] + as.numeric(0.5*(muc-mu1)%*%(solve(vc)%*%(yj-muc) + solve(v1)%*%(yj-mu1)))
            }
          }
        }
      }
    }
  }

  # output exponentiated w, exp(w)

  exp(LP)
}
