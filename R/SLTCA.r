#' SLTCA: Scalable and Robust Latent Trajectory Class Analysis Using Artificial Likelihood
#'
#' @description Conduct latent trajectory class analysis with longitudinal observations.
#' @param k Number of random initialization to start the algorithm.
#' @param num_class Number of latent classes in the fitted model.
#' @param dat Input data matrix.
#' @param id Column name in the data matrix `dat` for the patient id.
#' @param time Column name in the data matrix `dat` for the time of longitudinal observations.
#' @param num_obs Column name in the data matrix `dat` for the number of longitudinal observations (number of visits).
#' @param features A vector of column names in the data matrix `dat` for the longitudinal observations.
#' @param Y_dist A vector indicating the type of longitudinal observations. An element of Y_dist can be 'normal','bin', and 'poi' for continuous, binary and count data.
#' @param covx A vector of column names in the data matrix `dat` for baseline latent class risk factors.
#' @param ipw Column name in the data matrix `dat` for the inverse probability weights for missingness. ipw=1 if not specified.
#' @param stop Stopping criterion for the algorithm. stop can be either 'tau' based on posterior probabilities or 'par' based on point estimation.
#' @param tol A constant such that the algorithm stops if the stopping criterion is below this constant.
#' @param max Maximum number of iterations if the algorithm does not converge.
#' @param varest True or False: whether conduct variance estimation or not.
#' @param MSC Model selection criteria: 'AQIC','BQIC' or 'EQIC'.
#' @param balanced True or False: whether the longitudinal observations are equally spaced.
#' @param verbose Output progress of fitting the model.
#' @author Teng Fei. Email: <tfei@emory.edu>
#' @return A list with point estimates (alpha, beta0, beta1, phi, gamma), variance estimates (ASE), posterior membership probabilities (tau) and QICs (qic) of the latent trajectory class model. Point estimates and variance estimates are provided in matrix format, where columns represent latent classes and rows represent covariates or longitudinal features.
#' @references Hart, K.R., Fei, T. and Hanfelt, J.J. (2020), Scalable and robust latent trajectory class analysis using artificial likelihood. Biometrics. Accepted Author Manuscript <doi:10.1111/biom.13366>.
#' @examples
#'
#' dat <- simulation(500)
#'
#' res <- SLTCA(k=1,dat,num_class=2,"id","time","num_obs",paste("y.",1:6,sep=''),
#'              Y_dist=c('poi','poi','bin','bin','normal','normal'),
#'              "baselinecov",1,stop="tau",tol=0.005,max=50,
#'              varest=TRUE,balanced=TRUE,MSC='EQIC',verbose=TRUE)
#'
#' @importFrom stats as.formula binomial coef dist fitted gaussian poisson rbinom rmultinom rpois runif weights
#' @export

SLTCA <- function(k = 20,dat,num_class,id,time,num_obs,features,Y_dist,covx,ipw,stop,tol=0.005,max=50,varest=T,balanced=T,MSC='EQIC',verbose=T){

  #require(Matrix)
  #require(VGAM)
  #require(geepack)
  requireNamespace("Matrix")
  requireNamespace("VGAM")
  requireNamespace("geepack")

  IC = Inf

  if(MSC == 'AQIC'){
    for (i in 1:k){
      cat('random initialization',i,'\n')
      sol <- pointest(dat,num_class,id,time,num_obs,features,Y_dist,covx,ipw,stop,tol,max,varest,balanced,verbose)
      if (sol$qic[[1]] < IC){
        best_sol <- sol
        IC = sol$qic[[1]]
      }
    }
  }else if (MSC == 'BQIC'){
    for (i in 1:k){
      cat('random initialization',i,'\n')
      sol <- pointest(dat,num_class,id,time,num_obs,features,Y_dist,covx,ipw,stop,tol,max,varest,verbose)
      if (sol$qic[[2]] < IC){
        best_sol <- sol
        IC = sol$qic[[2]]
      }
    }
  }else if (MSC == 'EQIC'){
    for (i in 1:k){
      cat('random initialization',i,'\n')
      sol <- pointest(dat,num_class,id,time,num_obs,features,Y_dist,covx,ipw,stop,tol,max,varest,verbose)
      if (sol$qic[[3]] < IC){
        best_sol <- sol
        IC = sol$qic[[3]]
      }
    }
  }else{
    print('Error: MSC undefined.')
  }

  return(best_sol)
}
