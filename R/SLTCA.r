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
#' @param verbose Output progress of fitting the model.
#' @author Teng Fei. Email: tfei@emory.edu

SLTCA <- function(k = 20,dat,num_class,id,time,num_obs,features,Y_dist,covx,ipw,stop,tol=0.005,max=50,varest=T,MSC='EQIC',verbose=T){
  IC = Inf

  if(MSC == 'AQIC'){
    for (i in 1:k){
      sol <- pointest(dat,num_class,id,time,num_obs,features,Y_dist,covx,ipw,stop,tol,max,varest,verbose)
      if (sol$qic[[1]] < IC){
        cat(i)
        best_sol <- sol
        IC = sol$qic[[1]]
      }
    }
  }else if (MSC == 'BQIC'){
    for (i in 1:k){
      sol <- pointest(dat,num_class,id,time,num_obs,features,Y_dist,covx,ipw,stop,tol,max,varest,verbose)
      if (sol$qic[[2]] < IC){
        cat(i)
        best_sol <- sol
        IC = sol$qic[[2]]
      }
    }
  }else if (MSC == 'EQIC'){
    for (i in 1:k){
      sol <- pointest(dat,num_class,id,time,num_obs,features,Y_dist,covx,ipw,stop,tol,max,varest,verbose)
      if (sol$qic[[3]] < IC){
        cat(i)
        best_sol <- sol
        IC = sol$qic[[3]]
      }
    }
  }else{
    print('Error: MSC undefined.')
    break
  }

  return(best_sol)
}
