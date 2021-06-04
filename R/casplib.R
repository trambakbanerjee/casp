require(MASS)
require(far)

#' rmt.est
#'
#' Provides (1) efficient and bias corrected estimates of the leading eigenvalues
#' of Sigma under a spiked covariance model, and (2) asymptotic adjustment factors
#' for the sample eigenvectors.
#'
#' @param K the number of spikes
#' @param S the sample covariance matrix
#' @param mw sample size
#
#' @return
#' \enumerate{
#' \item l0.hat - bias corrected estimate of the unknown noise level
#' \item l.hat - bias corrected estimates of the leading K eignevalues of Sigma
#' \item pj - sample eigenvectors of \eqn{\mathbf{S}}
#' \item zeta - asymptotic adjustment factors for the eigenvectors of S
#' \item K - the number of spikes (provided as an input)
#' }
#'
#' @details This function is called by \code{\link{casp.checkloss}} and \code{\link{casp.linexloss}}
#'
#' @seealso \code{\link{eigen.est}}, \code{\link{casp.checkloss}}, \code{\link{casp.linexloss}}
#'
#' @examples
#' library(casp)
#' K = 4
#' S = diag(c(10,8,6,4,rep(1,6)))
#' mw = 50
#' rmt.out<- rmt.est(K,S,mw)
#'
#' @references
#' \enumerate{
#' \item Debashis Paul. Asymptotics of sample eigenstructure for a large dimensional spiked covariance
#' model. Statistica Sinica, pages 1617-1642, 2007.
#' \item Alexei Onatski. Asymptotics of the principal components estimator of large factor models
#' with weakly influential factors. Journal of Econometrics, 168(2):244-258, 2012.
#' \item Damien Passemier, Zhaoyuan Li, and Jianfeng Yao. On estimation of the noise variance
#' in high dimensional probabilistic principal component analysis. Journal of the Royal
#' Statistical Society: Series B (Statistical Methodology), 79(1):51-67, 2017.
#' }
#'
#' @export

rmt.est<- function(K,S,mw){

  n<- dim(S)[1]
  spd<- eigen(S,symmetric = TRUE)
  rho<- n/(mw-1)

  l.tilde<-spd$values
  l0.tilde<- ((n-K)^{-1})*sum(l.tilde[(K+1):n])

  out<- eigen.est(l.tilde,l0.tilde,rho,K,n)
  l0.hat<-out$l0.hat
  l.hat<-out$l.hat

  zeta<-matrix(0,K,1)
  for(k in 1:K){

    temp<-l.hat[k]/l0.hat
    a<- rho/(temp-1)^2
    b<- rho/(temp-1)
    zeta[k]<- sqrt((1-a)/(1+b))
  }

  return(list("l0.hat"=l0.hat,"l.hat"=l.hat,"zeta"=zeta,
              "pj"=spd$vectors,"K"=K))
}

#' eigen.est
#'
#' Provides efficient and bias corrected estimates of the leading eigenvalues
#' of \eqn{\mathbf{\Sigma}} under a spiked covariance model.
#'
#' @param l eigenvalues of the sample covariance matrix
#' @param l.0 an estimate of noise via maximum likelihood
#' @param K the number of spikes
#' @param n the dimension of the covariance matrix
#' @param rho the ratio \eqn{n/(mw-1)} where \eqn{mw} is the sample size
#
#' @return
#' \enumerate{
#' \item l0.hat - bias corrected estimate of the unknown noise level
#' \item l.hat - bias corrected estimates of the leading K eignevalues of \eqn{\mathbf{\Sigma}}
#' }
#'
#' @details This function is called by \code{\link{rmt.est}} and the estimate of noise, l.0, that
#' is used here is \eqn{\ell_0=(n-K)^{-1}\sum_{j=K+1}^{n}\ell_j}.
#'
#' @seealso \code{\link{rmt.est}}, \code{\link{casp.checkloss}}, \code{\link{casp.linexloss}}
#'
#' @examples
#' library(casp)
#' K = 4
#' n = 10
#' l = c(10,8,6,4,rep(1,6))
#' l.0 = ((n-K)^{-1})*sum(l[(K+1):n])
#' mw = 100
#' rho = n/(mw-1)
#' eigen.est(l,l.0,rho,K,n)
#'
#' @references
#' \enumerate{
#' \item Debashis Paul. Asymptotics of sample eigenstructure for a large dimensional spiked covariance
#' model. Statistica Sinica, pages 1617-1642, 2007.
#' \item Alexei Onatski. Asymptotics of the principal components estimator of large factor models
#' with weakly influential factors. Journal of Econometrics, 168(2):244-258, 2012.
#' \item Damien Passemier, Zhaoyuan Li, and Jianfeng Yao. On estimation of the noise variance
#' in high dimensional probabilistic principal component analysis. Journal of the Royal
#' Statistical Society: Series B (Statistical Methodology), 79(1):51-67, 2017.
#' }
#'
#' @export

eigen.est<-function(l,l.0,rho,K,n){

  x<- matrix(0,K,1)
  b<-l.0

  for(k in 1:K){
    a<- l[k]
    temp.1<-b+a-rho*b
    x[k]<- (temp.1+sqrt((temp.1^2)-4*a*b))/(2)
  }

  temp.2<- ((x/l.0)-1)^{-1}
  xi<-K+sum(temp.2[1:K])

  l0.hat<- l.0*(1+rho*xi/(n-K))
  l.hat<-matrix(0,K,1)
  for(k in 1:K){

    a<- l[k]/l0.hat
    temp.3<- a+1-rho
    if (temp.3^2-4*a<0){
      l.hat[k]<-(l0.hat/2)*(temp.3)
      print('imaginary roots')
    }
    if (temp.3^2-4*a>=0){

      l.hat[k]<- (l0.hat/2)*(temp.3+sqrt(temp.3^2-4*a))
    }
  }
  return(list("l0.hat"=l0.hat,"l.hat"=l.hat))
}

#' hfun.est
#'
#' Provides an estimate of \eqn{\mathbf{H}_{r,\alpha,\beta}} using equation (13) of
#' the casp paper (see the reference).
#'
#' @param r a value in \eqn{\{-1,0,1\}}
#' @param alpha a non-negative number
#' @param tau the positive scale hyper-parameter for the prior on the locations
#' @param beta the non-negative shape hyper-parameter for the prior on the locations
#' @param rmt the output from \code{\link{rmt.est}}
#' @param mx the sample size of the past observations \eqn{\mathbf{X}}
#'
#' @return
#' \enumerate{
#' \item h - an estimate of \eqn{H_{r,\alpha,\beta}}
#' }
#'
#' @details This function relies on the output from \code{\link{rmt.est}} and is called by \code{\link{casp.checkloss}},
#' and \code{\link{casp.linexloss}}. Please see equation (13) of the casp paper in the reference
#' to get more details about the estimation.
#'
#' @seealso \code{\link{rmt.est}}, \code{\link{casp.checkloss}}, \code{\link{casp.linexloss}}
#'
#' @examples
#' library(casp)
#' r = -1
#' alpha = 0
#' tau = 1
#' beta = 0.5
#' mx = 1
#' K = 4
#' S = diag(c(10,8,6,4,rep(1,6)))
#' mw = 50
#' rmt<- rmt.est(K,S,mw)
#' H.hat<- hfun.est(r,alpha,beta,rmt,mx,tau)
#'
#' @references
#' \enumerate{
#' \item Trambak Banerjee, Gourab Mukherjee, and Debashis Paul. Improved Shrinkage Prediction under a Spiked
#' Covariance Structure, 2021.
#' }
#'
#' @export

hfun.est<-function(r,alpha,beta,rmt,mx,tau){

  pj<- rmt$pj
  l0.hat<- rmt$l0.hat
  l.hat<-rmt$l.hat
  zeta<-rmt$zeta
  K<-rmt$K

  n<- dim(pj)[1]
  h<-matrix(0,n,n)

  a.1<- (l0.hat^(r+alpha))/((1+(1/(mx*tau))*l0.hat^{1-beta})^r)

  for(k in 1:K){

    b.1<- (l.hat[k]^(r+alpha))/((1+(1/(mx*tau))*l.hat[k]^{1-beta})^r)
    h<- (1/zeta[k]^2)*(b.1-a.1)*(pj[,k]%*%t(pj[,k]))+h

  }
  h<-h+a.1*diag(n)
  return(h)
}

#' casp.checkloss
#'
#' Main function for shrinkage prediction under check loss.
#'
#' @importFrom stats qnorm
#'
#' @param X a \eqn{mx\times n} matrix of past observations.
#' @param S \eqn{n\times n} sample covariance matrix based on mw samples
#' @param b.tilde this is a \eqn{n\times 1} vector of check loss parameters. Here
#' \eqn{\tilde{b}_i=b_i/(b_i+h+i)} where \eqn{\mathbf{b}} and \eqn{\mathbf{h}}
#' are the check loss parameters for the \eqn{n} coordinates.
#' @param tau the positive scale hyper-parameter for the prior on the locations
#' @param beta the non-negative shape hyper-parameter for the prior on the locations
#' @param eta the prior mean of the locations
#' @param mx the sample size of past observations \eqn{\mathbf{X}}
#' @param mw the size of the side information for calculating the sample covariance
#' matrix \eqn{\mathbf{S}}
#' @param m0 the sample size of the future observation. Usually this is set to 1
#' @param type if type = 1 then all the shrinkage factors are equal to 1. If type = 0, then
#' the shrinkage factors are estimated.
#'
#' @return
#' \enumerate{
#' \item q - shrinkage prediction under check loss for the \eqn{n} coordinates
#' \item f - estimated shrinkage factors (equal to 1 if type = 1)
#' }
#'
#' @details This function is based on Definition 2 of the casp paper, and relies on
#' \code{\link{rmt.est}} and \code{\link{hfun.est}}. The shrinkage factors
#' are estimated using the formulation given in Definition 3 and use
#' \code{\link{f.est}} in the background. Please see the casp paper in the reference
#' for more details about these estimation techniques. If \eqn{(\tau,\beta)} are unknown
#' then one may first use \code{\link{taubeta.casp.est}} to estimate them and then
#' use the estimated values in \code{\link{casp.checkloss}}.
#'
#' @seealso \code{\link{casp.linexloss}},\code{\link{f.est}},\code{\link{rmt.est}},\code{\link{taubeta.casp.est}}
#'
#' @examples
#' library(casp)
#' set.seed(42)
#' n = 10
#' mx = 5
#' S = diag(c(10,8,6,4,rep(1,n-4)))
#' X<- matrix(runif(mx*n),mx,n)
#' tau = 1
#' beta = 0.5
#' eta = rep(0,n)
#' mw = 100
#' m0 = 1
#' b.tilde = rep(0.5,n)
#' q.casp<- casp.checkloss(X,S,b.tilde,tau,beta,eta,mx,mw,m0,0)
#'
#' @references
#' \enumerate{
#' \item Trambak Banerjee, Gourab Mukherjee, and Debashis Paul. Improved Shrinkage Prediction under a Spiked
#' Covariance Structure, 2021.
#' }
#'
#' @export

casp.checkloss<- function(X,S,b.tilde,tau,beta,eta,mx,mw,m0,type){

  n<- dim(S)[1]
  K.est<- KN.test(S, mw, 0.05)
  rmt<- rmt.est(K.est$numOfSpikes,S,mw)

  f<- rep(1,n)
  if (type == 0){

    f<- f.est(rmt,beta,tau,mx)

  }

  h1<- hfun.est(1,-1,beta,rmt,mx,tau)
  h2<- hfun.est(1,0,beta,rmt,mx,tau)
  h3<- hfun.est(0,1,0,rmt,mx,tau)
  term1<- eta
  if(mx>1){
    term2<- f*(h1%*%(colMeans(X)-eta))
  }
  if(mx==1){
    term2<- f*(h1%*%(X-eta))
  }
  term3<- sqrt(mx^{-1}*diag(h2)+m0^{-1}*diag(h3))

  q<- term1+term2+term3*qnorm(b.tilde)

  return(list("q"=q,"f"=f,"h1"=h1,"h2"=h2,"h3"=h3))

}

#' casp.linexloss
#'
#' Main function for shrinkage prediction under Linex loss.
#'
#' @param X a \eqn{mx\times n} matrix of past observations.
#' @param S \eqn{n\times n} sample covariance matrix based on mw samples.
#' @param a this is a \eqn{n\times 1} vector of Linex loss parameter. See equation (4)
#' in the casp paper.
#' @param tau the positive scale hyper-parameter for the prior on the locations.
#' @param beta the non-negative shape hyper-parameter for the prior on the locations.
#' @param eta the prior mean of the locations.
#' @param mx the sample size of past observations \eqn{\mathbf{X}}.
#' @param mw the size of the side information for calculating the sample covariance
#' matrix \eqn{\mathbf{S}}.
#' @param m0 the sample size of the future observation. Usually this is set to 1.
#' @param type if type = 1 then all the shrinkage factors are equal to 1. If type = 0, then
#' the shrinkage factors are estimated.
#'
#' @return
#' \enumerate{
#' \item q - shrinkage prediction under Linex loss for the \eqn{n} coordinates
#' \item f - estimated shrinkage factors (equal to 1 if type = 1)
#' }
#'
#' @details This function is based on Definition 2 of the casp paper, and relies on
#' \code{\link{rmt.est}} and \code{\link{hfun.est}}. The shrinkage factors
#' are estimated using the formulation given in Definition 3 and use
#' \code{\link{f.est}} in the background. Please see the casp paper in the reference
#' for more details about these estimation techniques. If \eqn{(\tau,\beta)} are unknown
#' then one may first use \code{\link{taubeta.casp.est}} to estimate them and then
#' use the estimated values in \code{\link{casp.linexloss}}.
#'
#' @seealso \code{\link{casp.checkloss}},\code{\link{f.est}},\code{\link{rmt.est}},\code{\link{taubeta.casp.est}}
#'
#' @examples
#' library(casp)
#' set.seed(42)
#' n = 10
#' mx = 5
#' S = diag(c(10,8,6,4,rep(1,n-4)))
#' X<- matrix(runif(mx*n),mx,n)
#' tau = 1
#' beta = 0.5
#' eta = rep(0,n)
#' mw = 100
#' m0 = 1
#' a = rep(-1,n)
#' q.casp<- casp.linexloss(X,S,a,tau,beta,eta,mx,mw,m0,0)
#'
#' @references
#' \enumerate{
#' \item Trambak Banerjee, Gourab Mukherjee, and Debashis Paul. Improved Shrinkage Prediction under a Spiked
#' Covariance Structure, 2021.
#' }
#'
#' @export

casp.linexloss<- function(X,S,a,tau,beta,eta,mx,mw,m0,type){

  n<- dim(S)[1]
  K.est<- KN.test(S, mw, 0.05)
  rmt<- rmt.est(K.est$numOfSpikes,S,mw)

  f<- rep(1,n)
  if (type == 0){

    f<- f.est(rmt,beta,tau,mx)

  }

  h1<- hfun.est(1,-1,beta,rmt,mx,tau)
  h2<- hfun.est(1,0,beta,rmt,mx,tau)
  h3<- hfun.est(0,1,0,rmt,mx,tau)
  term1<- eta
  if(mx>1){
    term2<- f*(h1%*%(colMeans(X)-eta))
  }
  if(mx==1){
    term2<- f*(h1%*%(X-eta))
  }
  term3<- mx^{-1}*diag(h2)+m0^{-1}*diag(h3)

  q<- term1+term2-term3*(a/2)
  return(list("q"=q,"f"=f,"h1"=h1,"h2"=h2,"h3"=h3))

}

#' f.est
#'
#' Estimates the coordinate-wise shrinkage factors using the information in
#' the sample covariance matrix \eqn{\mathbf{S}}. See Definition (3) of
#' the casp paper (see the reference).
#'
#' @param rmt the output from \code{\link{rmt.est}}
#' @param tau the positive scale hyper-parameter for the prior on the locations
#' @param beta the non-negative shape hyper-parameter for the prior on the locations
#' @param mx the sample size of the past observations \eqn{\mathbf{X}}
#'
#' @return
#' \enumerate{
#' \item f - estimated shrinkage factors
#' }
#'
#' @details This function relies on the output from \code{\link{rmt.est}} and is called by \code{\link{casp.checkloss}},
#' and \code{\link{casp.linexloss}}. Please see Definition (3) of the casp paper in the reference
#' for more details.
#'
#' @seealso \code{\link{rmt.est}}, \code{\link{casp.checkloss}}, \code{\link{casp.linexloss}}
#'
#' @examples
#' library(casp)
#' tau = 1
#' beta = 0.5
#' mx = 1
#' K = 4
#' S = diag(c(10,8,6,4,rep(1,6)))
#' mw = 50
#' rmt<- rmt.est(K,S,mw)
#' f<- f.est(rmt,beta,tau,mx)
#'
#' @references
#' \enumerate{
#' \item Trambak Banerjee, Gourab Mukherjee, and Debashis Paul. Improved Shrinkage Prediction under a Spiked
#' Covariance Structure, 2021.
#' }
#'
#' @export

f.est<-function(rmt,beta,tau,mx){

  pj<- rmt$pj
  l0.hat<- rmt$l0.hat
  l.hat<-rmt$l.hat
  zeta<-rmt$zeta
  K<- rmt$K

  n<- dim(pj)[1]
  H<-matrix(0,n,n)

  a.1<- ((mx*tau)*l0.hat^beta)/(1+(1/(mx*tau))*l0.hat^{1-beta})

  for(k in 1:K){

    b.1<- (mx*tau*l.hat[k]^beta)/(1+(1/(mx*tau))*l.hat[k]^{1-beta})
    H<- (1/zeta[k]^2)*(b.1-a.1)*(pj[,k]%*%t(pj[,k]))+H

  }
  H<-H+a.1*diag(n)
  numerator<- diag(H)

  g<- l0.hat+(mx*tau)*l0.hat^beta
  M<-matrix(0,n,n)

  a.1<-1/(1+(1/(mx*tau))*l0.hat^{1-beta})

  for(k in 1:K){

    b.1<- 1/(1+(1/(mx*tau))*l.hat[k]^{1-beta})
    M<- (1/zeta[k]^4)*((b.1-a.1)^2)*(pj[,k]%*%t(pj[,k]))+M
  }

  R<- H+g*M
  denominator<- diag(R)
  return(numerator/denominator)

}

#' g.est
#'
#' Provides an estimate of (1) \eqn{\tau H_{-1,1+\beta,\beta}}, and
#' (2) \eqn{\tau^{-1} H_{1,-1-\beta,\beta}} using equation (13) of
#' the casp paper (see the reference).
#'
#' @param rmt the output from \code{\link{rmt.est}}.
#' @param tau the positive scale hyper-parameter for the prior on the locations.
#' @param beta the non-negative shape hyper-parameter for the prior on the locations.
#' @param mx the sample size of the past observations \eqn{\mathbf{X}}.
#'
#' @return
#' \enumerate{
#' \item G - an estimate of \eqn{\tau H_{-1,1+\beta,\beta}}
#' \item Ginv - an estimate of \eqn{\tau^{-1} H_{1,-1-\beta,\beta}}
#' }
#'
#' @details This function relies on the output from \code{\link{rmt.est}} and is called by \code{\link{casp.checkloss}},
#' and \code{\link{taubeta.casp.est}}. Please see equation (13) of the casp paper in the reference
#' to get more details about the estimation.
#'
#' @seealso \code{\link{rmt.est}}, \code{\link{taubeta.casp.est}}
#'
#' @examples
#' library(casp)
#' tau = 1
#' beta = 0.5
#' mx = 1
#' K = 4
#' S = diag(c(10,8,6,4,rep(1,6)))
#' mw = 50
#' rmt<- rmt.est(K,S,mw)
#' g.out<- g.est(rmt,beta,tau,mx)
#'
#' @references
#' \enumerate{
#' \item Trambak Banerjee, Gourab Mukherjee, and Debashis Paul. Improved Shrinkage Prediction under a Spiked
#' Covariance Structure, 2021.
#' }
#'
#' @export

g.est<-function(rmt,beta,tau,mx){

  pj<- rmt$pj
  l0.hat<- rmt$l0.hat
  l.hat<-rmt$l.hat
  zeta<-rmt$zeta
  K<- rmt$K

  n<- dim(pj)[1]
  G<-matrix(0,n,n)
  a.1<- l0.hat+(mx*tau)*l0.hat^beta
  for(k in 1:K){

    b.1<- l.hat[k]+(mx*tau*l.hat[k]^beta)
    G<- (1/zeta[k]^2)*(b.1-a.1)*(pj[,k]%*%t(pj[,k]))+G

  }
  G<-G+a.1*diag(n)

  Ginv<-matrix(0,n,n)
  a.2<- 1/(l0.hat+(mx*tau)*l0.hat^beta)
  for(k in 1:K){

    b.2<- 1/(l.hat[k]+(mx*tau*l.hat[k]^beta))
    Ginv<- (1/zeta[k]^2)*(b.2-a.2)*(pj[,k]%*%t(pj[,k]))+Ginv

  }
  Ginv<-Ginv+a.2*diag(n)
  return(list("G"=G,"Ginv"=Ginv))

}

#' Data-driven estimation of the hyper-parameters \eqn{(\tau,\beta)}
#'
#' Provides an estimate of the hyper-parameters \eqn{(\tau,\beta)} using equation (15) of
#' the casp paper (see the reference).
#'
#' @param grid.val a matrix with two columns. Each row of the matrix represents an
#' element of the two dimensional grid. The first component represents a likely value for \eqn{\tau}
#' while the second component is a likely value for \eqn{\beta}.
#' @param X a \eqn{mx\times n} matrix of past observations.
#' @param rmt the output from \code{\link{rmt.est}}.
#' @param mx the sample size of the past observations \eqn{\mathbf{X}}.
#'
#' @return
#' \enumerate{
#' \item An estimate of \eqn{(\tau,\beta)}
#' }
#'
#' @details This function relies on the output from \code{\link{rmt.est}} and calls \code{\link{g.est}}.
#' Please see Section 3.3 of the casp paper for more details.
#'
#' @seealso \code{\link{rmt.est}}, \code{\link{casp.checkloss}}, \code{\link{casp.linexloss}}
#'
#' @examples
#' library(casp)
#' set.seed(42)
#' n = 10
#' mx = 5
#' X = matrix(runif(mx*n),mx,n)
#' K = 4
#' S = diag(c(10,8,6,4,rep(1,n-4)))
#' mw = 50
#' rmt<- rmt.est(K,S,mw)
#' tau.grid<-c(0.2,0.3,0.4,0.5)
#' beta.grid<-c(0.15,0.25,0.35,0.5)
#' grid.val<- cbind(rep(tau.grid,each = length(beta.grid)),
#'                     rep(beta.grid,length(beta.grid)))
#' taubeta.estimated<- taubeta.casp.est(grid.val,rmt,mx,X)
#'
#' @references
#' \enumerate{
#' \item Trambak Banerjee, Gourab Mukherjee, and Debashis Paul. Improved Shrinkage Prediction under a Spiked
#' Covariance Structure, 2021.
#' }
#'
#' @export

taubeta.casp.est<-function(grid.val,rmt,mx,X){

  cv.taubeta<-matrix(0,nrow(grid.val),1)
  if (mx>1){

    X<-colMeans(X)
  }
  for (t in 1:nrow(grid.val)){

    ttau<- grid.val[t,1]
    bbeta<-grid.val[t,2]

    S.q<- g.est(rmt,bbeta,ttau,mx)$G
    Sinv.q<- g.est(rmt,bbeta,ttau,mx)$Ginv
    cv.taubeta[t]<- -0.5*log(det(S.q))-0.5*t(X)%*%Sinv.q%*%X

  }
  o<- order(-cv.taubeta)
  return(colMeans(grid.val[o[1:3],]))
}



