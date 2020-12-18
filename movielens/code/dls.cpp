// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
colvec vec_rmvn(arma::colvec mu, arma::mat lchol_sig) {
  int ndim = lchol_sig.n_cols;
  colvec z = randn<colvec>(ndim); 
  colvec nvec = mu + lchol_sig * z;
  
  return nvec;
} 

// [[Rcpp::export]]
colvec runif_range(int n, double low, double up) {
  colvec u = randu<colvec>(n); 
  colvec uvec = low + (up - low) * u;
  
  return uvec;
} 

// [[Rcpp::export()]]
arma::mat mat_corr (double phi, arma::mat distMat) {
  arma::mat res = exp(-phi * distMat);

  return(res);
}

// [[Rcpp::export()]]
arma::mat euclidean_dist (arma::mat x, arma::mat y) {
  int xsize = x.n_rows;
  int ysize = y.n_rows;
  arma::mat res(xsize, ysize); 
  
  for(int ii = 0; ii < xsize; ++ii){
    for(int jj = 0; jj < ysize; ++jj) {
      res(ii, jj) = sqrt(sum(pow((x.row(ii) - y.row(jj)), 2)));
    }
  }
  return(res);
}

// [[Rcpp::export]]
field<mat> format_y (arma::colvec y, arma::colvec groupsTrain) {
  colvec uniqueGroupsTrain = unique(groupsTrain);
  int ntrain = uniqueGroupsTrain.n_elem;

  int ii; uvec ids;
  field<mat> ylist(ntrain, 1);

  for (ii = 0; ii < ntrain; ii++) {
    ids = find(groupsTrain == uniqueGroupsTrain(ii));
    ylist(ii, 0) = y.rows(ids);        
  }
     
  return(ylist);
}

// [[Rcpp::export]]
field<mat> format_x (arma::mat x, arma::colvec groupsTrain) {
  colvec uniqueGroupsTrain = unique(groupsTrain);
  int ntrain = uniqueGroupsTrain.n_elem;
  
  int ii; uvec ids;
  field<mat> xlist(ntrain, 1);

  for (ii = 0; ii < ntrain; ii++) {
    ids = find(groupsTrain == uniqueGroupsTrain(ii));
    xlist(ii, 0) = x.rows(ids);    
  }
  
  return(xlist);
}

// [[Rcpp::export]]
field<mat> format_z (arma::mat z, arma::colvec groupsTrain) {
  colvec uniqueGroupsTrain = unique(groupsTrain);
  int ntrain = uniqueGroupsTrain.n_elem;

  int ii; uvec ids;
  field<mat> zlist(ntrain, 1);

  for (ii = 0; ii < ntrain; ii++) {
    ids = find(groupsTrain == uniqueGroupsTrain(ii));
    zlist(ii, 0) = z.rows(ids);
  }
   
  return(zlist);
}

// [[Rcpp::export()]]
field<mat> sample_ran_eff (field<mat> ylist, field<mat> fixefList, field<mat> ranefList, 
                           mat lmat, double errVar, mat fixSamp) {

  int nsample = ylist.n_rows;
  int nran = ranefList(0, 0).n_cols;
  int ii;

  mat tmp1, tmp2, postRanVar, postRanMean;
  field<mat> ranSampList (nsample, 1);
  
  for (ii = 0; ii < nsample; ii++) {
    tmp1 = ranefList(ii, 0) * lmat;
    tmp2 = inv_sympd(tmp1 * trans(tmp1) + errVar * eye(tmp1.n_rows, tmp1.n_rows));
    postRanVar = eye(nran, nran) - (trans(tmp1) * tmp2) * tmp1;    
    postRanMean = trans(tmp1) * (tmp2 * (ylist(ii, 0) - fixefList(ii, 0) * fixSamp));
    ranSampList(ii, 0) = postRanMean + (chol(postRanVar, "lower") * randn<colvec>(nran));
  }

  return(ranSampList);
}

// List test_sample_ran_eff (mat y, mat x, mat z, 
//                           mat lmat, double errVar, mat fixSamp, mat groupsTrain) {

//   field<mat> ylist = format_y(y, groupsTrain);
//   field<mat> fixefList = format_x(x, groupsTrain);  
//   field<mat> ranefList = format_z(z, groupsTrain);      
  
//   colvec uniqueGroups = unique(groupsTrain);
 
//   int nsample = ylist.n_rows;
//   cout << nsample << endl;
  
//   int nran = ranefList(0, 0).n_cols;
//   int ii;


//   mat tmp1, tmp2, postRanVar, postRanMean;
//   List ranSampList(nsample);
//   for (ii = 0; ii < nsample; ii++) {
//     tmp1 = ranefList(ii, 0) * lmat;
//     tmp2 = inv_sympd(tmp1 * trans(tmp1) + errVar * eye(tmp1.n_rows, tmp1.n_rows));
//     postRanVar = eye(nran, nran) - (trans(tmp1) * tmp2) * tmp1;    
//     postRanMean = trans(tmp1) * (tmp2 * (ylist(ii, 0) - fixefList(ii, 0) * fixSamp));
//     ranSampList[ii] = join_vert(postRanMean, (vectorise(postRanVar)));
//   }

//   return(ranSampList);
// }

// [[Rcpp::export()]]
mat sample_fix_eff (field<mat> ylist, field<mat> fixefList, field<mat> ranefList, 
                    double npart, mat dmat, double errVar, mat muBeta0, mat sigBetaInv0) {
  
  int nsample = ylist.n_rows;
  int nfix = fixefList(0, 0).n_cols;
  int ii;
  
  mat umat, umatInv, xtux, xtuy;
  xtux = zeros<mat>(nfix, nfix);
  xtuy = zeros<mat>(nfix, 1);

  for (ii = 0; ii < nsample; ii++) {
    umatInv = errVar * eye(ylist(ii, 0).n_rows, ylist(ii, 0).n_rows) + (ranefList(ii, 0) * dmat) * trans(ranefList(ii, 0));
    umat = inv_sympd(umatInv);
    xtux = xtux + trans(fixefList(ii, 0)) * (umat * fixefList(ii, 0));
    xtuy = xtuy + trans(fixefList(ii, 0)) * (umat * ylist(ii, 0));
  }
  mat postFixVar = inv_sympd(sigBetaInv0 + npart * xtux); //chol2inv(chol(sigBetaInv0 + npart * sumXtux))  
  colvec postFixMean = (postFixVar * (trans(sigBetaInv0) * muBeta0 + npart * xtuy));

  return(postFixMean + chol(postFixVar, "lower") * randn<colvec>(nfix));
}


// List test_sample_fix_eff (mat y, mat x, mat z, 
//                           double npart, mat dmat, double errVar, mat muBeta0, mat sigBetaInv0, mat groupsTrain) {

//   field<mat> ylist = format_y(y, groupsTrain);
//   field<mat> fixefList = format_x(x, groupsTrain);  
//   field<mat> ranefList = format_z(z, groupsTrain);      
  
//   colvec uniqueGroups = unique(groupsTrain);
 
//   int nsample = ylist.n_rows;
//   cout << nsample << endl;
  
//   int nfix = fixefList(0, 0).n_cols;
//   int ii;
  
//   mat umat, umatInv, xtux, xtuy;
//   xtux = zeros<mat>(nfix, nfix);
//   xtuy = zeros<mat>(nfix, 1);

//   for (ii = 0; ii < nsample; ii++) {
//     umatInv = errVar * eye(ylist(ii, 0).n_rows, ylist(ii, 0).n_rows) + (ranefList(ii, 0) * dmat) * trans(ranefList(ii, 0));
//     umat = inv_sympd(umatInv);
//     xtux = xtux + trans(fixefList(ii, 0)) * (umat * fixefList(ii, 0));
//     xtuy = xtuy + trans(fixefList(ii, 0)) * (umat * ylist(ii, 0));
//   }
//   mat postFixVar = inv_sympd(sigBetaInv0 + npart * xtux); //chol2inv(chol(sigBetaInv0 + npart * sumXtux))  
//   colvec postFixMean = (postFixVar * (trans(sigBetaInv0) * muBeta0 + npart * xtuy));

//   List res;
//   res["mean"] = postFixMean;
//   res["var"] = postFixVar;
  
//   return(res);
// }

// [[Rcpp::export]]
mat sample_lmat (field<mat> ylist,
                 field<mat> fixefList,
                 field<mat> ranefList,
                 double npart, 
                 double errVar,
                 mat fixSamp,
                 field<mat> ranSampList, 
                 mat priorMeanL, mat priorCovInvL) {

  int nran = ranefList(0, 0).n_cols;
  int ndim = nran * (nran + 1) / 2;
  int nsample = ranefList.n_rows;
  int ii, jj, cts1, cts2;

  field<mat> zb(nsample, 1);
  for (ii = 0; ii < nsample; ii++) {
    zb(ii, 0) = zeros<mat>(ranefList(ii, 0).n_rows, ndim);
    cts1 = 0; cts2 = 0;
    for (jj = 0; jj < nran; jj++) {
      cts2 = cts2 + (nran - jj);
      zb(ii, 0).cols(cts1, cts2 - 1) = ranSampList(ii, 0)(jj, 0) * ranefList(ii, 0).tail_cols(nran - jj);
      cts1 = cts2;
    }
  }
  
  mat zTr = zeros<mat>(ndim, 1); //Z^T * resids
  mat zTz = zeros<mat>(ndim, ndim); //Z^T Z
  for (ii = 0; ii < nsample; ii++) {
    zTz = zTz + trans(zb(ii, 0)) * zb(ii, 0);
    zTr = zTr + trans(zb(ii, 0)) * (ylist(ii, 0) - fixefList(ii, 0) * fixSamp);
  }
    
  mat covL = inv_sympd(npart / errVar * zTz + priorCovInvL);
  mat meanL = covL * (npart / errVar * zTr + priorCovInvL * priorMeanL);
  
  mat lvec, L;
  
  lvec = meanL + chol(covL, "lower") * randn<colvec>(ndim);

  L = zeros<mat>(nran, nran);
  cts1 = -1;
  for (ii = 0; ii < nran; ii++) {
    for (jj = ii; jj < nran; jj++) {
      cts1 = cts1 + 1;
      L(jj, ii) = lvec(cts1);        
    }
  }

  return(L);
}
  
// List test_sample_lmat (mat y,
//                        mat x,
//                        mat z,
//                        double npart, 
//                        double errVar,
//                        mat fixSamp,
//                        mat ranSampListMat, 
//                        mat priorMeanL, mat priorCovInvL, mat groupsTrain) {
  
//   field<mat> ylist = format_y(y, groupsTrain);
//   field<mat> fixefList = format_x(x, groupsTrain);    
//   field<mat> ranefList = format_z(z, groupsTrain);    

//   int nran = ranefList(0, 0).n_cols;
//   int ndim = nran * (nran + 1) / 2;
//   int nsample = ranefList.n_rows;
//   int ii, jj, cts1, cts2;
  
//   field<mat> ranSampList(ranSampListMat.n_rows, 1);
//   for (ii = 0; ii < ranSampListMat.n_rows; ii++) {
//     ranSampList(ii, 0) = trans(ranSampListMat.row(ii));
//   }


//   field<mat> zb(nsample, 1);
//   for (ii = 0; ii < nsample; ii++) {
//     zb(ii, 0) = zeros<mat>(ranefList(ii, 0).n_rows, ndim);
//     cts1 = 0; cts2 = 0;
//     for (jj = 0; jj < nran; jj++) {
//       cts2 = cts2 + (nran - jj);
//       zb(ii, 0).cols(cts1, cts2 - 1) = ranSampList(ii, 0)(jj, 0) * ranefList(ii, 0).tail_cols(nran - jj);
//       cts1 = cts2;
//     }
//   }
  
//   mat zTr = zeros<mat>(ndim, 1); //Z^T * resids
//   mat zTz = zeros<mat>(ndim, ndim); //Z^T Z
//   for (ii = 0; ii < nsample; ii++) {
//     zTz = zTz + trans(zb(ii, 0)) * zb(ii, 0);
//     zTr = zTr + trans(zb(ii, 0)) * (ylist(ii, 0) - fixefList(ii, 0) * fixSamp);
//   }
    
//   mat covL = inv_sympd(npart / errVar * zTz + priorCovInvL);
//   mat meanL = covL * (npart / errVar * zTr + priorCovInvL * priorMeanL);
  
//   mat lvec, L;
  
//   lvec = meanL + chol(covL, "lower") * randn<colvec>(ndim);

//   L = zeros<mat>(nran, nran);
//   cts1 = -1;
//   for (ii = 0; ii < nran; ii++) {
//     for (jj = ii; jj < nran; jj++) {
//       cts1 = cts1 + 1;
//       L(jj, ii) = lvec(cts1);        
//     }
//   }
  
//   List res;
//   res["mean"] = meanL;
//   res["cov"] = covL;
  
//   return(res);
// }

// [[Rcpp::export()]]
double sample_err_var (field<mat> ylist, field<mat> fixefList, field<mat> ranefList, 
                       double npart, mat lmat, mat fixSamp, field<mat> ranSampList, 
                       double sig0, double nu0) {  

  int nsample = fixefList.n_rows;
  int ii;
  int nobs = 0;
  for (ii = 0; ii < nsample; ii++) {
    nobs = nobs + fixefList(ii, 0).n_rows;
  }
  
  double aa = 0.5 * (npart * nobs + nu0);
   
  double trm1 = 0.0; double trm2 = 0.0; mat resids;
  for (ii = 0; ii < nsample; ii++) {
    resids = ylist(ii, 0) - fixefList(ii, 0) * fixSamp - ranefList(ii, 0) * lmat * ranSampList(ii, 0);
    trm1 = trm1 + accu(resids % resids);
  }
  trm2 = sig0 * nu0;
  
  double bb = 0.5 * (npart * trm1 + trm2);
  // cout << "aa " << aa << " bb " << bb << " mean " << bb / (aa - 1) << " var " << bb * bb / ((aa - 1) * (aa - 1) * (aa - 2)) << endl;
  
  double res = 1 / R::rgamma(aa, 1 / bb); // note that the default parametrization in R and R::rgamma DO NOT MATCH!
 
  return(res);
}
  
// List test_sample_err_var (mat y, mat x, mat z, 
//                           double npart, mat lmat, colvec fixSamp, mat ranSampListMat, 
//                           double sig0, double nu0, colvec groupsTrain) {  

//   field<mat> ylist = format_y(y, groupsTrain);
//   field<mat> fixefList = format_x(x, groupsTrain);  
//   field<mat> ranefList = format_z(z, groupsTrain);      
  
//   colvec uniqueGroups = unique(groupsTrain);
 
//   int nsample = uniqueGroups.n_elem;
//   int ii;
//   int nobs = 0;
//   for (ii = 0; ii < nsample; ii++) {
//     nobs = nobs + fixefList(ii, 0).n_rows;
//   }
   
//   cout << "here1" << endl;
  
//   field<mat> ranSampList(ranSampListMat.n_rows, 1);
//   for (ii = 0; ii < ranSampListMat.n_rows; ii++) {
//     ranSampList(ii, 0) = trans(ranSampListMat.row(ii));
//   }
  
//   double aa = 0.5 * (npart * nobs + nu0);
   
//   double trm1 = 0.0; double trm2 = 0.0; mat resids;
//   for (ii = 0; ii < nsample; ii++) {
//     resids = ylist(ii, 0) - fixefList(ii, 0) * fixSamp - ranefList(ii, 0) * lmat * ranSampList(ii, 0);
//     trm1 = trm1 + accu(resids % resids);
//   }
//   cout << "trm1 " << trm1 << endl;
//   trm2 = sig0 * nu0;
//   cout << "trm2 " << trm1 << endl;  
  
//   double bb = 0.5 * (npart * trm1 + trm2);
//   // cout << "aa " << aa << " bb " << bb << " mean " << bb / (aa - 1) << " var " << bb * bb / ((aa - 1) * (aa - 1) * (aa - 2)) << endl;
  
//   List res; res["aa"] = aa; res["bb"] = bb; 
  
//   return(res);
// }
  
// [[Rcpp::export()]]
List wasp_lme_sampler (arma::mat y, arma::mat x, arma::mat z,
                       colvec groupsTrain,
                       double npart,
                       arma::mat priorBetaMean,
                       arma::mat priorBetaCovInv, 
                       mat dmat0,
                       double errVar0,
                       arma::mat priorLMean, arma::mat priorLCovInv,
                       double sig0, double nu0, 
                       int niter, int nburn, int nthin) {

  int nfix = x.n_cols;
  int nran = z.n_cols;
  int cts, its, ii, vv1, vv2, vv3; 
  int ndraw = (niter - nburn) / nthin; 
  wall_clock timer;

  double sErrVar = errVar0;        
  mat dmat = dmat0;
  mat sLmat = chol(dmat0, "lower"); 
  colvec sampErrVar = zeros<colvec>(ndraw); sampErrVar.fill(999.0); 
  mat sampBetas = mat(nfix, ndraw); sampBetas.fill(999.0);
  mat sampLmat = mat(nran * (nran + 1) / 2, ndraw); 
 
  mat sBetas = zeros<mat>(nfix, 1);
  
  field<mat> ylist = format_y(y, groupsTrain);
  field<mat> fixefList = format_x(x, groupsTrain);
  field<mat> ranefList = format_z(z, groupsTrain);
  field<mat> sRans(ylist.n_rows, 1); 
  
  cts = -1; timer.tic();   
  for (its = 0; its <= niter; its++) {
    
    if ((its % 500) == 0) cout << "DLS iteration: " << its << endl;  

    // I step
    sRans = sample_ran_eff(ylist, fixefList, ranefList, sLmat, sErrVar, sBetas);

    // P step        
    sBetas = sample_fix_eff(ylist, fixefList, ranefList, npart, dmat, sErrVar, priorBetaMean, priorBetaCovInv);
    sLmat = sample_lmat(ylist, fixefList, ranefList, npart, sErrVar, sBetas, sRans, priorLMean, priorLCovInv);
    dmat = sLmat * trans(sLmat);    
    sErrVar = sample_err_var(ylist, fixefList, ranefList, npart, sLmat, sBetas, sRans, sig0, nu0);
                  
    if (its > nburn && its % nthin == 0) {
      cts = cts + 1;      
      sampBetas.col(cts) = sBetas;
      sampErrVar(cts, 0) = sErrVar;
      
      vv3 = 0;
      for (vv1 = 0; vv1 < nran; vv1++) {
        sampLmat(vv3, cts) = sLmat(vv1, vv1);                    
        vv3 = vv3 + 1;
      }
      
      for (vv1 = 0; vv1 < nran; vv1++) {
        for (vv2 = vv1 + 1; vv2 < nran; vv2++) {
          sampLmat(vv3, cts) = sLmat(vv2, vv1);            
          vv3 = vv3 + 1;
        }
      }      
    }                           
    
  }
  double totTime = timer.toc();
  
  List res;
  res["beta"] = sampBetas.t();
  res["errVar"] = sampErrVar;
  res["lmat"] = sampLmat.t();  
  res["iter"] = its;
  res["time"] = totTime;  
  
  return(res);
}
