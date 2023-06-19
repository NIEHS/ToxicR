#include <analysis_of_deviance.h>
#include <mcmc_analysis.h>
#include <bmd_calculate.h>
#include <IDPriorMCMC.h>
 
#include <chrono>
#include <algorithm>
#include <vector>
#include <limits>
#include <set>
#include <numeric>

void inverse_transform_dose(continuous_model_result *model){
  if (model){
    model->bmd = sinh(model->bmd); 
    for (int i = 0; i< model->dist_numE; i ++){
      double temp = double(i)/double(model->dist_numE); 
      model->bmd_dist[i] = sinh(model->bmd_dist[i]);     // BMD @ probability
    }
  }
}

void inverse_transform_dose(bmd_analysis_MCMC *b){
  
  if (b){
    for(unsigned int i= 0; i < b->samples;  i++){
      b->BMDS[i] = sinh(b->BMDS[i]); 
    }
  }
  
}

Eigen::MatrixXd quadraticRegression(Eigen::MatrixXd Y_N, Eigen::MatrixXd X){
  
  Eigen::MatrixXd mX = Eigen::MatrixXd::Zero(Y_N.rows(), 3); 
  Eigen::MatrixXd W  = Eigen::MatrixXd::Zero(Y_N.rows(), Y_N.rows());
  for (int i = 0; i < mX.rows(); i++)
  {  
    W(i, i) = Y_N.cols() == 3? pow(1/Y_N(i,1),2) * Y_N(i,2) : 1;
    for (int j = 0; j < 3; j++) {
      switch(j){
      case 2:
        mX(i, j) = X(i,0)*X(i,0);
        break; 
      case 1:
        mX(i, j) = X(i,0); 
        break; 
      default: 
        mX(i, j) = 1 ;
      break;  
      } 
    }
  }

  Eigen::MatrixXd betas = mX.transpose()*W*mX;
  betas = betas.inverse()*mX.transpose()*W*Y_N.col(0);
  return betas;
}

Eigen::MatrixXd powerSearchRegression(Eigen::MatrixXd Y_N, Eigen::MatrixXd X){
  
  double min = 0; 
  double sum ; 
  
  Eigen::MatrixXd mX = Eigen::MatrixXd::Zero(Y_N.rows(), 2); 
  Eigen::MatrixXd W  = Eigen::MatrixXd::Zero(Y_N.rows(), Y_N.rows());
  Eigen::MatrixXd E  = mX; 
  Eigen::MatrixXd betas;
  Eigen::MatrixXd rbetas(3,1);
  
  for (double pows = 1.0; pows < 17; pows += 0.5){
 
    for (int i = 0; i < mX.rows(); i++)
    {  
      W(i, i) = Y_N.cols() == 3? pow(1/Y_N(i,1),2) * Y_N(i,2) : 1;
      for (int j = 0; j < 2; j++) {
        switch(j){
        case 1:
          mX(i, j) = pow(X(i,0),pows);
          break; 
        default: 
          mX(i, j) = 1 ;
          break;  
        } 
      }
    }
    
    betas = mX.transpose()*W*mX;
    betas = betas.inverse()*mX.transpose()*W*Y_N.col(0);
    
    E = (Y_N.col(0) - mX * betas);
    E = E.array() * E.array(); 
    E = W*E; 
    sum = E.array().sum();
    
    if (pows == 1.0 || sum < min){
      rbetas(0,0) = betas(0,0);  
      rbetas(1,0) = betas(1,0); 
      rbetas(2,0) = pows; 
      min    =  sum; 
    }

  }
  return rbetas;
}

Eigen::MatrixXd init_funl_nor(Eigen::MatrixXd Y_N, Eigen::MatrixXd X, Eigen::MatrixXd prior){
  
  
  std::vector<double> vec(X.data(), X.data() + X.rows() * X.cols());
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()),	vec.end());
  //udoses = vec; // this should be the unique dose group
  
  Eigen::MatrixXd betas = quadraticRegression(Y_N,  X);
  prior(0,1) = betas(0,0); 
  double max_d = vec[vec.size()-1]; 
  double max_r = (betas(0,0)+betas(1,0)*max_d + betas(2,0)*max_d*max_d);
  prior(1,1)   = (betas(0,0)+betas(1,0)*max_d + betas(2,0)*max_d*max_d - prior(0,1))/max_d; 
  prior(2,1)   = max_r; 
  prior(3,1)   = 0.5;
  prior(4,1)   = 1;
  prior(5,1)   = 0.75;
  prior(6,1)   = 1;
  
  for (int i = 0; i < 7; i++){
    if (prior(i,1) < prior(i,3)) prior(i,1) = prior(i,3); 
    if (prior(i,1) > prior(i,4)) prior(i,1) = prior(i,4);
  }
  
  
  return prior; 
}

Eigen::MatrixXd init_test4_nor(Eigen::MatrixXd Y_N, Eigen::MatrixXd X, Eigen::MatrixXd prior){

  double minDose = X.minCoeff();
  double maxDose = X.maxCoeff();
  double init = 0;
  int nmin = 0, nmax = 0;

  for (int i = 0; i < X.rows(); i++){
    if (X(i,0)==minDose){
      nmin++;
      init += Y_N(i,0);
    }
  }
  init *= init/double(nmin);
  prior(0,1) = init;

  init = 0;
  for (int i = 0; i < X.rows(); i++){
    if (X(i,0)==maxDose){
      nmax++;
      init += Y_N(i,0);
    }
  }
  init *= init/double(nmax);

  prior(2,1)   =  init / prior(0,1);
  prior(1,1)   = 0.0001*maxDose;
  prior(3,1)   = 5;

  //make sure the starting point is within bounds; if not, put on boundary
  for(int i = 0; i < 4; i++){
	  if (prior(i,1) < prior(i,3)) prior(i,1) = prior(i,3);
	  if (prior(i,1) > prior(i,4)) prior(i,1) = prior(i,4);
  }
  //cerr << prior << endl;
  return prior;
}

Eigen::MatrixXd init_test4_lognor(Eigen::MatrixXd Y_LN, Eigen::MatrixXd X, Eigen::MatrixXd prior){
  Y_LN.col(0) = exp(Y_LN.col(0).array());
  if (Y_LN.cols() ==3 ){
    Y_LN.col(1) = exp(Y_LN.col(1).array());
  }
  return init_test4_nor(Y_LN,  X, prior);

}

Eigen::MatrixXd init_test5_nor(Eigen::MatrixXd Y_N, Eigen::MatrixXd X, Eigen::MatrixXd prior){

  double minDose = X.minCoeff();
  double maxDose = X.maxCoeff();
  double init = 0;
  int nmin = 0, nmax = 0;

  for (int i = 0; i < X.rows(); i++){
    if (X(i,0)==minDose){
      nmin++;
      init += Y_N(i,0);
    }
  }
  init *= init/double(nmin);
  prior(0,1) = init;

  init = 0;
  for (int i = 0; i < X.rows(); i++){
    if (X(i,0)==maxDose){
      nmax++;
      init += Y_N(i,0);
    }
  }
  init *= init/double(nmax);

  prior(2,1)   =  init / prior(0,1);
  prior(1,1)   = 0.0001*maxDose;
  prior(3,1)   = 5;
  prior(4,1)   = 2;

  //make sure the starting point is within bounds; if not, put on boundary
  for(int i = 0; i < 5; i++){
	  if (prior(i,1) < prior(i,3)) prior(i,1) = prior(i,3);
	  if (prior(i,1) > prior(i,4)) prior(i,1) = prior(i,4);
  }
  //cerr << prior << endl;
  return prior;
}

Eigen::MatrixXd init_test5_lognor(Eigen::MatrixXd Y_LN, Eigen::MatrixXd X, Eigen::MatrixXd prior){
  Y_LN.col(0) = exp(Y_LN.col(0).array());
  if (Y_LN.cols() ==3 ){
    Y_LN.col(1) = exp(Y_LN.col(1).array());
  }
  return init_test5_nor(Y_LN,  X, prior);

}


Eigen::MatrixXd init_hill_nor(Eigen::MatrixXd Y_N, Eigen::MatrixXd X, Eigen::MatrixXd prior){
  
  
  std::vector<double> vec(X.data(), X.data() + X.rows() * X.cols());
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()),	vec.end());
  //udoses = vec; // this should be the unique dose group
  double minDose = X.minCoeff(); 
  double maxDose = X.maxCoeff(); 
  double init = 0; 
  int nmin = 0, nmax = 0;  
  
  for (int i = 0; i < X.rows(); i++){
    if (X(i,0)==minDose){
      nmin++; 
      init += Y_N(i,0); 
    }
  }
  init *= init/double(nmin); 
  
  Eigen::MatrixXd betas = quadraticRegression(Y_N,  X);
  prior(0,1) = init; 
  init = 0;
  for (int i = 0; i < X.rows(); i++){
    if (X(i,0)==maxDose){
      nmax++; 
      init += Y_N(i,0); 
    }
  }
  init *= init/double(nmin); 
  
  prior(1,1)   =  (init - prior(0,1))/(maxDose-minDose); 
  prior(2,1)   = 0;//0.0001*maxDose; 
  prior(3,1)   = 10;
  
  if (prior(0,1) < prior(0,3)) prior(0,1) = prior(0,3); 
  if (prior(0,1) > prior(0,4)) prior(0,1) = prior(0,4);
  
  if (prior(1,1) < prior(1,3)) prior(1,1) = prior(1,3); 
  if (prior(1,1) > prior(1,4)) prior(1,1) = prior(1,4);
  //cerr << prior << endl; 
  return prior; 
}


Eigen::MatrixXd init_pow_nor(Eigen::MatrixXd Y_N, Eigen::MatrixXd X, Eigen::MatrixXd prior){
  
  
  std::vector<double> vec(X.data(), X.data() + X.rows() * X.cols());
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()),	vec.end());
  //udoses = vec; // this should be the unique dose group
  
  Eigen::MatrixXd betas = powerSearchRegression(Y_N,  X);
  prior(0,1)   = betas(0,0); 
  prior(1,1)   = betas(1,0);  
  prior(2,1)   = betas(2,0);
  
  if (prior(0,1) < prior(0,3)) prior(0,1) = prior(0,3); 
  if (prior(0,1) > prior(0,4)) prior(0,1) = prior(0,4);
  
  if (prior(1,1) < prior(1,3)) prior(1,1) = prior(1,3); 
  if (prior(1,1) > prior(1,4)) prior(1,1) = prior(1,4);
  
  if (prior(2,1) < prior(1,3)) prior(2,1) = prior(1,3); 
  if (prior(2,1) > prior(1,4)) prior(2,1) = prior(1,4);
  
  return prior; 
}

Eigen::MatrixXd init_hill_lognor(Eigen::MatrixXd Y_LN, Eigen::MatrixXd X, Eigen::MatrixXd prior){
  Y_LN.col(0) = exp(Y_LN.col(0).array());
  if (Y_LN.cols() ==3 ){
    Y_LN.col(1) = exp(Y_LN.col(1).array());
  }
  return init_hill_nor(Y_LN,  X, prior); 
  
}


Eigen::MatrixXd init_exp_nor(Eigen::MatrixXd Y_N, Eigen::MatrixXd X, Eigen::MatrixXd prior){
  
  std::vector<double> vec(X.data(), X.data() + X.rows() * X.cols());
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()),	vec.end());
  //udoses = vec; // this should be the unique dose group
  
  Eigen::MatrixXd betas = quadraticRegression(Y_N,  X);
  prior(0,1) = betas(0,0); 
  double max_d = vec[vec.size()-1]; 
  double max_r = (betas(0,0)+betas(1,0)*max_d + betas(2,0)*max_d*max_d);
  prior(2,1)   = log(0.001);  
  double temp = max_r/prior(0,1);
  
  temp =  -(temp-exp(prior(2,1)))/(exp(prior(2,1))-1.0);
  
  prior(1,1) = 0.05; 
  prior(3,1)   = 2.5; 
  
  if (prior(0,1) < prior(0,3)) prior(0,1) = prior(0,3); 
  if (prior(0,1) > prior(0,4)) prior(0,1) = prior(0,4);
  
  if (prior(1,1) < prior(1,3)) prior(1,1) = prior(1,3); 
  if (prior(1,1) > prior(1,4)) prior(1,1) = prior(1,4);
  
  
  return prior; 
}

Eigen::MatrixXd init_exp_lognor(Eigen::MatrixXd Y_LN, Eigen::MatrixXd X, Eigen::MatrixXd prior){
 //right here
  Y_LN.col(0) = exp(Y_LN.col(0).array());
  if (Y_LN.cols() ==3 ){
      Y_LN.col(1) = exp(Y_LN.col(1).array());
  }
  return init_exp_nor(Y_LN, X, prior); 
}

Eigen::MatrixXd init_poly(Eigen::MatrixXd Y, Eigen::MatrixXd tX, 
                          Eigen::MatrixXd prior, int deg = 2){

  Eigen::MatrixXd X = Eigen::MatrixXd::Ones(tX.rows(),deg+1);
  Eigen::MatrixXd W = Eigen::MatrixXd::Identity(tX.rows(),tX.rows());
 
  for (int i = 0; i < X.rows(); i++){
    if (Y.cols()>1){
      W(i,i) = Y(i,2)/Y(i,1)*Y(i,1); // Weights: \sigma^2/N
    }
    for (int j = 1; j < X.cols(); j++){
      X(i,j) = pow(tX(i,0),j);  
    }
  }
  Eigen::MatrixXd B = Eigen::MatrixXd::Ones(deg+1,1);
  B = X.transpose()*W*X;
  B = B.inverse()*X.transpose()*W*Y.col(0); 
  for(int i = 0; i < B.rows(); i++){
    if ( B(i,0) < prior(i,3) ){
      prior(i,1) = prior(i,3); 
    } else if (B(i,0) > prior(i,4)){
      prior(i,1) = prior(i,4); 
    }else{
      prior(i,1) = B(i,0);
    }
  } 
  
  return prior; 
}

bool convertSStat(Eigen::MatrixXd Y, Eigen::MatrixXd X,
                  Eigen::MatrixXd *SSTAT, Eigen::MatrixXd *SSTAT_LN,
                  Eigen::MatrixXd *UX){
  bool convert = true; 
  
  if (Y.cols() == 1 ){
    // check to see if it can be converted into sufficient statistics4
    int temp = 0; 
    // go through each row to see if there are duplicates
    for (int i = 0; i < X.rows(); i++){
      for (int j = 0 ; j < X.rows(); j++){
        if (X(i,0) == X(j,0)){
          temp++; 
        }
      }
      if( temp == 1){
        convert = false; 
      }
      temp = 0; 
    }
    
    if (convert){
      // we can convert the data
       *SSTAT    = createSuffStat( Y, X, false);
       *SSTAT_LN = createSuffStat(Y,X,true); 
        std::vector<double> uniqueX = unique_list(X );
        Eigen::MatrixXd temp_X(uniqueX.size(),1);
        for (int i = 0; i < uniqueX.size(); i++){
          temp_X(i,0) = uniqueX[i]; 
        }
        *UX = temp_X; 
        return true; 
 
     }
  
  }else{
    *SSTAT    = createSuffStat( Y, X, false);
    *SSTAT_LN = createSuffStat(Y,X,true); 
  }
  
  
    
  return false; 
}
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
  unsigned int numRows = matrix.rows()-1;
  unsigned int numCols = matrix.cols();
  
  if( rowToRemove < numRows )
    matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);
  
  matrix.conservativeResize(numRows,numCols);
}


void removeCol(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols()-1;
  
  if( colToRemove < numCols )
    matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);
  
  matrix.conservativeResize(numRows,numCols);
}

bmd_analysis create_bmd_analysis_from_mcmc(unsigned int burnin, mcmcSamples s,double max_d){
  bmd_analysis rV;
  rV.MAP          = s.map; 
  rV.MAP_ESTIMATE = s.map_estimate; 
  rV.COV          = s.map_cov; 
  rV.MAP_BMD      = 0; 
  int total = 0; 
  int bad   = 0; 
 
  std::vector<double> v; 
  for (int i = burnin; i < s.BMD.cols(); i++){
    total ++;
    if ( isfinite(s.BMD(0,i)) && s.BMD(0,i) < 10*max_d){ // always look at 5x max dose tested
	        v.push_back(s.BMD(0,i));   // get rid of the burn in samples
    }else{
      bad++; 
    }
  }
  
  double sum = std::accumulate(v.begin(), v.end(), 0.0); 
//= sum/v.size(); //average of the non-infinite bmds
  std::vector<double>  prob;
  std::vector<double> bmd_q; 
  double inf_prob =  double(bad)/double(total); // bad observations 
  if (v.size() > 0){
    std::sort(v.begin(), v.end());
   for (double k = 0.004; k <= 0.9999; k += 0.005){
    	    prob.push_back(k*(1.0-inf_prob)); 
          int idx = int(k*double(v.size()));
          idx = idx == 0? 0: idx-1; 
    	    bmd_q.push_back(v[idx]);
   }
 
    // fix numerical quantile issues.
    for (int i = 1; i < bmd_q.size(); i++){
      if (bmd_q[i] <= bmd_q[i-1]){
    	  bmd_q[i] = bmd_q[i-1] + 1e-6;
//         for (int kk = i; kk <  bmd_q.size(); kk++){
//            bmd_q[kk] = bmd_q[kk-1] + 1e-6;
//         }
      } 
 
    }
  
    if (prob.size() > 10 && *min_element(bmd_q.begin(), bmd_q.end())  < 1e8 
                         && bmd_q[0] > 0 ){  
        rV.BMD_CDF = bmd_cdf(prob,bmd_q);
    }
    rV.MAP_BMD = rV.BMD_CDF.inv(0.5/(1.0-inf_prob));

  }
   // approximate median; 
  return rV; 
}

void transfer_mcmc_output(mcmcSamples a, bmd_analysis_MCMC *b){

  if (b){
    b->samples = a.samples.cols(); 
    b->nparms  = a.samples.rows(); 

    for(unsigned int i= 0; i < a.BMD.cols();  i++){
      b->BMDS[i] = a.BMD(0,i); 
      for(unsigned int j = 0; j < a.samples.rows(); j++){
        b->parms[i + j*a.BMD.cols()] = a.samples(j,i); 
      }
    }
  }

}

void bmd_range_find(continuousMA_result *res, 
					double *range){
 // assume the minimum BMD for the MA is always 0
 range[0] = 0.0; 
 double current_max = 0.0; 
 for (int j = 10; j > 1; j--){
	 for (int i = 0; i < res->nmodels;i++){
		int temp_idx = res->models[i]->dist_numE -j; 
		
		// make sure we are not dealing with an infinite value
		// or not a number
		if (isfinite(res->models[i]->bmd_dist[temp_idx]) && 
			  !isnan(res->models[i]->bmd_dist[temp_idx])){
			if ( res->models[i]->bmd_dist[temp_idx] > current_max){
				current_max = res->models[i]->bmd_dist[temp_idx]; 
			}
		}
		 
	 }
 }
 // if we don't find a max then the upper limit is NAN
 range[1] = current_max == 0.0 ? std::numeric_limits<double>::quiet_NaN():current_max; 
  
}

void transfer_continuous_model(bmd_analysis a, continuous_model_result *model){
	if (model){
	  model->nparms = a.COV.rows(); 
		model->max = a.MAP; 
		model->bmd = a.MAP_BMD; 
		for (int i = 0; i< model->dist_numE; i ++){
				double temp = double(i)/double(model->dist_numE); 
				model->bmd_dist[i] = a.BMD_CDF.inv(temp);     // BMD @ probability
				model->bmd_dist[model->dist_numE + i] = temp; // probability 
		}
		for (int i = 0; i < model->nparms; i++){
				model->parms[i] = a.MAP_ESTIMATE(i,0); 
				for (int j = 0; j < model->nparms; j++){
					model->cov[i + j*model->nparms] = a.COV(i,j); 
				}
		}
	}
}