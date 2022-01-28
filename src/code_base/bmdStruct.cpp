#ifdef R_COMPILATION
//necessary things to run in R    
  #include <RcppEigen.h>
  #include <RcppGSL.h>
#else 
  #include <Eigen/Dense>
#endif

#include "bmdStruct.h"

void del_continuous_analysis(continuous_analysis a){
  if (a.Y)  delete[] a.Y; 
  if (a.n_group) delete[] a.n_group; 
  if (a.sd) delete[] a.sd; 
  if (a.doses) delete[] a.doses; 
  if (a.prior) delete[] a.prior; 
  a.Y = NULL;  
  a.n_group = NULL; 
  a.sd = NULL; 
  a.doses = NULL; 
  a.prior = NULL; 
}


dichotomous_model_result * new_dichotomous_model_result(int model,
                                                        int parms,
                                                        int dist_numE){
  
  dichotomous_model_result *dmodel = new dichotomous_model_result; 
  dmodel->model = model; 
  dmodel->nparms = parms; 
  dmodel->max = -1.0*std::numeric_limits<double>::infinity(); 
  dmodel->dist_numE = dist_numE; 
  dmodel->parms = new double[parms]; 
  dmodel->cov   = new double[parms*parms]; 
  dmodel->bmd_dist = new double[dist_numE*2]; 
 
  return dmodel;  
}

void delete_dichotomous_model_result(dichotomous_model_result * dmodel)
                                                           {
  if (dmodel){
    delete[] dmodel->parms; delete[] dmodel->cov;   
    delete[] dmodel->bmd_dist; 
    delete dmodel;
  }
  return ;  
}

dichotomousMA_result * new_dichotomousMA_result(int nmodels,
                                                int dist_numE)
{
  dichotomousMA_result * dres = new dichotomousMA_result; 
  dres->nmodels = nmodels; 
  dres->models  = new  dichotomous_model_result * [nmodels];
  dres->dist_numE = dist_numE; 
  dres->post_probs = new double[nmodels];
  dres->bmd_dist   = new double[dist_numE*2];

  return dres; 
}
  
  
  
void delete_dichotomousMA_result(dichotomousMA_result *res){
  if (res){
    for (int i = 0; i < res->nmodels; i++){
       if (res->models[i]){
          delete_dichotomous_model_result(res->models[i]);
       }
    }
    
    delete[] res->models; 
    delete[] res->post_probs;
    delete[] res->bmd_dist;
    delete res; 
  }
  
  return; 
}
  
/////////////////////////////////////////////////////////
bmd_analysis_MCMC * new_mcmc_analysis(int model,
                                       int parms, 
                                       unsigned int samples){
  bmd_analysis_MCMC *rV = new bmd_analysis_MCMC; 
  rV->model = model; 
  rV->samples = samples; 
  rV->BMDS    = new double[samples]; 
  rV->parms   = new double[samples*parms]; 
  
  return rV; 
  
}



void del_mcmc_analysis(bmd_analysis_MCMC *an){
  if (an){
    delete[] an->BMDS; 
    delete[] an->parms; 
    delete an; 
  }
}

continuous_model_result * new_continuous_model_result(int model,
													  unsigned int n_parm,
                                                      unsigned int n_elm){
  continuous_model_result *rV = new continuous_model_result; 
  rV->dist_numE = n_elm; 
  rV->cov = new double[n_parm*n_parm];
  rV->parms = new double[n_parm*n_parm]; 
  rV->nparms = n_parm; 
  rV->bmd_dist = new double[n_elm*2]; // n_elm x 2 
  return rV; 
}

void del_continuous_model_result(continuous_model_result * cm){
	
	if(cm){	
		if (cm->cov){
			delete[] cm->cov; 
		}
		if (cm->parms){
			delete[] cm->parms; 
		}
		if (cm->bmd_dist){
			delete[] cm->bmd_dist; 
		}
		delete cm; 
	}
	
}

void del_continuousMA_analysis(continuousMA_analysis CMA){
  if (CMA.priors){
    for (int i = 0; i < CMA.nmodels; i++){
      delete[] CMA.priors[i]; 
    } 
    delete[] CMA.nparms;
    delete[] CMA.modelPriors; 
    delete[] CMA.priors; 
    delete[] CMA.actual_parms; 
    delete[] CMA.models; 
    delete[] CMA.disttype; 
    delete[] CMA.prior_cols; 
  }
}

void cp_prior(Eigen::MatrixXd temp, double *priors){
  
  for (int i = 0; i < temp.rows(); i++){
    for (int j = 0; j < temp.cols(); j++){
      priors[i+j*temp.rows()] = temp(i,j);
    }
  }
  return; 
}
