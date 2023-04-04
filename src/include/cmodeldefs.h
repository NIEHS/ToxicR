

#ifndef CONTINUOUS_MODEL_DEFINE
#define CONTINUOUS_MODEL_DEFINE

typedef  int contbmd; 


#define CONTINUOUS_BMD_ABSOLUTE     1
#define CONTINUOUS_BMD_STD_DEV      2
#define CONTINUOUS_BMD_REL_DEV      3
#define CONTINUOUS_BMD_POINT        4
#define CONTINUOUS_BMD_EXTRA        5
#define CONTINUOUS_BMD_HYBRID_EXTRA 6
#define CONTINUOUS_BMD_HYBRID_ADDED 7
#define CONTINUOUS_BMD_EMPTY        0.0


#define PROFILE_INEQUALITY   1000
#define PROFILE_EQUALITY     2000


enum est_method {est_mle = 1, est_laplace=2, est_mcmc=3}; 
enum dich_model {d_hill =1, d_gamma=2,d_logistic=3, d_loglogistic=4,
                 d_logprobit=5, d_multistage=6,d_probit=7,
                 d_qlinear=8,d_weibull=9}; 

enum cont_model {generic = 0, hill = 6,exp_3 = 3,exp_5=5,power=8, funl = 10, polynomial = 666,
exp_aerts = 11, invexp_aerts = 12, gamma_aerts = 13, invgamma_aerts = 14, hill_aerts = 15, lomax_aerts = 16,
invlomax_aerts = 17, lognormal_aerts = 18, logskew_aerts = 19, invlogskew_aerts = 20, logistic_aerts = 21,
probit_aerts = 22, LMS = 23, gamma_efsa = 24};
enum distribution {normal = 1, normal_ncv = 2, log_normal = 3}; 

enum prior_iidtype {iid_normal = 1, iid_lognormal = 2, iid_mle = 0,
					iid_cauchy = 3, iid_gamma = 4, iid_pert = 5};

#endif
