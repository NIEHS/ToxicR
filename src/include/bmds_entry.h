#pragma once

#ifndef bmds_nc_entryH
#define bmds_nc_entryH

#if defined WIN32 || defined _WINDOWS
#ifndef R_COMPILATION
#ifdef BMDS_MODELS_EXPORTS
#define BMDS_ENTRY_API __declspec(dllexport)
#else
#define BMDS_ENTRY_API __declspec(dllimport)
#endif // BMDS_MODELS_EXPORTS

typedef BSTR BMDS_MODEL_ID;
#else
#define BMDS_ENTRY_API
#ifndef _stdcall
#define _stdcall
#endif
typedef char *BMDS_MODEL_ID;
#endif // defined WIN32 || defined _WINDOWS

#else
#define BMDS_ENTRY_API
#define _stdcall
typedef char *BMDS_MODEL_ID;
#endif // defined WIN32 || defined _WINDOWS
#define BMDS_BLANK_VALUE -9999

static const int NUM_PRIOR_COLS = 5;

// The "pack(4)" pragma is required because Excel requires four-byte aligned structures,
// but Visual Studio uses eight-byte boundaries by default.

#ifndef R_COMPILATION
#pragma pack(4)
#endif

#define NUM_LIKELIHOODS_OF_INTEREST 5
#define NUM_TESTS_OF_INTEREST 4

// Indices for likelihoods of interest
#define LK_A1 0
#define LK_A2 1
#define LK_A3 2
#define LK_R 3
#define LK_FIT 4

enum VarType_t
{
  eVarTypeNone = 0,
  eConstant = 1,
  eModeled = 2
};

enum CModelID_t
{
  eExp2 = 2,
  eExp3 = 3,
  eExp4 = 4,
  eExp5 = 5,
  eHill = 6,
  ePoly = 7,
  ePow = 8
};

enum DModelID_t
{
  eDHill = 1,
  eGamma = 2,
  eLogistic = 3,
  eLogLogistic = 4,
  eLogProbit = 5,
  eMultistage = 6,
  eProbit = 7,
  eQLinear = 8,
  eWeibull = 9
};

enum BMDSPrior_t
{
  eNone = 0,
  eNormal = 1,
  eLognormal = 2
};

enum BMRType_t
{
  eAbsoluteDev = 1,
  eStandardDev = 2,
  eRelativeDev = 3,
  ePointEstimate = 4,
  eExtra = 5, // Not used
  eHybrid_Extra = 6,
  eHybrid_Added = 7
};

enum RiskType_t
{
  eExtraRisk = 1,
  eAddedRisk = 2
};

struct BMDS_C_Options_t
{
  double bmr;
  double alpha;
  double background;
  double tailProb; // Valid only for hybrid bmr type
  int bmrType;
  int degree;           // Valid for polynomial type models; for exponential, identifies the submodel
  int adverseDirection; // Direction of adversity: 0=auto, 1=up, -1=down
  int restriction;      // Restriction on parameters for certain models
  VarType_t varType;
  bool bLognormal;    // Valid only for continuous models
  bool bUserParmInit; // Use specified priors instead of calculated values
};

struct BMDS_D_Options_t
{
  double bmr;
  double alpha;
  double background;
  int bmrType;
  int degree; // Polynomial degree for the multistage model
};

struct BMDS_D_Opts1_t
{
  double bmr;
  double alpha;
  double background;
};

struct BMDS_D_Opts2_t
{
  int bmrType;
  int degree; // Polynomial degree for the multistage model
};

// Likelihoods of interest
struct LLRow_t
{
  double ll; // Log-likelihood
  double aic;
  int model;  // Data model number for test
  int nParms; // Count of model parameters
};

// Tests of interest (model deviance tests)
struct TestRow_t
{
  double deviance; // -2*log-likelihood ratio
  double pvalue;   // test p-value
  int testNumber;
  int df; // test degrees of freedom
};

// Indices for tests of interest
#define TI_1 0
#define TI_2 1
#define TI_3 2
#define TI_4 3

struct gofRow
{
  double dose;
  double estProb;  // Model-estimated probability for dose
  double expected; // Expected dose-response according to the model
  double observed;
  double size;
  double scaledResidual;
  double ebLower; // Error bar lower bound
  double ebUpper; // Error bar upper bound
};
typedef struct gofRow GoFRow_t;

typedef struct dGoF
{
  double chiSquare;
  double pvalue;
  GoFRow_t *pzRow;
  int df;
  int n; // Number of rows
} dGoF_t;

struct ContinuousDeviance_t
{
  LLRow_t *llRows;
  TestRow_t *testRows;
};

struct cGoFRow_t
{
  double dose;
  double obsMean;
  double obsStDev;
  double calcMedian;
  double calcGSD;
  double estMean; // Expected dose-response according to the model
  double estStDev;
  double size;
  double scaledResidual;
  double ebLower; // Error bar lower bound
  double ebUpper; // Error bar upper bound
};

struct cGoF_t
{
  double chiSquare;
  double pvalue;
  cGoFRow_t *pzRow;
  int df;
  int n; // Number of rows
};

struct BMD_C_ANAL
{
  BMDS_MODEL_ID model_id;
  double *PARMS;
  ContinuousDeviance_t deviance;
  cGoFRow_t *gofRow; // Goodness of Fit
  bool *boundedParms;
  double MAP;
  double BMD;
  double BMDL;
  double BMDU;
  double AIC;
  double BIC_Equiv; // BIC equivalent for Bayesian runs
  double ll_const;  // LL "additive" constant term
  double *aCDF;     // Array of cumulative density function values for BMD
  int nCDF;         // Requested number of aCDF elements to return
  int nparms;
  bool bAdverseUp;
};

// Result structure for individual dichotomous models
struct BMD_ANAL
{
  BMDS_MODEL_ID model_id;
  double MAP; // Equals the -LL for frequentist runs
  double BMD;
  double BMDL;
  double BMDU;
  double AIC;
  double BIC_Equiv; // BIC equivalent for Bayesian runs
  double *PARMS;
  double *aCDF; // Array of cumulative density function values for BMD
  struct DichotomousDeviance_t *deviance;
  dGoF_t *gof; // Goodness of Fit
  bool *boundedParms;
  int nparms;
  int nCDF; // Requested number of aCDF elements to return
  double *covM;
};

struct PRIOR
{
  double type; // 0= None (frequentist), 1=  normal (Bayesian), 2= log-normal (Bayesian)
  double initalValue;
  double stdDev; // Only used for type= 1 or 2
  double minValue;
  double maxValue;
};

// This holds model-specific details for model averaging
struct MA_ModelInfo
{
  double priorWeight;
  PRIOR *priors;
  int modelID;
  int degree;
  int nparms;
  int restriction;
  VarType_t varType;
  bool bLognormal;
  MA_ModelInfo()
  {
    bLognormal = false;
  }
};

struct MA_ModelOut
{
  double bmd;
  double bmdl;
  double bmdu;
  double post_prob;
};

struct MA_Results
{
  double avgBMD;
  double avgBMDL;
  double avgBMDU;
  MA_ModelOut *pModels;
};

struct MA_PRIORS
{
  double *mean_logistic;
  double *sd_logistic;
  double *mean_probit;
  double *sd_probit;
  double *mean_loglogit;
  double *sd_loglogit;
  double *mean_logprobit;
  double *sd_logprobit;
  double *mean_weibull;
  double *sd_weibull;
  double *mean_gamma;
  double *sd_gamma;
  double *mean_mult2;
  double *sd_mult2;
  double *mean_qlinear;
  double *sd_qlinear;
  double *mean_hill;
  double *sd_hill;
};

// Result structure for model averaging
struct MA_ANALYSIS
{
  double BMD;
  double BMDL;
  double BMDU;
  double *bmd;
  double *bmdl;
  double *bmdu;
  double *post_probs;
};

enum BMDSInputType_t
{
  unused = 0,
  eCont_2 = 1, // Individual dose-responses
  eCont_4 = 2, // Summarized dose-responses
  eDich_3 = 3, // Regular dichotomous dose-responses
  eDich_4 = 4  // Dichotomous d-r with covariate (e.g., nested)
};

struct BMDSInputData_t
{
  double dose;
  double response; // Mean value for summary data
  double groupSize;
  double col4; // stddev for cont_4 or covariate for dich_4
};

struct DichotomousDeviance_t
{
  double llFull;     // Full model log-likelihood
  double llReduced;  // Reduced model log-likelihood
  double devFit;     // Fit model deviance
  double devReduced; // Reduced model deviance
  double pvFit;      // Fit model p-value
  double pvReduced;  // Reduced model p-value
  int nparmFull;
  int nparmFit;
  int dfFit;
  int nparmReduced;
  int dfReduced;
};

#ifndef R_COMPILATION
#pragma pack()
#endif

// Declare public C interface functions to run analyses

extern "C"
{

  BMDS_ENTRY_API
  int _stdcall run_cmodel(CModelID_t *p_model, BMD_C_ANAL *returnV, BMDSInputType_t *p_inputType,
                          BMDSInputData_t *dataIn, PRIOR *priorsIn, BMDS_C_Options_t *options, int *p_n);

  BMDS_ENTRY_API
  int _stdcall run_dmodel2(DModelID_t *p_m_id, BMD_ANAL *returnV,
                           BMDSInputType_t *p_inputType, BMDSInputData_t *dataIn,
                           PRIOR *priorsIn, BMDS_D_Opts1_t *opt1, BMDS_D_Opts2_t *opt2, int *p_n);

  BMDS_ENTRY_API
  void _stdcall bmd_MA(BMDSInputType_t inputType, BMDSInputData_t *dataIn,
                       MA_PRIORS *priors, double *p_m_probs,
                       BMDS_D_Opts1_t *opt1, BMDS_D_Opts2_t *opt2, int n, double *post_p,
                       double *ma_bmd, double *bmd, double *bmdl, double *bmdu);

  BMDS_ENTRY_API
  int _stdcall run_cMA(BMDSInputType_t *p_inputType, BMDSInputData_t *dataIn,
                       MA_ModelInfo *modelInfo, int *p_nModels,
                       BMDS_C_Options_t *opts, int *p_nData, MA_Results *maOut);

} // extern "C"

#endif // ifndef bmds_nc_entryH
