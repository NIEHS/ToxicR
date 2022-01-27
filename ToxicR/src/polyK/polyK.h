/*
 * $URL: https://svn2.niehs.nih.gov/svn/ntp-cus/dev/tdmse/branches/b_2.5.0.0/src/gov/nih/niehs/tdmse/reporting/pathology/P8_Calculator.java $
 * National Toxicology Program (NIEHS)
 */


/**
 * Calculator helper class for P8 and P10 reports
 * @author $Author: mercerm $
 * 
 *   
 public static void main(String[] args)  {
   TDMSE_PolyK calc = new TDMSE_PolyK();
   PolyKPrepareClass subsetVars = new PolyKPrepareClass();
   subsetVars.setFilename("C:\\Data\\TDMSE Custom\\PolyK Test\\BCE Male Liver HepAd PW4.txt");
   subsetVars.setDebugFlag(false);   // Change to true to show debug messages
   subsetVars.prepare();
   int ntrt = subsetVars.getNumDoseLevels();
   double pval = calc.polyk_mod(subsetVars, ntrt, POLY_3_TEST, -1.0);
   System.out.println("Poly-3 P-value = " + pval);
   pval = calc.polyk_mod(subsetVars, ntrt, POLY_1PT5_TEST, -1.0);
   System.out.println("Poly-1.5 P-value = " + pval);
   pval = calc.polyk_mod(subsetVars, ntrt, POLY_6_TEST, -1.0);
   System.out.println("Poly-6 P-value = " + pval);
  }
 * @version $Revision: 9150 $ $Date: 2011-05-16 14:08:23 -0400 (Mon, 16 May 2011) $
 * @since 1.0
 */
#include <math.h>
#include <limits>
#include "polyK_setup.h"

#ifndef _POLYK_H_DEF
#define _POLYK_H_DEF
namespace PolyK
{

static  double             BAD_PVALUE                      = 3.2;
static double                    PVALUE_NO_STATS_CALCULATED      = 2.8;

// Global constants
static const int               POLY_3_TEST                     = 0;
static const int               POLY_6_TEST                     = 1;
static const int               POLY_1PT5_TEST                  = 2;
static double            MIN_PVALUE                      = 0.001;
// Constants for gser and gcf functions
static int               ITMAX                           = 100;
static double            EPS                             = 3.0e-7;
// Constants for gammln function
static double            COF[]                           = { 76.18009173e0,
                                                                     -86.50532033e0,
                                                                     24.01409822e0,
                                                                     -1.231739516e0,
                                                                     0.120858003e-2,
                                                                     -.536382e-5};
static double            STP                             = 2.50662827465e0;
static double            FPF                             = 5.5;
static int               NUM_COF                         = 6;
// Constants for poly-k test
static double            POLYK_EPSILON                   = 1e-15;
// Stupid JTest constants
static double            TESTSTAT_EPSILON                = 1e-10;

class TDMSE_PolyK
{
    public: 
      
    double polyk_mod(PolyKPrepareClass subsetVars,
                               int ntrt,
                               int polyKTestType,
                               double cochArmFishPValue); 
 
  
  /**
   * Calculate the chisq probability for the given test statistic and degrees of freedom
   */
private: 
  double chisq(double x2, int idf);

  /**
   * Helper function for chisq calculation
   */
  double gammp(double a, double x);

  /**
   * Helper function for chisq calculation
   */
  double gser(double a, double x);

  /**
   * Helper function for chisq calculation
   */
  double gcf(double a, double x);

  /**
   * Helper function for chisq calculation
   */
  double gammln(double xx);

};

}
#endif
