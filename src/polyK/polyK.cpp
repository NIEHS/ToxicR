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
#include <cmath>
#include "polyK.h"
namespace PolyK
{
/**************************************************************************/

/**************************************************************************/
/**
 * Compute the p-value for the Poly-k test of Bailer and Portier
 * The Equation is:
 *
 *     T1 = (ABS(A) - CF) / B
 *
 *  where A = Sum(a*d*Pi) - (Sum(a*d)*Sum(a*Pi)/Sum(a)
 *        B = S * sqrt(Sum(a*d*d) - (Sum(a*d))*(Sum(a*d))/Sum(a)
 *        CF = 0.5 * Max of (nStar[i]/n[i] * (d[i] - dBar) - nStar[i-1]/n[i-1] * [d[i-1] - dBar))
 *
 *  where dBar = Sum (a*d) / Sum(a)
 *        a[i] = nStar[i] * nStar[i] / n[i]
 *        d[i] = dose value at i
 *        nStar[i] = poly-k adjusted sample size (= Sum(t[i,j]**3)  where t = fraction of duration of study)
 *        Pi[i] = n[i]/nStar[i]
 *        S = sqrt(Sum(Sum((r[i,j]-rBar[i])**2)) / Sum(n[i]-1))
 *        r[i,j] = y[i,j] - piHat * w[i,j]
 *        piHat[i] = Sum(y[i,j])/nStar[i]
 *
 *
 * Input variable are:
 *   subsetVars = Input data summarized in terms of number of animals with and without tumors for each day one or more animals were removed
 *   ntrt = Number of groups, including control
 *   polyKTestType = One of POLY_3_TEST, POLY_6_TEST, POLY_1PT5_TEST
 *   cochArmFishPValue = P-value from the Cochran-Armitage test (rarely if ever used)
 */

double TDMSE_PolyK::polyk_mod(PolyKPrepareClass subsetVars,
                              int ntrt,
                              int polyKTestType,
                              double cochArmFishPValue)
{
  bool debugFlag = false;   // Change to true to show debug messages
  
  double pvalue = BAD_PVALUE;
  double dBar;
  double dBartop;
  double dBarbot;
  double front = 0.0;
  double back = 0.0;
  double cf = 0.0;
  double sum_adpi;
  double sum_ad;
  double sum_api;
  double sum_a;
  double sum_add;
  double sum6;
  double sum7;
  double sum8;
  double sum9;
  double sum10;
  double sum11;
  double y_n;
  double y_nk;
  double y_ak;
  double y_pk;
  double y_dk;
  double base;
  double sumj_bn0;
  double sumj_am3;
  double sumj_am6;
  double xbn, xam = 0;
  double top, den_p, p, q, bot;
  double term1;
  double term2;
  double term3;
  double topcc;
  double test_stat;
  std::vector<double> ai(ntrt+1);// =    new double[ntrt + 1];
  std::vector<double> ai2(ntrt + 1);// =   new double[ntrt + 1];
  std::vector<double> scale(ntrt + 1);// = new double[ntrt + 1];
  double timeValue = 0;
  double exp1 = 0;
  double exp2 = 0;
  
  int quitnow;
  int ndf;
  dBar = 0;
  dBartop = 0;
  dBarbot = 0;

  try
  {
    for(int ib = 1; ib <= ntrt; ib++)
    {
      ai[ib] = 0.0;
      ai2[ib] = 0.0;
      scale[ib] = subsetVars.getWeight(ib);   /* These are the dose levels scaled from 0 to 1 */
    }
    
    if (ntrt == 2)
    {
      scale[1] = 0.0;
      scale[2] = 1.0;
    }
    
    if(subsetVars.getN_subj(1) == 0)
    {
      pvalue = BAD_PVALUE;
      return pvalue;
    }
    
    /* Check for conditions where no calculation is performed */
    quitnow = 0;
    for(int dcheck = 1; dcheck <= ntrt; dcheck++)
    {
      if(subsetVars.getN_subj(dcheck) == 0) quitnow = 1;
    }
    
    if(quitnow == 1)
    {
      pvalue = BAD_PVALUE;
      return pvalue;
    }
    
    // Calculate dbar used to calculate CF, also calculate ai values
    for(int xyz = 1; xyz <= ntrt; xyz++)
    {
      if(polyKTestType == POLY_3_TEST)
      {
        ai[xyz] = (subsetVars.getPoly3(xyz) * subsetVars.getPoly3(xyz)) / subsetVars.getN_subj(xyz);
      }
      else if(polyKTestType == POLY_1PT5_TEST)
      {
        ai[xyz] = (subsetVars.getPoly15(xyz) * subsetVars.getPoly15(xyz)) / subsetVars.getN_subj(xyz);
      }
      else if(polyKTestType == POLY_6_TEST)
      {
        ai[xyz] = (subsetVars.getPoly6(xyz) * subsetVars.getPoly6(xyz)) / subsetVars.getN_subj(xyz);
      }
      
      dBartop += ai[xyz] * scale[xyz];
      dBarbot += ai[xyz];
    } // label 5
    dBar = dBartop / dBarbot;
    // Calculate CF value
    if(ntrt == 2) // Note this is mathematically the same as the trend correction if we assume getScale(2) - getScale(1) = 1
    {
      if(polyKTestType == POLY_3_TEST)
      {
        cf = (ai[1] / (ai[1] + ai[2]))
        * (subsetVars.getPoly3(2) / subsetVars.getN_subj(2))
        + (ai[2] / (ai[1] + ai[2]))
        * (subsetVars.getPoly3(1) / subsetVars.getN_subj(1));
      }
      else if(polyKTestType == POLY_1PT5_TEST)
      {
        cf = (ai[1] / (ai[1] + ai[2]))
        * (subsetVars.getPoly15(2) / subsetVars.getN_subj(2))
        + (ai[2] / (ai[1] + ai[2]))
        * (subsetVars.getPoly15(1) / subsetVars.getN_subj(1));
      }
      else if(polyKTestType == POLY_6_TEST)
      {
        cf = (ai[1] / (ai[1] + ai[2]))
        * (subsetVars.getPoly6(2) / subsetVars.getN_subj(2))
        + (ai[2] / (ai[1] + ai[2]))
        * (subsetVars.getPoly6(1) / subsetVars.getN_subj(1));
      }
    }
    else
    {
      cf = 0;
      for(int xyz2 = 2; xyz2 <= ntrt; xyz2++)
      {
        if(polyKTestType == POLY_3_TEST)
        {
          front = (subsetVars.getPoly3(xyz2) / subsetVars.getN_subj(xyz2))
          * (scale[xyz2] - dBar);
          back = (subsetVars.getPoly3(xyz2 - 1) / subsetVars.getN_subj(xyz2 - 1))
            * (scale[xyz2 - 1] - dBar);
        }
        else if(polyKTestType == POLY_1PT5_TEST)
        {
          front = (subsetVars.getPoly15(xyz2) / subsetVars.getN_subj(xyz2))
          * (scale[xyz2] - dBar);
          back = (subsetVars.getPoly15(xyz2 - 1) / subsetVars.getN_subj(xyz2 - 1))
            * (scale[xyz2 - 1] - dBar);
        }
        else if(polyKTestType == POLY_6_TEST)
        {
          front = (subsetVars.getPoly6(xyz2) / subsetVars.getN_subj(xyz2))
          * (scale[xyz2] - dBar);
          back = (subsetVars.getPoly6(xyz2 - 1) / subsetVars.getN_subj(xyz2 - 1))
            * (scale[xyz2 - 1] - dBar);
        }
        ai2[xyz2] = front - back;
        if(ai2[xyz2] > cf)
        {
          cf = ai2[xyz2];
        }
      } // label 6
    }
 
    cf /= 2;
    
    sum_adpi = 0;
    sum_ad = 0;
    sum_api = 0;
    sum_a = 0;
    sum_add = 0;
    sum6 = 0;
    sum7 = 0;
    sum8 = 0;
    sum9 = 0;
    sum10 = 0;
    sum11 = 0;
    y_n = 0;
    y_nk = 0;
    y_ak = 0;
    y_pk = 0;
    y_dk = 0;
    
    // Calculate denominator for S value
    base = (double)subsetVars.getN_subj_tot() - (double)ntrt;
     
    if(base <= 0) return pvalue;
    
    xbn = 0;
    xam = 0;
    
    for(int k = 1; k <= ntrt; k++)
    {
      sumj_bn0 = 0;
      sumj_am3 = 0;
      sumj_am6 = 0;
      
      /* For each removal age, get the number of animals removed with and without a tumor */
      for(int j = 1; j <= subsetVars.getNage(); j++)
      {
        xbn = subsetVars.getTumorAnimals(j, k);  // Number of tumor animals
        xam = subsetVars.getNonTumorAnimals(j, k);  // Number of non-tumor animals
        sumj_bn0 += xbn;
        // Set up for summing
        timeValue = subsetVars.getT(j);
        switch(polyKTestType)
        {
        case POLY_3_TEST:
          exp1 = 3;
          exp2 = 6;
          break;
        case POLY_1PT5_TEST:
          exp1 = 1.5;
          exp2 = 3;
          break;
        case POLY_6_TEST:
          exp1 = 6;
          exp2 = 12;
          break;
        default:
          break;
        }
        // Do summing
        sumj_am3 += xam * pow(timeValue, exp1);
        sumj_am6 += xam * pow(timeValue, exp2);

      } // label 10
      
      y_n = subsetVars.getN_subj(k);
      y_nk = sumj_bn0 + sumj_am3;
      if(y_n <= 0 || y_nk <= 0) return pvalue;
      
      y_ak = y_nk * y_nk / y_n;
      y_pk = sumj_bn0 / y_nk;
      y_dk = scale[k];
      
      sum_adpi += y_ak * y_pk * y_dk;
      sum_ad += y_ak * y_dk;
      sum_api += y_ak * y_pk;
      sum_a += y_ak;
      sum_add += y_ak * y_dk * y_dk;
      sum6 += sumj_bn0;
      sum7 += sumj_am3;
      sum8 += sumj_am6;
      sum9 += sumj_bn0 * sumj_bn0 / y_n;
      sum10 += sumj_am3 * sumj_am3 / y_n;
      sum11 += sumj_bn0 * sumj_am3 / y_n;
      
    } // label 20
    
    if(sum_a <= 0) return pvalue;
    
    top = sum_adpi - (sum_ad * sum_api / sum_a); // This is A in the big equation
    den_p = sum6 + sum7;
    if(den_p <= 0) return pvalue;
    
    p = sum6 / den_p;
    q = sum7 / den_p;
    term1 = q * q * sum6 + p * p * sum8;
    term2 = q * q * sum9 + p * p * sum10 - 2 * p * q * sum11;
    term3 = (term1 - term2) / base;
    
    bot = term3 * (sum_add - sum_ad * sum_ad / sum_a); // This is B in the big equation
    
    // If the bot value is <= 0 we default to the Cochran-Armitage (or Fisher's Exact) p-value
    if(bot <= POLYK_EPSILON) return cochArmFishPValue;
    
    topcc = fabs(top) - cf;
    
    test_stat = topcc * topcc / bot;
    ndf = 1;
    
   
    pvalue = chisq(test_stat, ndf);
    
    pvalue /= 2;
    if(topcc < 0) pvalue = 1 - pvalue;
    if(pvalue < MIN_PVALUE) pvalue = MIN_PVALUE;
    if(fabs(test_stat) < TESTSTAT_EPSILON) pvalue = 0.5;
    if(top < 0) pvalue = -1 * pvalue;
  }
  catch (std::exception &e)
  {
  //  System.out.println("Exception in polyk_mod = " + e);
    pvalue = BAD_PVALUE;
  }
  return pvalue;
}

/**
 * Calculate the chisq probability for the given test statistic and degrees of freedom
 */
double TDMSE_PolyK::chisq(double x2, int idf)
{
  double p;
  double d2;
  double x22;
  
  if(fabs(x2) < TESTSTAT_EPSILON && idf == 1) return 0.5;
  d2 = idf / 2.0;
  x22 = x2 / 2.0;
  p = 1.0 - gammp(d2, x22);
  
  // This code is an alternative to the above using apache library
  //        DistributionFactory factory = DistributionFactory.newInstance();
  //        ChiSquaredDistribution chiSq = factory.createChiSquareDistribution(idf);
  //        p = 1 - chiSq.cummulativeProbability(x2);
  
  return p;
}

/**
 * Helper function for chisq calculation
 */
double TDMSE_PolyK::gammp(double a, double x)
{
  double returnValue;
  if(x < a + 1.0)
  {
    returnValue = gser(a, x);
  }
  else
  {
    double gammcf = gcf(a, x);
    returnValue = 1.0 - gammcf;
  }
  

  return returnValue;
}

/**
 * Helper function for chisq calculation
 */
double TDMSE_PolyK::gser(double a, double x)
{
  double gamser;
  
  double gln = gammln(a);
  double ap = a;
  double sum = 1.0 / a;
  double del = sum;
  
  for(int n = 1; n <= ITMAX; n++)
  {
    ap += 1.0;
    del *= x / ap;
    sum += del;
    if(fabs(del) < (fabs(sum) * EPS)) break;
  }
  
  if(x <= 0) x = 0.01;
  
  gamser = sum * exp(-x + a * log(x) - gln);
  return gamser;
}

/**
 * Helper function for chisq calculation
 */
double TDMSE_PolyK::gcf(double a, double x)
{
  double gammcf;
  double ana;
  double anf;
  
  double g = 0.0;
  double gln = gammln(a);
  double gold = 0.0;
  double a0 = 1.0;
  double a1 = x;
  double b0 = 0.0;
  double b1 = 1.0;
  double fac = 1.0;
  
  for(int n = 1; n <= ITMAX; n++)
  {
    ana = (double)n - a;
    a0 = (a1 + a0 * ana) * fac;
    b0 = (b1 + b0 * ana) * fac;
    anf = (double)n * fac;
    a1 = x * a0 + anf * a1;
    b1 = x * b0 + anf * b1;
    if(fabs(a1) > std::numeric_limits<double>::lowest())
    {
      fac = 1.0 / a1;
      g = b1 * fac;
      if(fabs((g - gold) / g) < EPS) break;
      gold = g;
    }
  }
  
  gammcf = exp(-x + a * log(x) - gln) * g;
  
  return gammcf;
}

/**
 * Helper function for chisq calculation
 **/
double TDMSE_PolyK::gammln(double xx)
{
  double returnValue;
  
  double x = xx - 1.0;
  double tmp = x + FPF;
  tmp = (x + 0.5) * log(tmp) - tmp;
  double ser = 1.0;
  
  for(int j = 0; j < NUM_COF; j++)
  {
    x += 1.0;
    ser += COF[j] / x;
  }
  
  returnValue = tmp + log(STP * ser);
  
  return returnValue;
}

}
  
