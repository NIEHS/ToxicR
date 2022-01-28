#include <vector>
#include <array>
#include <list>
#include <algorithm>
#include <RcppEigen.h>
#include <RcppGSL.h>

#include "polyK_animal.h"

#ifndef __POLYK_SETUP_H
#define __POLYK_SETUP_H

namespace PolyK{


/**
 * Helper class used by P8_Calculator to calculate dose and age vectors for each sex
 * 
 * @author $Author: mercerm $
 * @version $Revision: 9150 $ $Date: 2011-05-16 14:08:23 -0400 (Mon, 16 May
 *          2011) $
 * @since 1.0
 */
static  int MAX_DOSE_LEVELS = 20;
static  int MAX_AGES = 4000;

class PolyKPrepareClass {
  /**
   * Helper class to hold information for one animal
   * @author harriss2
   *
   */

private: 
  // Constants
  
  // Values passed from P8_Calculator
  int m_MaxTime;
  int m_NAge;
  Eigen::VectorXi m_NumAnimals; 
 // int  m_NumAnimals*;
  Eigen::VectorXd m_TimeArray;
  // double m_TimeArray*;
  Eigen::VectorXd m_Poly15Denom; 
  //double m_Poly15Denom*;
  Eigen::VectorXd m_Poly3Denom; 
  //double m_Poly3Denom*;
  Eigen::VectorXd m_Poly6Denom; 
  //double m_Poly6Denom*;
  Eigen::VectorXd m_Scale; 
  // double m_Scale*;
  Eigen::MatrixXi m_TumorAnimals; 
 // int  m_TumorAnimals*;
  Eigen::MatrixXi m_NonTumorAnimals; 
 // int  m_NonTumorAnimals*;

  std::list<Animal> m_AnimalList;
  

  
  int m_NumDoseLevels;
  
  double dosediffmax;
  double tmax;
  double dmax; // Max dose value
  public:
  
  /**
   * Create a new instance of PolyKPrepareClass
   */
  PolyKPrepareClass()
  {
	  Eigen::VectorXi temp_im(MAX_DOSE_LEVELS + 1);
	  Eigen::VectorXd temp_ti(MAX_AGES+1); 
	  Eigen::VectorXd temp_m(MAX_DOSE_LEVELS+ 1);   
	  Eigen::MatrixXi temp_i(MAX_AGES+1,MAX_DOSE_LEVELS+ 1);        

	  m_Scale = temp_m*0;
	  m_NumAnimals = temp_im*0; 
	  m_TimeArray  = temp_ti*0;
	  m_Poly15Denom = temp_m*0;
	  m_Poly3Denom = temp_m*0;  
	  m_Poly6Denom = temp_m*0; 
	  m_TumorAnimals = temp_i*0;
	  m_NonTumorAnimals = temp_i*0; 
 
  }
  
  ~PolyKPrepareClass() {
    
    // Initialize arrays
  
    
  }


  /**
   * Get the weighting value (scaled dose value) for the treatment group index
   * @param idx Treatment group index
   * @return Scaled dose level to apply to treatment group
   */
  double getWeight(int idx) {
    return m_Scale[idx];
  };
  
  /**
   * Get the number of distinct ages
   * @return Number of distinct ages
   */
  int getNage() {
    return m_NAge;
  };
  
  /**
   * Get the number of distinct dose levels
   * @return Number of treatment groups
   */
  int getNumDoseLevels() {
    return m_NumDoseLevels;
  };
  
  /**
   * Get the scaled time on dose (days on study/study length)
   * @param idx The age index (1-based)
   * @return Days on dose
   */
  double getT(int idx) {
    return m_TimeArray[idx];
  };
  
  /**
   * Get the number of animals in a dose group
   * @param idx Dose group index (1-based)
   * @return Number of animals in the dose group
   */
  int getN_subj(int idx) {
    return m_NumAnimals[idx];
  };
  
  /**
   * Get the total number of animals across the dose groups
   * @return Total animal count
   */
  int getN_subj_tot() {
    return m_AnimalList.size();
  };
  
  /**
   * Get the poly-1.5 incidence rate for the specified dose
   * @param idx Dose index (1-based)
   * @return Poly-1.5 adjusted tumor rate
   */
  double getPoly15(int idx) {
    return m_Poly15Denom[idx];
  };
  
  /**
   * Get the poly-3 incidence rate for the specified dose
   * @param idx Dose index (1-based)
   * @return Poly-3 adjusted tumor rate
   */
  double getPoly3(int idx) {
    return m_Poly3Denom[idx];
  };
  
  /**
   * Get the poly-6 incidence rate for the specified dose
   * @param idx Dose index (1-based)
   * @return Poly-6 adjusted tumor rate
   */
  double getPoly6(int idx) {
    return m_Poly6Denom[idx];
  };
  
  /**
   * Get the number of animals with tumor for the given time and dose index
   * @param ageIdx The age index (1-based)
   * @param doseIdx Dose index (1-based)
   * @return Number of tumor animals
   */
  int getTumorAnimals(int ageIdx, int doseIdx)  {
    return m_TumorAnimals(ageIdx,doseIdx);
  };
  
  /**
   * Get the number of animals without tumor for the given time and dose index
   * @param ageIdx The age index (1-based)
   * @param doseIdx Dose index (1-based)
   * @return Number of non-tumor animals
   */
  int getNonTumorAnimals(int ageIdx, int doseIdx)  {
    return m_NonTumorAnimals(ageIdx,doseIdx);
  };
  
  /**
   * Read the text input file
   * @return TRUE if successful
   */
public:
  bool SetupStudy(std::vector<double> dose, std::vector<int> tumor, 
                  std::vector<int> daysOnStudy); 
  
  /**
   * Do intermediate calculations needed for P8_Calculator. compute the dose
   * and age vectors for a given sex
   */
   void prepare();
};


}
#endif
