

//necessary things to run in R    
#include <RcppEigen.h>

#include <vector>
#include <array>
#include <list>
#include <algorithm>
#include "polyK_setup.h"
#include "polyK_animal.h"


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
namespace PolyK{


bool PolyKPrepareClass::SetupStudy(std::vector<double> dose, std::vector<int> tumor, 
                                   std::vector<int> daysOnStudy) {
  Animal an;

               
   try{
     if (dose.size() != tumor.size() ||
         dose.size() != daysOnStudy.size()){
       throw std::runtime_error("Error With PolyKPrepare Init"); 
     }
     for(int i = 0; i < dose.size(); i++){
       an.set(dose[i], tumor[i], daysOnStudy[i]);
       m_AnimalList.push_back(an);
     }
   }catch (std::exception ex) {
     return false;
   }
   return true;
} 


void PolyKPrepareClass::prepare() {
  
  // These arrays are indexed by dose level
  std::vector<double> tumorSum(MAX_DOSE_LEVELS);//* = new double[MAX_DOSE_LEVELS];
  
  // initialize output variables
  dmax = -1;
  
  // Assume animals are in dose order initially
  
  double prevDose = -1.0;
  double maxDose = -1.0;
  int doseIdx = 0;
  m_MaxTime = -1;
  
  // Loop through the animals to get N values for each dose group
  prevDose = -1;
  
  for (auto &an : m_AnimalList){ 
    
    if (an.getDose() > prevDose) {
      doseIdx++;
      m_Scale[doseIdx] = an.getDose();
      m_NumAnimals[doseIdx] = 0;
    }
    an.setDoseIdx(doseIdx);
    if (an.getDose() > maxDose){
      maxDose = an.getDose();
    }
    if (an.getDaysOnStudy() > m_MaxTime){
      m_MaxTime = an.getDaysOnStudy();
    }
    
    m_NumAnimals[doseIdx] = m_NumAnimals[doseIdx] + 1;
    prevDose = an.getDose();
  }
  
  // Loop through the animals to get poly-k values for each dose group
  prevDose = -1;
  doseIdx = 0;
  for (auto &an : m_AnimalList){ 
    
    if (an.getDose() > prevDose){
      doseIdx++;
    }
    
    if (an.getTumor() == 1) {
      tumorSum[doseIdx] = tumorSum[doseIdx] + 1;
      m_Poly15Denom[doseIdx] = m_Poly15Denom[doseIdx] + 1;
      m_Poly3Denom[doseIdx] = m_Poly3Denom[doseIdx] + 1;
      m_Poly6Denom[doseIdx] = m_Poly6Denom[doseIdx] + 1;
    } else {
      m_Poly15Denom[doseIdx] = m_Poly15Denom[doseIdx] + pow(an.getDaysOnStudy() / (double)m_MaxTime, 1.5);
      m_Poly3Denom[doseIdx]  = m_Poly3Denom[doseIdx]  + pow(an.getDaysOnStudy() / (double)m_MaxTime, 3);
      m_Poly6Denom[doseIdx]  = m_Poly6Denom[doseIdx]  + pow(an.getDaysOnStudy() / (double)m_MaxTime, 6);
    }
    prevDose = an.getDose();
    
  }
  
  m_NumDoseLevels = doseIdx;
  
  // Calculate poly-k incidence rates for each dose
  for (int i = 1; i <= m_NumDoseLevels; i++) {
    // Scale the dose values to a 0-1 range
    m_Scale[i] = m_Scale[i] / maxDose;
  //  Rcpp::Rcout << "Poly-3 Denom # " <<  i << " = " << m_Poly3Denom[i] <<std::endl;
  //  Rcpp::Rcout << "TumorSum #" << i << " = " << tumorSum[i] << std::endl;
   // Rcpp::Rcout << "Scale #" << i << " = " << m_Scale[i] << std::endl;
  }
  
  // Sort animals by age
  m_AnimalList.sort(); 
  // Collections.sort(m_AnimalList);
  
  // Loop through animals to fill in ages
  int prevAge = -1;
  int ageIdx = 0;
 
  for (auto &an : m_AnimalList){ 
   
    if (an.getDaysOnStudy() > prevAge) {
      ageIdx++;
      m_TimeArray[ageIdx] = an.getDaysOnStudy() / (double)m_MaxTime;
    }
    prevAge = an.getDaysOnStudy();
    doseIdx = an.getDoseIdx();
    if (an.getTumor() == 1) {
      m_TumorAnimals(ageIdx,doseIdx) = m_TumorAnimals(ageIdx,doseIdx) + 1;
    } else {
      m_NonTumorAnimals(ageIdx,doseIdx) = m_NonTumorAnimals(ageIdx,doseIdx) + 1;
    }
  }
  m_NAge = ageIdx;


}
}


