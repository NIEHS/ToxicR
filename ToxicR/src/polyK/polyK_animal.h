
#ifndef __POLYK_ANIMAL_H
#define __POLYK_ANIMAL_H

class Animal{
private: 
  double dose;
  int doseIdx;
  int tumor;
  int daysOnStudy;
public:
  Animal(){
    
  };
  
  Animal(double doseVal, int tumorVal, int daysVal) {
    dose = doseVal;
    tumor = tumorVal;
    daysOnStudy = daysVal;
    
  };
  
  double getDose() {
    return dose;
  };
  
  int getTumor() {
    return tumor;
  };
  
  int getDaysOnStudy() {
    return daysOnStudy;
  };
  
  void setDoseIdx(int doseIndex) {
    doseIdx = doseIndex;
  };
  
  int getDoseIdx() {
    return doseIdx;
  };
  
  // The animals are sorted by age
  /**
   * Compare this object to another Animal object
   * @param o Other Animal object for comparison
   * @return 0 if equal
   */
public: 
  void set(double dosed, int tumort, int daysOnStudyt){
    dose= dosed; 
    tumor= tumort; 
    daysOnStudy= daysOnStudyt;
  };

     
  bool  operator>(const Animal& b){
    return this->daysOnStudy > b.daysOnStudy; 
  };
  bool  operator>=(const Animal& b){
    return this->daysOnStudy >= b.daysOnStudy; 
  };
  bool  operator<(const Animal& b) const{
    return this->daysOnStudy < b.daysOnStudy; 
  };
  bool operator()(const Animal& A, const Animal& B) const
  {
       return A.daysOnStudy < B.daysOnStudy; 
       //Compare the 2 locations, return true if loc1 is less than loc2
  }
  bool  operator<=(const Animal& b){
    return this->daysOnStudy <= b.daysOnStudy; 
 };
  bool  operator==(const Animal& b){
    return this->daysOnStudy == b.daysOnStudy; 
  };
  
  void operator = (const Animal &A ) { 
    this->dose        = A.dose;
    this->daysOnStudy = A.daysOnStudy;
    this->doseIdx     = A.doseIdx; 
    this->tumor       = A.tumor; 
  };
 
};
#endif
