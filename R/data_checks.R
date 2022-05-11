#Copyright 2021  NIEHS <matt.wheeler@nih.gov>
#   
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
#and associated documentation files (the "Software"), to deal in the Software without restriction, 
#including without limitation the rights to use, copy, modify, merge, publish, distribute, 
#sublicense, and/or sell copies of the Software, and to permit persons to whom the Software 
#is furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies 
#or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
#CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
#OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
.check_for_na <- function(data){
 

  is_a_na <- rep(FALSE,nrow(data)) 
  for (ii in 1:ncol(data)){
    is_a_na <- is_a_na + is.na(data[,ii]) + is.infinite(data[,ii]) + is.nan(data[,ii])
  }
  
  if ( sum(is_a_na) > 0){
    warning("Infinite or NA value was found.  These data row(s)
             were removed from the analysis.")
  }
  
  data <- data[!(is_a_na > 0),,drop=F]

  if (nrow(data) < 3){
    stop("Less than three viable data rows were found.")
  }
  return(!(is_a_na > 0))
}
 
 .check_negative_response <- function(Y){

   if (sum(Y[,1] <=0) >=1  ){
      return(TRUE)
   }else{
      return(FALSE)
   }
 }