#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector countsToProbs(NumericVector x){
  // purpose : Takes a vector of counts and returns a vector of the same
  //           length where each entry is the proportion of the total counts
  //           of that entry in x
  
  int xLen = x.size();
  
  int sum = 0;
  for (int i=0; i < xLen; i++){sum = sum + x[i];}
  
  NumericVector output(xLen);
  for (int i=0; i < xLen; i++){output[i] = x[i]/sum;}
  return output;
}

// [[Rcpp::export]]
NumericVector vectorToCounts(NumericVector x, NumericVector breaks){
  // purpose : Takes a vector of values, and a vector of breaks and returns
  //           a vector which indicates how many of the observations are in 
  //           each of the intervals defined by breaks.
  // inputs  : x      - the vector of values which we want to count through
  //           breaks - an ordered (increasing) vector of break points, defining
  //                    the boundaries of the length(breaks) - 1 intervals.
  // output  :
  // initialise output vector:
  int intervalNum;
  int xLen;
  intervalNum = breaks.size() - 1;
  xLen = x.size();
  NumericVector output(intervalNum);
  
  for (int i=0; i<xLen; i++){
    // for each element in x:
    for (int j=0; j<intervalNum; j++){
      // for each break:
      if ((x[i]>breaks[j]) & (x[i]<=breaks[j+1])){
        output[j]++;
      }
    }
  }
  
  return output;
}//vectorToCounts

// [[Rcpp::export]]
NumericVector weightedSelection(NumericVector x, NumericVector probs){
  // purpose : Takes a vector of values, and a vector of probs and returns
  //           a vector which contains each element of x with the probability
  //           specific by probs
  // inputs  : x      - the vector of values which we want to select from
  //           breaks - an ordered vector where the ith entry is the probability
  //                    of keeping entry i of x.
  // output  : A numeric vector of selected values from x
  
  // declare vars:
  int xLen;
  xLen = x.size();
  NumericVector storage(xLen);
  NumericVector deviate = runif(xLen);
  int counter;
  counter = 0;
  
  for (int i=0; i<xLen; i++){
    if (deviate[i]<probs[i]){
      storage[counter] = x[i];
      counter++;
    }
  }
  
  NumericVector output(counter);
  for (int i=0; i<counter; i++){
    output[i] = storage[i];
  }
  
  return output;
}//weightedSelection


// [[Rcpp::export]]
NumericVector divv(NumericVector x){
  // purpose : Inspired by the R function 'diff'. For each entry in position i
  //           returns x[i+1]/x[i]
  // inputs  : x      - the vector of values which we want to iterate over
  // output  : A NumericVector of ratios
  
  int xLen;
  xLen = x.size();
  
  NumericVector output(xLen-1);
  
  for (int i=1; i<xLen; i++){
    output[i-1] = x[i]/x[i-1];
  }
  
  return output;
}//divv

NumericVector fetchFromLogical(NumericVector x, LogicalVector I){
  // purpose : returns a NumericVector with elements of x corresponding to when
  //            elements of I are true 
  int xLen;
  xLen = x.size();
  NumericVector output;
  for (int i=0; i<xLen; i++){
    if (I[i]){
      output.push_back(x[i]);
    }
  }
  return output;
}//fetchFromLogical

NumericVector fetchFromIndex(NumericVector x, int lower, int upper){
  // purpose : returns elements lower:upper of x
  NumericVector output;
  for (int i=lower; i<upper; i++){output.push_back(x[i]);}
  return output;
}//fetchFromIndex

NumericVector repCpp(NumericVector x, NumericVector y){
  // purpose : like the rep function in R
  NumericVector output;
  int xLen = x.size();
  
  for (int i=0; i < xLen; i++){
    int jMax = y[i];
    for (int j=0; j<jMax; j++){
      output.push_back(x[i]);
    }
  }
  
  return output;
}//repCpp

NumericVector cRcpp(NumericVector x, NumericVector y){
  // purpose : like the c function in R
  NumericVector output;
  int xLen = x.size();
  int yLen = y.size();
  for (int i=0; i < xLen; i++){output.push_back(x[i]);}
  for (int i=0; i < yLen; i++){output.push_back(y[i]);}
  return output;
}//repCpp

// 
// //[[Rcpp::export]]
// NumericVector sampleStateIPMrcpp(NumericVector previousState,
//                                  Function survFunc,
//                                  NumericVector survPars,
//                                  Function growthSamp,
//                                  NumericVector growthPars,
//                                  Function reprFunc,
//                                  NumericVector reprPars,
//                                  Function offNumSamp,
//                                  NumericVector offNumPars,
//                                  Function offSizeSamp,
//                                  NumericVector offSizePars,
//                                  double Schild,
//                                  NumericVector breaks,
//                                  LogicalVector oneSex,
//                                  double shift){
//   // purpose : the sampleStateIPM function in Cpp code
//   // inputs  : SEE sampleStateIPM
//   // output  : A NumericVector of counts of individuals in each size class at
//   //           the next time step
// 
//   Schild = exp(Schild)/(1+exp(Schild));
//   shift = exp(shift)/(1+exp(shift));
// 
//   if (oneSex){
//     Schild = Schild/2;
//   }
// 
//   int D = previousState.size();
// 
//   NumericVector sizes = fetchFromIndex(breaks, 0, D);
// 
//   for (int i=1; i<=D; i++){
//     sizes[i-1] = shift*(breaks[i] - breaks[i-1]) + breaks[i-1];
//   }
// 
//   NumericVector survProbs = survFunc(sizes, survPars);
//   NumericVector survSelectFrom = repCpp(sizes, previousState);
//   NumericVector allSurvProbs = repCpp(survProbs, previousState);
//   NumericVector survSizes = weightedSelection(survSelectFrom, allSurvProbs);
//   NumericVector grownSizes = growthSamp(survSizes, growthPars);
//   NumericVector reprProbs = reprFunc(grownSizes, reprPars);
//   NumericVector reprSizes = weightedSelection(grownSizes, reprProbs);
//   NumericVector offNumbers = offNumSamp(reprSizes, offNumPars);
//   NumericVector parentSizes = repCpp(reprSizes, offNumbers);
//   NumericVector bornSizes = offSizeSamp(parentSizes, offSizePars);
//   int childNum = bornSizes.size();
//   NumericVector vecSchild = repCpp(Schild, childNum);
//   NumericVector childSizes = weightedSelection(bornSizes, vecSchild);
//   NumericVector combinedSizes = cRcpp(childSizes, grownSizes);
//   NumericVector output = vectorToCounts(combinedSizes, breaks);
//   return output;
// }//sampleStateIPMrcpp
