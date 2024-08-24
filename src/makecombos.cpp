#include <Rcpp.h>
using namespace Rcpp;

// Helper function to call R's findInterval function from the "stats" package
NumericVector callFindInterval(NumericVector x, NumericVector vec) {
  Function findInterval("findInterval", "stats");
  return findInterval(x, vec);
}

// [[Rcpp::export]]
DataFrame model_hichip(NumericVector x, Nullable<NumericVector> y = R_NilValue, NumericVector borders = NumericVector::create(), bool yvals = true) {
  NumericVector y_vec;
  
  if (y.isNotNull()) {
    y_vec = as<NumericVector>(y);
  } else {
    y_vec = NumericVector(1, NA_REAL);
  }

  // Input validation: Check if x and y have the same length if yvals is true
  if (yvals && x.size() != y_vec.size()) {
    stop("x and y must have the same length when yvals is true");
  }
  
  int n = x.size();
  int nBins = borders.size() + 1;
  
  // Find bin indices for x values based on borders
  NumericVector bin = callFindInterval(x, borders);
  NumericVector sumofy(nBins, 0.0);
  NumericVector meanofx(nBins, NA_REAL);  // Initialize with NA_REAL
  NumericVector sumofx(nBins, 0.0);
  NumericVector countofx(nBins, 0.0);
  
  // Calculate sums and counts
  for (int i = 0; i < n; i++) {
    int idx = bin[i];
    if (idx >= 0 && idx < nBins) {  // Ensure bin index is within bounds
      if (yvals) {
        sumofy[idx] += y_vec[i];
      } else {
        sumofy[idx] += 1.0;
      }
      sumofx[idx] += x[i];
      countofx[idx] += 1.0;
    }
  }
  
  // Calculate mean values
  for (int i = 0; i < nBins; i++) {
    if (countofx[i] > 0) {
      meanofx[i] = sumofx[i] / countofx[i];
    }
  }
  
  // Calculate probabilities based on sumofy
  double sumofsumofy = sum(sumofy);
  NumericVector pvals = sumofy / sumofsumofy;
  
  // Create and return DataFrame with results
  return DataFrame::create(
    _["meanofx"] = meanofx,
    _["sumofy"] = sumofy,
    _["pvals"] = pvals,
    _["sumofx"] = sumofx,
    _["countofx"] = countofx
  );
}

// [[Rcpp::export]]
DataFrame makecombos(CharacterVector chrom, DataFrame chrpeaks, double mindist = 0, double maxdist = 100000000) {
  // Input validation: Ensure chrpeaks contains required columns
  if (!chrpeaks.containsElementNamed("chr") || !chrpeaks.containsElementNamed("start") ||
      !chrpeaks.containsElementNamed("end") || !chrpeaks.containsElementNamed("score")) {
    stop("chrpeaks must contain columns: chr, start, end, and score");
  }
  
  StringVector chrpeaks_chr = chrpeaks["chr"];
  IntegerVector chrpeaks_start = chrpeaks["start"];
  IntegerVector chrpeaks_end = chrpeaks["end"];
  NumericVector chrpeaks_score = chrpeaks["score"];
  
  int npeaks = chrpeaks.nrows();
  int estimatedSize = npeaks * (npeaks - 1) / 2; // Upper bound of possible combinations
  
  // Pre-allocate vectors based on the estimated size
  NumericVector score1(estimatedSize);
  NumericVector score2(estimatedSize);
  IntegerVector dist(estimatedSize);
  
  int validCombos = 0;
  
  // Calculate distances and filter combinations based on the distance criteria
  for (int i = 0; i < npeaks; ++i) {
    for (int j = i + 1; j < npeaks; ++j) {
      int peakDist = std::abs((chrpeaks_start[i] + chrpeaks_end[i]) / 2 - 
                              (chrpeaks_start[j] + chrpeaks_end[j]) / 2);
      
      if (peakDist >= mindist && peakDist <= maxdist) {
        score1[validCombos] = chrpeaks_score[i];
        score2[validCombos] = chrpeaks_score[j];
        dist[validCombos] = peakDist;
        validCombos++;
      }
    }
  }
  
  Rcout << "Valid combination number: " << validCombos << std::endl;
  
  // Trim vectors to the actual number of valid combinations and create DataFrame
  return DataFrame::create(
    Named("score1") = score1[Range(0, validCombos - 1)],
    Named("score2") = score2[Range(0, validCombos - 1)],
    Named("dist") = dist[Range(0, validCombos - 1)]
  );
}
