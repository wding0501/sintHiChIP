#include <Rcpp.h>
#include <cstdint>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// Helper function to call R's findInterval function from the "stats" package
NumericVector callFindInterval(NumericVector x, NumericVector vec) {
  Function findInterval("findInterval", "stats");
  return findInterval(x, vec);
}

// [[Rcpp::export]]
DataFrame model_hichip(NumericVector x, Nullable<NumericVector> y = R_NilValue, 
                       NumericVector borders = NumericVector::create(), 
                       bool yvals = true) {
  // Convert nullable y to vector
  NumericVector y_vec;
  if (y.isNotNull()) {
    y_vec = as<NumericVector>(y);
  } else {
    y_vec = NumericVector(1, NA_REAL);
  }
  
  // Input validation
  if (yvals && y.isNotNull() && x.size() != y_vec.size()) {
    stop("x and y must have the same length when yvals is true");
  }
  
  // Initialize variables
  int n = x.size();
  int nBins = borders.size() + 1;
  
  // Find bin indices for x values based on borders
  NumericVector bin = callFindInterval(x, borders);
  NumericVector sumofy(nBins, 0.0);
  NumericVector meanofx(nBins, NA_REAL);
  NumericVector sumofx(nBins, 0.0);
  NumericVector countofx(nBins, 0.0);
  
  // Calculate sums and counts for each bin
  for (int i = 0; i < n; i++) {
    int idx = static_cast<int>(bin[i]); // Ensure proper type conversion
    if (idx >= 0 && idx < nBins) {
      if (yvals && y.isNotNull()) {
        sumofy[idx] += y_vec[i];
      } else {
        sumofy[idx] += 1.0;
      }
      sumofx[idx] += x[i];
      countofx[idx] += 1.0;
    }
  }
  
  // Calculate mean values for each bin
  for (int i = 0; i < nBins; i++) {
    if (countofx[i] > 0) {
      meanofx[i] = sumofx[i] / countofx[i];
    }
  }
  
  // Calculate probabilities
  double sumofsumofy = sum(sumofy);
  NumericVector pvals(nBins);
  if (sumofsumofy > 0) {
    pvals = sumofy / sumofsumofy;
  } else {
    std::fill(pvals.begin(), pvals.end(), 0.0);
  }
  
  // Return results as DataFrame
  return DataFrame::create(
    _["meanofx"] = meanofx,
    _["sumofy"] = sumofy,
    _["pvals"] = pvals,
    _["sumofx"] = sumofx,
    _["countofx"] = countofx
  );
}


// [[Rcpp::export]]
DataFrame makecombos(CharacterVector chrom, DataFrame chrpeaks, 
                     int64_t mindist = 0, int64_t maxdist = 100000000) {
  
  // === Input validation ===
  if (!chrpeaks.containsElementNamed("chr") || 
      !chrpeaks.containsElementNamed("start") ||
      !chrpeaks.containsElementNamed("end") || 
      !chrpeaks.containsElementNamed("score")) {
    stop("chrpeaks must contain columns: chr, start, end, and score");
  }
  
  // === Extract data and convert to C++ vectors ===
  StringVector chrpeaks_chr = chrpeaks["chr"];
  IntegerVector chrpeaks_start = chrpeaks["start"];
  IntegerVector chrpeaks_end = chrpeaks["end"];
  NumericVector chrpeaks_score = chrpeaks["score"];
  
  int npeaks = chrpeaks.nrows();
  
  if (npeaks < 2) {
    // Rcout << "Valid combination number: 0" << std::endl;
    return DataFrame::create(
      Named("score1") = NumericVector(0),
      Named("score2") = NumericVector(0),
      Named("dist") = IntegerVector(0)
    );
  }
  
  // === Precompute midpoints and sorting ===
  std::vector<int> midpoints(npeaks);
  std::vector<double> scores(npeaks);
  std::vector<int> indices(npeaks);
  
  // Calculate midpoints and build indices
  for (int i = 0; i < npeaks; ++i) {
    midpoints[i] = (chrpeaks_start[i] + chrpeaks_end[i]) / 2;
    scores[i] = chrpeaks_score[i];
    indices[i] = i;
  }
  
  // Sort indices by midpoint
  std::sort(indices.begin(), indices.end(), 
            [&midpoints](int a, int b) {
              return midpoints[a] < midpoints[b];
            });
  
  // === Memory allocation strategy ===
  long long peak_span = static_cast<long long>(midpoints[indices[npeaks-1]]) - 
                        static_cast<long long>(midpoints[indices[0]]);
  double avg_density = static_cast<double>(npeaks) / (peak_span + 1);
  
  long long avg_partners = std::min(static_cast<long long>(avg_density * maxdist), 
                                    static_cast<long long>(npeaks));
  long long estimated_combos = std::min(static_cast<long long>(npeaks) * avg_partners / 4, 
                                        static_cast<long long>(npeaks) * (npeaks - 1) / 2);
  
  // Limit initial capacity to avoid memory issues
  int initial_capacity = std::min(estimated_combos, 50000000LL);
  
  std::vector<double> score1_vec;
  std::vector<double> score2_vec;
  std::vector<int> dist_vec;
  
  score1_vec.reserve(initial_capacity);
  score2_vec.reserve(initial_capacity);
  dist_vec.reserve(initial_capacity);
  
  // === Main loop with binary search ===
  long long total_valid = 0;
  
  for (int i = 0; i < npeaks - 1; ++i) {
    int idx_i = indices[i];
    int mid_i = midpoints[idx_i];
    double score_i = scores[idx_i];
    
    // Binary search for valid range start
    int64_t min_mid = static_cast<int64_t>(mid_i) + mindist;  
    int j_start = i + 1;
    
    if (mindist > 0) {
      int left = i + 1, right = npeaks;
      while (left < right) {
        int mid = (left + right) / 2;
        if (static_cast<int64_t>(midpoints[indices[mid]]) < min_mid) {  
          left = mid + 1;
        } else {
          right = mid;
        }
      }
      j_start = left;
    }
    
    // Binary search for valid range end
    int64_t max_mid = static_cast<int64_t>(mid_i) + maxdist;  
    int j_end = npeaks;
    
    int left = j_start, right = npeaks;
    while (left < right) {
      int mid = (left + right) / 2;
      if (static_cast<int64_t>(midpoints[indices[mid]]) <= max_mid) {  
        left = mid + 1;
      } else {
        right = mid;
      }
    }
    j_end = left;
    
    // Check if vector expansion is needed
    int batch_size = j_end - j_start;
    if (batch_size <= 0) continue;
    
    long long new_total = total_valid + batch_size;
    if (new_total > static_cast<long long>(score1_vec.capacity())) {
      // Expand capacity with upper limit protection
      size_t new_capacity = std::min(score1_vec.capacity() * 2, 
                                     score1_vec.capacity() + 20000000);
      score1_vec.reserve(new_capacity);
      score2_vec.reserve(new_capacity);
      dist_vec.reserve(new_capacity);
    }
    
    // Add valid combinations
    for (int j = j_start; j < j_end; ++j) {
      int idx_j = indices[j];
      int dist = midpoints[idx_j] - mid_i;
      
      score1_vec.push_back(score_i);
      score2_vec.push_back(scores[idx_j]);
      dist_vec.push_back(dist);
    }
    
    total_valid = new_total;
    
    // Progress reporting and user interrupt checking
    if (npeaks > 10000 && i % (npeaks / 10) == 0) {
      Rcpp::checkUserInterrupt();
    }
  }
  
  // Rcout << "Valid combination number: " << total_valid << std::endl;
  
  if (total_valid == 0) {
    return DataFrame::create(
      Named("score1") = NumericVector(0),
      Named("score2") = NumericVector(0),
      Named("dist") = IntegerVector(0)
    );
  }
  
  // Convert to R vectors
  NumericVector r_score1(total_valid);
  NumericVector r_score2(total_valid);
  IntegerVector r_dist(total_valid);
  
  std::copy(score1_vec.begin(), score1_vec.end(), r_score1.begin());
  std::copy(score2_vec.begin(), score2_vec.end(), r_score2.begin());
  std::copy(dist_vec.begin(), dist_vec.end(), r_dist.begin());
  
  return DataFrame::create(
    Named("score1") = r_score1,
    Named("score2") = r_score2,
    Named("dist") = r_dist
  );
}