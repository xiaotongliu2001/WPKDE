#' Find peaks using the estimated values
#'
#' @param k Output of the `kdeC` function, containing estimated values.
#' @param filter A numeric value used to filter out results with estimated values less than the given `filter` argument.
#' @param select A numeric value specifying the number of peaks to retain, selecting the K peaks with the largest estimated values.
#'
#' @return A three-column matrix (`markMat`) where:
#'   - Column 1: x-coordinates of the peaks
#'   - Column 2: y-coordinates of the peaks
#'   - Column 3: Corresponding estimated values of the peaks.
#' @export
findPeak = function (k, filter, select) {
  estimate = k$estimate
  evalpointsX = k$evalpointsX
  evalpointsY = k$evalpointsY
  if (missing(filter)) {
    filter <- 0
  }

  if (missing(estimate)) {
    stop("can not miss the first argument 'estimate'.\n")
  }
  if (length(estimate[, 1]) < 3 || length(estimate[1, ]) <
      3) {
    stop("first argument 'estimate' should not be smaller than 3*3 matrix.\n")
  }
  markMat <- matrix(0, nrow = length(estimate[, 1]), ncol = length(estimate[1,
  ]))
  for (r in 2:(length(estimate[, 1]) - 1)) {
    for (c in 2:(length(estimate[1, ]) - 1)) {
      if (estimate[r, c] > estimate[r - 1, c - 1] && estimate[r,
                                                              c] > estimate[r - 1, c] && estimate[r, c] >
          estimate[r - 1, c + 1] && estimate[r, c] > estimate[r,
                                                              c - 1] && estimate[r, c] > estimate[r, c + 1] &&
          estimate[r, c] > estimate[r + 1, c - 1] && estimate[r,
                                                              c] > estimate[r + 1, c] && estimate[r, c] >
          estimate[r + 1, c + 1] && estimate[r, c] >=
          filter) {
        markMat[r, c] <- estimate[r, c]
      }
    }
  }

  if (missing(select)) {
    select <- sum(markMat>0)
  }

  if(select < sum(markMat>0)){
    elem = as.vector(markMat)
    sorted_vec = sort(elem, decreasing = TRUE)
    n_largest = sorted_vec[select]

    markMat[markMat<n_largest] = 0
  }

  peaks_index = which(markMat > 0, arr.ind = TRUE)
  colnames(markMat) = 1:dim(markMat)[2]
  rownames(markMat) = 1:dim(markMat)[1]
  inten <- mapply(function(x, y) markMat[as.character(x), as.character(y)], peaks_index[,1], peaks_index[,2])

  peaks = cbind(evalpointsX[peaks_index[,1]],evalpointsY[peaks_index[,2]],inten)

  return(peaks)
}

