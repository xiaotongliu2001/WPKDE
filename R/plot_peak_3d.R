#' Plot of the 3D data points with peaks highlighted in green
#'
#' This function creates an interactive 3D scatter plot of data points and highlights the peaks
#' that are within a specified tolerance distance from any data point.
#'
#'
#' @param dat A numeric matrix or data frame with at least three columns representing x, y, and z coordinates of data points.
#' @param peaks A numeric matrix or data frame with at least two columns representing the x and y coordinates of peak candidates.
#' @param x.range A numeric vector of length 2 specifying the x-axis range to include.
#' @param y.range A numeric vector of length 2 specifying the y-axis range to include.
#' @param tol A numeric value specifying the tolerance threshold: only peaks within this Euclidean distance from a data point are retained.
#' @import plotly
#' @import RANN
#' @export
plot_peak_3d = function(dat, peaks,x.range=NA, y.range=NA, tol=1e-5){
  p = plot_ly(x=dat[,1],y=dat[,2],z=dat[,3],type = 'scatter3d',mode = 'markers', marker = list(size = 1, color = 'black'),name="Data points")
  p = p %>%
    layout(
      scene = list(
        xaxis = list(range = x.range),
        yaxis = list(range = y.range)
      )
    )

  if(!any(is.na(x.range))){
    peaks = peaks[peaks[,1]<x.range[2]&peaks[,1]>x.range[1],]
  }
  if(!any(is.na(y.range))){
    peaks = peaks[peaks[,2]<y.range[2]&peaks[,2]>y.range[1],]
  }
  if(length(peaks)==0){
    cat("No peaks were detected within the specified range.\n")
    print(p)
  }
  else{
    nearest = nn2(dat[,1:2], peaks[,1:2,drop = FALSE], k = 1)
    cloest_peak = nearest$nn.dists<tol
    cloest_peak_index = nearest$nn.idx[cloest_peak]

    if(length(cloest_peak_index)==0){
      cat("No peaks were retained within the specified tolerance.\n")
    }
    else{
      cat(sprintf("Found %d peaks in total, %d peaks were retained\n", nrow(peaks), length(cloest_peak_index)))
      p = p %>% add_trace(x = dat[cloest_peak_index,1], y = dat[cloest_peak_index,2], z = dat[cloest_peak_index,3], mode = 'markers', marker = list(size = 4, color = 'green'),name="Peaks")
    }
    print(p)
  }
}
