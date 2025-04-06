#' Two-dimensional fast weighted kernel density estimation
#'
#' @param x Data points in the format of an n x 2 matrix.
#' @param H Bandwidth, a vector containing 2 numeric values.
#' @param gridsize Number of points for each direction, a vector containing 2 integer values.
#' @param cutNum Number of pieces to be cut for each direction, a vector containing 2 integer values.
#' @param w Weight, a vector corresponding to parameter 'x'.
#'
#' @return A list containing three elements:
#'   \item{estimate}{The estimated values of the kernel density.}
#'   \item{evalpointsX}{The evaluation points along the X direction.}
#'   \item{evalpointsY}{The evaluation points along the Y direction.}
#' @import Rcpp
#' @useDynLib WPKDE
#' @export
kdeC<-function(x,H,gridsize,cutNum,w){
  if(missing(x)){
    stop("can not miss the first argument 'x'.\n")
  }
  if(missing(H)){
    warning("argument 'H' is missing - it has been set as default (0.01,0.01).\n")
    gridsize<-c(0.01,0.01)
  }
  if(missing(gridsize)){
    warning("argument 'gridsize' is missing - it has been set as default (200,50).\n")
    gridsize<-c(200,50)
  }
  if(missing(cutNum)){
    cutNum<-rep(1,2)
  }
  if(missing(w)){
    w<-rep(1,length(x)/2)
  }
  xmin<-min(x[,1])
  xmax<-max(x[,1])
  ymin<-min(x[,2])
  ymax<-max(x[,2])
  wmin<-min(w)
  xlen<-(xmax-xmin)/cutNum[1]    #length of each cutted piece in x-axis direction
  ylen<-(ymax-ymin)/cutNum[2]    #length of each cutted piece in y-axis direction
  xGridLength<-xlen/(gridsize[1]-1)   #length of each grid in x-axis direction
  yGridLength<-ylen/(gridsize[2]-1)   #length of each grid in y-axis direction
  extendNumX<-floor(2*H[1]/xGridLength)  #calculate the number of grids to be extended in x-axis direction
  extendNumY<-floor(2*H[2]/yGridLength)  #calculate the number of grids to be extended in y-axis direction

  evalpointsX<-numeric(0)
  evalpointsY<-numeric(0)

  if(cutNum[1]==1 & cutNum[2]==1){
    #only 1 piece, do not need to cut
    #estimate<-matrix(nrow = gridsize[1],ncol = gridsize[2])
    result<-.C("portal",vecmatrix = as.double(x),rows = as.integer(length(x)/2),vecH = as.double(H),vecGridSize = as.integer(gridsize),vecweight = as.double(w),vecestimate = as.double(rep(0,gridsize[1]*gridsize[2])),vecevalpoints = as.double(rep(0,gridsize[1]+gridsize[2])))
    estimate<-matrix(result$vecestimate,nrow = gridsize[1],ncol = gridsize[2])
    evalpointsX<-result$vecevalpoints[1:gridsize[1]]
    evalpointsY<-result$vecevalpoints[(gridsize[1]+1):(gridsize[1]+gridsize[2])]
  }
  else if(cutNum[1]>1 & cutNum[2]==1){
    #only to cut in the x-axis direction
    estimate<-matrix(0,nrow = (gridsize[1]*cutNum[1]-cutNum[1]+1),ncol = gridsize[2])
    tempEstimate<-matrix(0,nrow = gridsize[1]+2*extendNumX,ncol = gridsize[2])
    for(i in 1:cutNum[1]){
      if(i ==1){
        index<-which(x[,1]>=xmin & x[,1]<=xmin+xlen+extendNumX*xGridLength)
        selData<-x[index,1:2]
        selData<-rbind(selData,c(xmin,ymin))
        selData<-rbind(selData,c(xmin+xlen+extendNumX*xGridLength,ymax))
        selW<-w[index]
        selW<-c(selW,wmin)
        selW<-c(selW,wmin)
        selGridSize<-gridsize
        selGridSize[1]<-gridsize[1]+extendNumX
        result<-.C("portal",vecmatrix = as.double(selData),rows = as.integer(length(selData)/2),vecH = as.double(H),vecGridSize = as.integer(selGridSize),vecweight = as.double(selW),vecestimate = as.double(rep(0,selGridSize[1]*selGridSize[2])),vecevalpoints = as.double(rep(0,selGridSize[1]+selGridSize[2])))
        #tempEstimate<-matrix(result$vecestimate,nrow = selGridSize[1],ncol = selGridSize[2])
        tempEstimate[1:selGridSize[1],1:selGridSize[2]]<-result$vecestimate
        estimate[1:gridsize[1],1:selGridSize[2]]<-tempEstimate[1:gridsize[1],1:selGridSize[2]]
        evalpointsX<-result$vecevalpoints[1:gridsize[1]]
        evalpointsY<-result$vecevalpoints[(selGridSize[1]+1):(selGridSize[1]+gridsize[2])]
      }
      else if(i>1 & i<cutNum[1]){
        index<-which(x[,1]>=xmin+(i-1)*xlen-extendNumX*xGridLength & x[,1]<=xmin+i*xlen+extendNumX*xGridLength)
        selData<-x[index,1:2]
        selData<-rbind(selData,c(xmin+(i-1)*xlen-extendNumX*xGridLength,ymin))
        selData<-rbind(selData,c(xmin+i*xlen+extendNumX*xGridLength,ymax))
        selW<-w[index]
        selW<-c(selW,wmin)
        selW<-c(selW,wmin)
        selGridSize<-gridsize
        selGridSize[1]<-gridsize[1]+extendNumX*2
        result<-.C("portal",vecmatrix = as.double(selData),rows = as.integer(length(selData)/2),vecH = as.double(H),vecGridSize = as.integer(selGridSize),vecweight = as.double(selW),vecestimate = as.double(rep(0,selGridSize[1]*selGridSize[2])),vecevalpoints = as.double(rep(0,selGridSize[1]+selGridSize[2])))
        xFirstIndex<-gridsize[1]*(i-1)-i+2
        xLastIndex<-gridsize[1]*i-i+1
        tempEstimate[1:selGridSize[1],1:selGridSize[2]]<-result$vecestimate
        estimate[xFirstIndex:xLastIndex,1:selGridSize[2]]<-tempEstimate[(extendNumX+1):(extendNumX+gridsize[1]),1:selGridSize[2]]
        evalpointsX<-c(evalpointsX,result$vecevalpoints[(extendNumX+2):(extendNumX+gridsize[1])])
      }
      else if(i == cutNum[1]){
        index<-which(x[,1]>=xmin+(i-1)*xlen-extendNumX*xGridLength & x[,1]<=xmax)
        selData<-x[index,1:2]
        selData<-rbind(selData,c(xmin+(i-1)*xlen-extendNumX*xGridLength,ymin))
        selData<-rbind(selData,c(xmax,ymax))
        selW<-w[index]
        selW<-c(selW,wmin)
        selW<-c(selW,wmin)
        selGridSize<-gridsize
        selGridSize[1]<-gridsize[1]+extendNumX
        result<-.C("portal",vecmatrix = as.double(selData),rows = as.integer(length(selData)/2),vecH = as.double(H),vecGridSize = as.integer(selGridSize),vecweight = as.double(selW),vecestimate = as.double(rep(0,selGridSize[1]*selGridSize[2])),vecevalpoints = as.double(rep(0,selGridSize[1]+selGridSize[2])))
        xFirstIndex<-gridsize[1]*(i-1)-i+2
        tempEstimate[1:(gridsize[1]+extendNumX),1:selGridSize[2]]<-result$vecestimate
        estimate[xFirstIndex:(gridsize[1]*cutNum[1]-cutNum[1]+1),1:selGridSize[2]]<-tempEstimate[(extendNumX+1):(extendNumX+gridsize[1]),1:gridsize[2]]
        evalpointsX<-c(evalpointsX,result$vecevalpoints[(extendNumX+2):(extendNumX+gridsize[1])])
      }
    }
  }
  else if(cutNum[1]==1 & cutNum[2]>1){
    #only to cut in the y-axis direction
    estimate<-matrix(0,nrow = gridsize[1],ncol = gridsize[2]*cutNum[2]-cutNum[2]+1)
    tempEstimate<-matrix(0,nrow=gridsize[1],ncol=gridsize[2]+extendNumY*2)
    for(i in 1:cutNum[2]){
      if(i ==1){
        index<-which(x[,2]>=ymin & x[,2]<=ymin+ylen+extendNumY*yGridLength)
        selData<-x[index,1:2]
        selData<-rbind(selData,c(xmin,ymin))
        selData<-rbind(selData,c(xmax,ymin+ylen+extendNumY*yGridLength))
        selW<-w[index]
        selW<-c(selW,wmin)
        selW<-c(selW,wmin)
        selGridSize<-gridsize
        selGridSize[2]<-gridsize[2]+extendNumY
        result<-.C("portal",vecmatrix = as.double(selData),rows = as.integer(length(selData)/2),vecH = as.double(H),vecGridSize = as.integer(selGridSize),vecweight = as.double(selW),vecestimate = as.double(rep(0,selGridSize[1]*selGridSize[2])),vecevalpoints = as.double(rep(0,selGridSize[1]+selGridSize[2])))
        #tempEstimate<-matrix(result$vecestimate,nrow = selGridSize[1],ncol = selGridSize[2])
        tempEstimate[1:selGridSize[1],1:selGridSize[2]]<-result$vecestimate
        estimate[1:selGridSize[1],1:gridsize[2]]<-tempEstimate[1:selGridSize[1],1:gridsize[2]]
        evalpointsX<-result$vecevalpoints[1:gridsize[1]]
        evalpointsY<-result$vecevalpoints[(gridsize[1]+1):(gridsize[1]+gridsize[2])]
      }
      else if(i>1 & i<cutNum[2]){
        index<-which(x[,2]>=ymin+(i-1)*ylen-extendNumY*yGridLength & x[,2]<=ymin+i*ylen+extendNumY*yGridLength)
        selData<-x[index,1:2]
        selData<-rbind(selData,c(xmin,ymin+(i-1)*ylen-extendNumY*yGridLength))
        selData<-rbind(selData,c(xmax,ymin+i*ylen+extendNumY*yGridLength))
        selW<-w[index]
        selW<-c(selW,wmin)
        selW<-c(selW,wmin)
        selGridSize<-gridsize
        selGridSize[2]<-gridsize[2]+extendNumY*2
        result<-.C("portal",vecmatrix = as.double(selData),rows = as.integer(length(selData)/2),vecH = as.double(H),vecGridSize = as.integer(selGridSize),vecweight = as.double(selW),vecestimate = as.double(rep(0,selGridSize[1]*selGridSize[2])),vecevalpoints = as.double(rep(0,selGridSize[1]+selGridSize[2])))
        yFirstIndex<-gridsize[2]*(i-1)-i+2
        yLastIndex<-gridsize[2]*i-i+1
        tempEstimate[1:selGridSize[1],1:selGridSize[2]]<-result$vecestimate
        estimate[1:selGridSize[1],yFirstIndex:yLastIndex]<-tempEstimate[1:selGridSize[1],(extendNumY+1):(extendNumY+gridsize[2])]
        evalpointsY<-c(evalpointsY,result$vecevalpoints[(gridsize[1]+extendNumY+2):(gridsize[1]+extendNumY+gridsize[2])])
      }
      else if(i == cutNum[2]){
        index<-which(x[,2]>=ymin+(i-1)*ylen-extendNumY*yGridLength & x[,2]<=ymax)
        selData<-x[index,1:2]
        selData<-rbind(selData,c(xmin,ymin+(i-1)*ylen-extendNumY*yGridLength))
        selData<-rbind(selData,c(xmax,ymax))
        selW<-w[index]
        selW<-c(selW,wmin)
        selW<-c(selW,wmin)
        selGridSize<-gridsize
        selGridSize[2]<-gridsize[2]+extendNumY
        result<-.C("portal",vecmatrix = as.double(selData),rows = as.integer(length(selData)/2),vecH = as.double(H),vecGridSize = as.integer(selGridSize),vecweight = as.double(selW),vecestimate = as.double(rep(0,selGridSize[1]*selGridSize[2])),vecevalpoints = as.double(rep(0,selGridSize[1]+selGridSize[2])))
        yFirstIndex<-gridsize[2]*(i-1)-i+2
        tempEstimate[1:selGridSize[1],1:(extendNumY+gridsize[2])]<-result$vecestimate
        estimate[1:selGridSize[1],yFirstIndex:(gridsize[2]*cutNum[2]-cutNum[2]+1)]<-tempEstimate[1:selGridSize[1],(extendNumY+1):(extendNumY+gridsize[2])]
        evalpointsY<-c(evalpointsY,result$vecevalpoints[(gridsize[1]+extendNumY+2):(gridsize[1]+extendNumY+gridsize[2])])
      }
    }
  }
  else if(cutNum[1]>1 & cutNum[2]>1){
    #need to cut in both direction
    estimate<-matrix(0,nrow = gridsize[1]*cutNum[1]-cutNum[1]+1,ncol = gridsize[2]*cutNum[2]-cutNum[2]+1)
    tempEstimate<-matrix(0,nrow = gridsize[1]+extendNumX*2,ncol=gridsize[2]+extendNumY*2)
    for(j in 1:cutNum[2]){
      if(j == 1){
        for(i in 1:cutNum[1]){
          if(i == 1){
            index<-which(x[,1]>=xmin & x[,1]<=(xmin+xlen+extendNumX*xGridLength) & x[,2]>=ymin & x[,2]<=(ymin+ylen+extendNumY*yGridLength))
            selData<-x[index,1:2]
            selData<-rbind(selData,c(xmin+xlen+extendNumX*xGridLength,ymin+ylen+extendNumY*yGridLength))
            selData<-rbind(selData,c(xmin,ymin))
            selW<-w[index]
            selW<-c(selW,wmin)
            selW<-c(selW,wmin)
            selGridSize<-gridsize
            selGridSize[1]<-gridsize[1]+extendNumX
            selGridSize[2]<-gridsize[2]+extendNumY
            result<-.C("portal",vecmatrix = as.double(selData),rows = as.integer(length(selData)/2),vecH = as.double(H),vecGridSize = as.integer(selGridSize),vecweight = as.double(selW),vecestimate = as.double(rep(0,selGridSize[1]*selGridSize[2])),vecevalpoints = as.double(rep(0,selGridSize[1]+selGridSize[2])))
            #tempEstimate<-matrix(result$vecestimate,nrow = selGridSize[1],ncol = selGridSize[2])
            tempEstimate[1:selGridSize[1],1:selGridSize[2]]<-result$vecestimate
            estimate[1:gridsize[1],1:gridsize[2]]<-tempEstimate[1:gridsize[1],1:gridsize[2]]
            evalpointsX<-result$vecevalpoints[1:gridsize[1]]
            evalpointsY<-result$vecevalpoints[(selGridSize[1]+1):(selGridSize[1]+gridsize[2])]
          }
          else if(i>1 & i<cutNum[1]){
            index<-which(x[,1]>=(xmin+(i-1)*xlen-extendNumX*xGridLength) & x[,1]<=(xmin+i*xlen+extendNumX*xGridLength) & x[,2]>=ymin & x[,2]<=(ymin+ylen+extendNumY*yGridLength))
            selData<-x[index,1:2]
            selData<-rbind(selData,c(xmin+(i-1)*xlen-extendNumX*xGridLength,ymin))
            selData<-rbind(selData,c(xmin+i*xlen+extendNumX*xGridLength,ymin+ylen+extendNumY*yGridLength))
            selW<-w[index]
            selW<-c(selW,wmin)
            selW<-c(selW,wmin)
            selGridSize<-gridsize
            selGridSize[1]<-gridsize[1]+extendNumX*2
            selGridSize[2]<-gridsize[2]+extendNumY
            result<-.C("portal",vecmatrix = as.double(selData),rows = as.integer(length(selData)/2),vecH = as.double(H),vecGridSize = as.integer(selGridSize),vecweight = as.double(selW),vecestimate = as.double(rep(0,selGridSize[1]*selGridSize[2])),vecevalpoints = as.double(rep(0,selGridSize[1]+selGridSize[2])))
            xFirstIndex<-gridsize[1]*(i-1)-i+2
            xLastIndex<-gridsize[1]*i-i+1
            tempEstimate[1:selGridSize[1],1:selGridSize[2]]<-result$vecestimate
            estimate[xFirstIndex:xLastIndex,1:gridsize[2]]<-tempEstimate[(extendNumX+1):(extendNumX+gridsize[1]),1:gridsize[2]]
            evalpointsX<-c(evalpointsX,result$vecevalpoints[(extendNumX+2):(extendNumX+gridsize[1])])
          }
          else if(i == cutNum[1]){
            index<-which(x[,1]>=xmin+(i-1)*xlen-extendNumX*xGridLength & x[,1]<=xmax & x[,2]>=ymin & x[,2]<=(ymin+ylen+extendNumY*yGridLength))
            selData<-x[index,1:2]
            selData<-rbind(selData,c(xmin+(i-1)*xlen-extendNumX*xGridLength,ymin))
            selData<-rbind(selData,c(xmax,ymin+ylen+extendNumY*yGridLength))
            selW<-w[index]
            selW<-c(selW,wmin)
            selW<-c(selW,wmin)
            selGridSize<-gridsize
            selGridSize[1]<-gridsize[1]+extendNumX
            selGridSize[2]<-gridsize[2]+extendNumY
            result<-.C("portal",vecmatrix = as.double(selData),rows = as.integer(length(selData)/2),vecH = as.double(H),vecGridSize = as.integer(selGridSize),vecweight = as.double(selW),vecestimate = as.double(rep(0,selGridSize[1]*selGridSize[2])),vecevalpoints = as.double(rep(0,selGridSize[1]+selGridSize[2])))
            xFirstIndex<-gridsize[1]*(i-1)-i+2
            tempEstimate[1:selGridSize[1],1:selGridSize[2]]<-result$vecestimate
            estimate[xFirstIndex:(gridsize[1]*cutNum[1]-cutNum[1]+1),1:gridsize[2]]<-tempEstimate[(extendNumX+1):(extendNumX+gridsize[1]),1:gridsize[2]]
            evalpointsX<-c(evalpointsX,result$vecevalpoints[(extendNumX+2):(extendNumX+gridsize[1])])
          }
        }
      }
      else if(j>1 & j<cutNum[2]){
        for(i in 1:cutNum[1]){
          if(i == 1){
            index<-which(x[,1]>=xmin & x[,1]<=(xmin+xlen+extendNumX*xGridLength) & x[,2]>=(ymin+(j-1)*ylen-extendNumY*yGridLength) & x[,2]<=(ymin+j*ylen+extendNumY*yGridLength))
            selData<-x[index,1:2]
            selData<-rbind(selData,c(xmin+xlen+extendNumX*xGridLength,ymin+j*ylen+extendNumY*yGridLength))
            selData<-rbind(selData,c(xmin,ymin+(j-1)*ylen-extendNumY*yGridLength))
            selW<-w[index]
            selW<-c(selW,wmin)
            selW<-c(selW,wmin)
            selGridSize<-gridsize
            selGridSize[1]<-gridsize[1]+extendNumX
            selGridSize[2]<-gridsize[2]+extendNumY*2
            result<-.C("portal",vecmatrix = as.double(selData),rows = as.integer(length(selData)/2),vecH = as.double(H),vecGridSize = as.integer(selGridSize),vecweight = as.double(selW),vecestimate = as.double(rep(0,selGridSize[1]*selGridSize[2])),vecevalpoints = as.double(rep(0,selGridSize[1]+selGridSize[2])))
            #tempEstimate<-matrix(result$vecestimate,nrow = selGridSize[1],ncol = selGridSize[2])
            yFirstIndex<-gridsize[2]*(j-1)-j+2
            yLastIndex<-gridsize[2]*j-j+1
            tempEstimate[1:selGridSize[1],1:selGridSize[2]]<-result$vecestimate
            estimate[1:gridsize[1],yFirstIndex:yLastIndex]<-tempEstimate[1:gridsize[1],(extendNumY+1):(extendNumY+gridsize[2])]
            evalpointsY<-c(evalpointsY,result$vecevalpoints[(selGridSize[1]+extendNumY+2):(selGridSize[1]+extendNumY+gridsize[2])])
          }
          else if(i>1 & i<cutNum[1]){
            index<-which(x[,1]>=(xmin+(i-1)*xlen-extendNumX*xGridLength) & x[,1]<=(xmin+i*xlen+extendNumX*xGridLength) & x[,2]>=(ymin+(j-1)*ylen-extendNumY*yGridLength) & x[,2]<=(ymin+j*ylen+extendNumY*yGridLength))
            selData<-x[index,1:2]
            selData<-rbind(selData,c(xmin+i*xlen+extendNumX*xGridLength,ymin+j*ylen+extendNumY*yGridLength))
            selData<-rbind(selData,c(xmin+(i-1)*xlen-extendNumX*xGridLength,ymin+(j-1)*ylen-extendNumY*yGridLength))
            selW<-w[index]
            selW<-c(selW,wmin)
            selW<-c(selW,wmin)
            selGridSize<-gridsize
            selGridSize[1]<-gridsize[1]+extendNumX*2
            selGridSize[2]<-gridsize[2]+extendNumY*2
            result<-.C("portal",vecmatrix = as.double(selData),rows = as.integer(length(selData)/2),vecH = as.double(H),vecGridSize = as.integer(selGridSize),vecweight = as.double(selW),vecestimate = as.double(rep(0,selGridSize[1]*selGridSize[2])),vecevalpoints = as.double(rep(0,selGridSize[1]+selGridSize[2])))
            #tempEstimate<-matrix(result$vecestimate,nrow = selGridSize[1],ncol = selGridSize[2])
            yFirstIndex<-gridsize[2]*(j-1)-j+2
            yLastIndex<-gridsize[2]*j-j+1
            xFirstIndex<-gridsize[1]*(i-1)-i+2
            xLastIndex<-gridsize[1]*i-i+1
            tempEstimate[1:selGridSize[1],1:selGridSize[2]]<-result$vecestimate
            estimate[xFirstIndex:xLastIndex,yFirstIndex:yLastIndex]<-tempEstimate[(extendNumX+1):(extendNumX+gridsize[1]),(extendNumY+1):(extendNumY+gridsize[2])]
          }
          else if(i == cutNum[1]){
            index<-which(x[,1]>=(xmin+(i-1)*xlen-extendNumX*xGridLength) & x[,1]<=xmax & x[,2]>=(ymin+(j-1)*ylen-extendNumY*yGridLength) & x[,2]<=(ymin+j*ylen+extendNumY*yGridLength))
            selData<-x[index,1:2]
            selData<-rbind(selData,c(xmax,ymin+j*ylen+extendNumY*yGridLength))
            selData<-rbind(selData,c(xmin+(i-1)*xlen-extendNumX*xGridLength,ymin+(j-1)*ylen-extendNumY*yGridLength))
            selW<-w[index]
            selW<-c(selW,wmin)
            selW<-c(selW,wmin)
            selGridSize<-gridsize
            selGridSize[1]<-gridsize[1]+extendNumX
            selGridSize[2]<-gridsize[2]+extendNumY*2
            result<-.C("portal",vecmatrix = as.double(selData),rows = as.integer(length(selData)/2),vecH = as.double(H),vecGridSize = as.integer(selGridSize),vecweight = as.double(selW),vecestimate = as.double(rep(0,selGridSize[1]*selGridSize[2])),vecevalpoints = as.double(rep(0,selGridSize[1]+selGridSize[2])))
            #tempEstimate<-matrix(result$vecestimate,nrow = selGridSize[1],ncol = selGridSize[2])
            yFirstIndex<-gridsize[2]*(j-1)-j+2
            yLastIndex<-gridsize[2]*j-j+1
            xFirstIndex<-gridsize[1]*(i-1)-i+2
            xLastIndex<-gridsize[1]*i-i+1
            tempEstimate[1:selGridSize[1],1:selGridSize[2]]<-result$vecestimate
            estimate[xFirstIndex:xLastIndex,yFirstIndex:yLastIndex]<-tempEstimate[(extendNumX+1):(extendNumX+gridsize[1]),(extendNumY+1):(extendNumY+gridsize[2])]
          }
        }
      }
      else if(j == cutNum[2]){
        for(i in 1:cutNum[1]){
          if(i == 1){
            index<-which(x[,1]>=xmin & x[,1]<=(xmin+xlen+extendNumX*xGridLength) & x[,2]>=(ymin+(j-1)*ylen-extendNumY*yGridLength) & x[,2]<=ymax)
            selData<-x[index,1:2]
            selData<-rbind(selData,c(xmin+xlen+extendNumX*xGridLength,ymax))
            selData<-rbind(selData,c(xmin,ymin+(j-1)*ylen-extendNumY*yGridLength))
            selW<-w[index]
            selW<-c(selW,wmin)
            selW<-c(selW,wmin)
            selGridSize<-gridsize
            selGridSize[1]<-gridsize[1]+extendNumX
            selGridSize[2]<-gridsize[2]+extendNumY
            result<-.C("portal",vecmatrix = as.double(selData),rows = as.integer(length(selData)/2),vecH = as.double(H),vecGridSize = as.integer(selGridSize),vecweight = as.double(selW),vecestimate = as.double(rep(0,selGridSize[1]*selGridSize[2])),vecevalpoints = as.double(rep(0,selGridSize[1]+selGridSize[2])))
            #tempEstimate<-matrix(result$vecestimate,nrow = selGridSize[1],ncol = selGridSize[2])
            yFirstIndex<-gridsize[2]*(j-1)-j+2
            yLastIndex<-gridsize[2]*j-j+1
            tempEstimate[1:selGridSize[1],1:selGridSize[2]]<-result$vecestimate
            estimate[1:gridsize[1],yFirstIndex:yLastIndex]<-tempEstimate[1:gridsize[1],(extendNumY+1):(extendNumY+gridsize[2])]
            evalpointsY<-c(evalpointsY,result$vecevalpoints[(selGridSize[1]+extendNumY+2):(selGridSize[1]+extendNumY+gridsize[2])])
          }
          else if(i>1 & i<cutNum[1]){
            index<-which(x[,1]>=(xmin+(i-1)*xlen-extendNumX*xGridLength) & x[,1]<=(xmin+i*xlen+extendNumX*xGridLength) & x[,2]>=(ymin+(j-1)*ylen-extendNumY*yGridLength) & x[,2]<=ymax)
            selData<-x[index,1:2]
            selData<-rbind(selData,c(xmin+(i-1)*xlen-extendNumX*xGridLength,ymin+(j-1)*ylen-extendNumY*yGridLength))
            selData<-rbind(selData,c(xmin+i*xlen+extendNumX*xGridLength,ymax))
            selW<-w[index]
            selW<-c(selW,wmin)
            selW<-c(selW,wmin)
            selGridSize<-gridsize
            selGridSize[1]<-gridsize[1]+extendNumX*2
            selGridSize[2]<-gridsize[2]+extendNumY
            result<-.C("portal",vecmatrix = as.double(selData),rows = as.integer(length(selData)/2),vecH = as.double(H),vecGridSize = as.integer(selGridSize),vecweight = as.double(selW),vecestimate = as.double(rep(0,selGridSize[1]*selGridSize[2])),vecevalpoints = as.double(rep(0,selGridSize[1]+selGridSize[2])))
            xFirstIndex<-gridsize[1]*(i-1)-i+2
            xLastIndex<-gridsize[1]*i-i+1
            yFirstIndex<-gridsize[2]*(j-1)-j+2
            yLastIndex<-gridsize[2]*j-j+1
            tempEstimate[1:selGridSize[1],1:selGridSize[2]]<-result$vecestimate
            estimate[xFirstIndex:xLastIndex,yFirstIndex:yLastIndex]<-tempEstimate[(extendNumX+1):(extendNumX+gridsize[1]),(extendNumY+1):(extendNumY+gridsize[2])]
          }
          else if(i == cutNum[1]){
            index<-which(x[,1]>=xmin+(i-1)*xlen-extendNumX*xGridLength & x[,1]<=xmax & x[,2]>=(ymin+(j-1)*ylen-extendNumY*yGridLength) & x[,2]<=ymax)
            selData<-x[index,1:2]
            selData<-rbind(selData,c(xmin+(i-1)*xlen-extendNumX*xGridLength,ymin+(j-1)*ylen-extendNumY*yGridLength))
            selData<-rbind(selData,c(xmax,ymax))
            selW<-w[index]
            selW<-c(selW,wmin)
            selW<-c(selW,wmin)
            selGridSize<-gridsize
            selGridSize[1]<-gridsize[1]+extendNumX
            selGridSize[2]<-gridsize[2]+extendNumY
            result<-.C("portal",vecmatrix = as.double(selData),rows = as.integer(length(selData)/2),vecH = as.double(H),vecGridSize = as.integer(selGridSize),vecweight = as.double(selW),vecestimate = as.double(rep(0,selGridSize[1]*selGridSize[2])),vecevalpoints = as.double(rep(0,selGridSize[1]+selGridSize[2])))
            xFirstIndex<-gridsize[1]*(i-1)-i+2
            xLastIndex<-gridsize[1]*i-i+1
            yFirstIndex<-gridsize[2]*(j-1)-j+2
            yLastIndex<-gridsize[2]*j-j+1
            tempEstimate[1:selGridSize[1],1:selGridSize[2]]<-result$vecestimate
            estimate[xFirstIndex:xLastIndex,yFirstIndex:yLastIndex]<-tempEstimate[(extendNumX+1):(extendNumX+gridsize[1]),(extendNumY+1):(extendNumY+gridsize[2])]
          }
        }
      }
    }
  }
  result<-list(estimate=estimate,evalpointsX=evalpointsX,evalpointsY=evalpointsY)
  return(result)
}
