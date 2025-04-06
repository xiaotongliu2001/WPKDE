#include<stdlib.h>
#include<stdio.h>
#include<math.h>

struct gridX* makeGridKs(float **matrix, int rows, float * H, float tol, int * gridsize);
struct suppX* makeSupp(float **matrix, int rows, float * H, float tol);
struct gridPts* findGridPts(struct gridX * gridx, struct suppX * suppx, int rows);
struct fhatList * kde(float **matrix, int rows, float * H, int * gridsize, float *w);
float * dmvnorm(float **x, int rows, float * mean, float * sigma);

struct gridX
{
	float * gridx1d;
	float * gridx2d;
	int gridsize[2];
};
struct suppX
{
	float ** xmin;
	float ** xmax;
};
struct gridPts
{
	int ** xmin;
	int ** xmax;
};
struct fhatList
{
	float ** x;
	float * evalpoints1d;    //same with gridx1d in gridX
	float * evalpoints2d;
	float ** estimate;
	float * H;
};

void portal(double * vecmatrix, int * rows, double * vecH, int * vecGridSize, double * vecweight, double * vecestimate,double *vecevalpoints)
{
	//attention:pass matrix as an argument, it will change the actual parameter of this matrix
	//matrix of the data points
	float **matrix;
	matrix = (float **)malloc(sizeof(float *) * (*rows));
	for (int j = 0; j < (*rows); j++) {
		matrix[j] = (float *)malloc(sizeof(float) * 2);
		matrix[j][0] = vecmatrix[j];
		matrix[j][1] = vecmatrix[j + (*rows)];
	}
	
	//weight vector w
	float * w = (float *)malloc(sizeof(float)*(*rows));
	for(int i = 0; i< (*rows);++i){
		w[i] = vecweight[i];
	}
	
	float H[2];
	H[0] = vecH[0];  					//H[0] = 0.0000001;
	H[1] = vecH[1];						//H[1] = 50;
	//int gridsize[2] = { 200,50 };
	int gridsize[2];
	gridsize[0] = vecGridSize[0];
	gridsize[1] = vecGridSize[1];

	struct fhatList * fhatList = kde(matrix, *rows, H, gridsize, w);
	for (int i = 0; i < vecGridSize[1]; ++i) {     //cols
		for (int j = 0; j < vecGridSize[0]; ++j) {  //rows
			vecestimate[i*vecGridSize[0]+j]=fhatList->estimate[j][i];
		}
	}
	for(int i = 0; i < vecGridSize[0]; ++i){
		vecevalpoints[i] = fhatList->evalpoints1d[i];
	}
	for(int j = 0; j < vecGridSize[1]; ++j){
		vecevalpoints[vecGridSize[0]+j] = fhatList->evalpoints2d[j];
	}

	free(fhatList);

}

struct gridX* makeGridKs(float **matrix,int rows,float * H,float tol,int * gridsize ) 
//Generates a 2D grid based on the minimum and maximum values found in a given matrix. The grid is expanded by a specified tolerance and step size.
//Returns: A pointer to a gridX struct containing the generated grid's x and y coordinates and grid size.
{
	float xmin[2] = {0,0}, xmax[2] = {0,0};
	//calculate the max and min value for each dimension
	for (int i = 0; i < rows; i++)
	{
		if (i == 0)
		{
			xmin[0] = matrix[0][0];
			xmax[0] = matrix[0][0];
			xmin[1] = matrix[0][1];
			xmax[1] = matrix[0][1];
		}
		if (matrix[i][0] < xmin[0])
			xmin[0] = matrix[i][0];
		if (matrix[i][0] > xmax[0])
			xmax[0] = matrix[i][0];
		if (matrix[i][1] < xmin[1])
			xmin[1] = matrix[i][1];
		if (matrix[i][1] > xmax[1])
			xmax[1] = matrix[i][1];			 
	}
	xmin[0] = xmin[0] - tol*H[0];
	xmin[1] = xmin[1] - tol*H[1];
	xmax[0] = xmax[0] + tol*H[0];
	xmax[1] = xmax[1] + tol*H[1];

	struct gridX * gridx = (struct gridX*)malloc(sizeof(struct gridX));
	float *gridx1d = (float *)malloc(sizeof(float)*gridsize[0]);
	float *gridx2d = (float *)malloc(sizeof(float)*gridsize[1]);
	float stepsize[2];
	stepsize[0] = (xmax[0] - xmin[0]) / (gridsize[0] - 1);
	stepsize[1] = (xmax[1] - xmin[1]) / (gridsize[1] - 1);

	for (int i = 0; i < gridsize[0]; ++i)
	{
		gridx1d[i] = xmin[0] + i*stepsize[0];
	}
	for (int i = 0; i < gridsize[1]; ++i)
	{
		gridx2d[i] = xmin[1] + i*stepsize[1];
	}
	gridx->gridx1d = gridx1d;
	gridx->gridx2d = gridx2d;
	gridx->gridsize[0] = gridsize[0];
	gridx->gridsize[1] = gridsize[1];
	return gridx;
}

struct suppX* makeSupp(float **matrix, int rows, float * H, float tol)
{
	struct suppX * suppx = (struct suppX*)malloc(sizeof(struct suppX));
	float ** xmin = (float **)malloc(sizeof(float *) * rows);
	float ** xmax = (float **)malloc(sizeof(float *) * rows);
	for (int j = 0; j < rows; j++) {
		xmin[j] = (float *)malloc(sizeof(float) * 2);
		xmax[j] = (float *)malloc(sizeof(float) * 2);
	}

	for (int i = 0; i < rows; ++i) {
		xmin[i][0] = matrix[i][0] - tol*H[0];
		xmin[i][1] = matrix[i][1] - tol*H[1];
		xmax[i][0] = matrix[i][0] + tol*H[0];
		xmax[i][1] = matrix[i][1] + tol*H[1];
	}
	suppx->xmin = xmin;
	suppx->xmax = xmax;
	return suppx;
}

struct gridPts* findGridPts(struct gridX * gridx, struct suppX * suppx,int rows)
{
	struct gridPts * gridpts = (struct gridPts*)malloc(sizeof(struct gridPts));
	int ** gridptsmin = (int **)malloc(sizeof(int *) * rows);
	int ** gridptsmax = (int **)malloc(sizeof(int *) * rows);
	for (int j = 0; j < rows; j++) { 
		gridptsmin[j] = (int *)malloc(sizeof(int) * 2);
		gridptsmax[j] = (int *)malloc(sizeof(int) * 2);
	}
	
	int count;
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < 2; ++j) {
			count = 0;
			if (j == 0) {
				//find index of last element of gridx smaller than min support 
				while (gridx->gridsize[j]>count && suppx->xmin[i][j] >= gridx->gridx1d[count] ) {
					++count;
				}
				if (count == 0)
					gridptsmin[i][j] = 0;
				else
					gridptsmin[i][j] = count-1;

				// find index of first element gridx greater than max support 
				count = 0;
				while (gridx->gridsize[j] > count && suppx->xmax[i][j] >= gridx->gridx1d[count]) {
					++count;
				}
				gridptsmax[i][j] = count-1;     //count should be equal or larger than 1
			}
			else {
				//find index of last element of gridx smaller than min support 
				while (gridx->gridsize[j]>count && suppx->xmin[i][j] >= gridx->gridx2d[count]) {
					++count;
				}
				if (count == 0)
					gridptsmin[i][j] = 0;
				else
					gridptsmin[i][j] = count - 1;

				// find index of first element gridx greater than max support 
				count = 0;
				while (gridx->gridsize[j] > count && suppx->xmax[i][j] >= gridx->gridx2d[count]) {
					++count;
				}
				gridptsmax[i][j] = count - 1;     //count should be equal or larger than 1
			}
		}
	}
	gridpts->xmin = gridptsmin;
	gridpts->xmax = gridptsmax;
	return gridpts;
}

struct fhatList * kde(float **matrix, int rows, float * H, int * gridsize,float *w)
{
	//sum(w) should be equal to rows
	float wsum = 0;
	for (int i = 0; i < rows; i++) {
		wsum = wsum + w[i];
	}
	if (wsum != rows)
	{
		float scale = rows / wsum;
		for (int i = 0; i < rows; ++i) {
			w[i] = w[i] * scale;
		}
	}

	float * H1 = (float *)malloc(sizeof(float) * 2);
	H1[0] = (float)H[0];
	H1[1] = (float)H[1];
	struct gridX * gridx = makeGridKs(matrix, rows, H1, 0, gridsize);
	struct suppX * suppx = makeSupp(matrix, rows, H1, 1);
	struct gridPts * gridpts = findGridPts(gridx, suppx, rows);

	float **fhatgrid = (float **)calloc(gridx->gridsize[0], sizeof(float *));
    for (int j = 0; j < gridx->gridsize[0]; j++) {
        fhatgrid[j] = (float *)calloc(gridx->gridsize[1], sizeof(float));
    }
                                                                        
	float ** evalpts;
	for (int i = 0; i < rows; i++) {
		//compute evaluation points
		int xminind = gridpts->xmin[i][0];
		int xmaxind = gridpts->xmax[i][0];
		int yminind = gridpts->xmin[i][1];
		int ymaxind = gridpts->xmax[i][1];		//record the influenced points' first and last index in both x-axis and y-axis  

		int xindnum = xmaxind - xminind + 1;
		int yindnum = ymaxind - yminind + 1;	//record the influenced points' sum number in both x-axis and y-axis
		int ptsnum = xindnum*yindnum; 			//record the total number of the influenced points
		
		evalpts = (float **)malloc(sizeof(float *) * ptsnum);
		for (int j = 0; j < ptsnum; j++) {
			evalpts[j] = (float *)malloc(sizeof(float) * 2);
			int xind = j / yindnum + xminind;       //record the evaluated point's index  
			int yind = j%yindnum + yminind;
			evalpts[j][0] = gridx->gridx1d[xind];
			evalpts[j][1] = gridx->gridx2d[yind];
		}
		
		//call function dmvnorm
		float * fhat = dmvnorm(evalpts,ptsnum,matrix[i],H);

		// place vector of density estimate values 'fhat' onto grid 'fhat.grid'
		for (int m = 0; m < xindnum; m++) {
			for (int n = 0; n < yindnum; n++) {
				fhatgrid[xminind + m][yminind + n] = fhatgrid[xminind + m][yminind + n] + fhat[m*yindnum + n] * w[i];
			}
		}
		//release memory
		for (int j = 0; j < ptsnum; j++)
			free(evalpts[j]);
		free(evalpts);
		free(fhat);
	}
	//for (int i = 0; i < gridx->gridsize[0]; ++i) {
		//for (int j = 0; j < gridx->gridsize[1]; ++j)
			//fhatgrid[i][j] = fhatgrid[i][j] / rows;
	//}
	struct fhatList * fhatlist = (struct fhatList*)malloc(sizeof(struct fhatList));
	fhatlist->x = matrix;
	fhatlist->evalpoints1d = gridx->gridx1d;
	fhatlist->evalpoints2d = gridx->gridx2d;
	fhatlist->estimate = fhatgrid;
	fhatlist->H = H;

	return fhatlist;
}

float * dmvnorm(float **x,int rows, float * mean, float * sigma)   //'rows' record the rows of x and is equal to xindnum*yindnum														   
{																   //pass the evaluated point's x-coordinate and y-coordinate to mean, pass H to sigma
	float tmp[2][rows];
	float * logretval = (float *)malloc(sizeof(float) * rows);
	float scale = 2*M_PI*sigma[0]*sigma[1];

	for (int i = 0; i < rows; i++) {
		tmp[0][i] = (x[i][0]- mean[0]) / sigma[0];
		tmp[1][i] = (x[i][1] - mean[1]) / sigma[1];
		logretval[i] = exp(- 0.5*(tmp[0][i]*tmp[0][i] + tmp[1][i]*tmp[1][i])) / scale;
	}

	return logretval;    //logretval is a matrix of 1*rows, nrow = 1 and ncol = rows
}