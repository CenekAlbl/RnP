#include <cstring>
#include <mex.h>
#include <rnp.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs<2||nlhs<4){
		mexErrMsgTxt("Minimum call: [R,w,C,t] = R6P2lin(X,u)\n X: 3x6 matrix of 3D points \n u: 2x6 matrix of 2D projections\n");
		return;
	}
	const mwSize	*dims = mxGetDimensions(prhs[0]);
	if(dims[0]!=3 || dims[1]!=6){
		mexErrMsgTxt("First argument must be a 3x6 matrix of 3D points");
	}
	dims = mxGetDimensions(prhs[1]);
	if(dims[0]!=2 || dims[1]!=6){
		mexErrMsgTxt("Second argument must be a 2x6 matrix of 2D projections");
	}
	double *ud;
	if(nrhs>2){
		dims = mxGetDimensions(prhs[2]);
		if(dims[0]!=2 || dims[1]!=6){
			mexErrMsgTxt("Third argument must be a 2x6 matrix of 2D distorted projections");
		}
		ud = mxGetPr(prhs[2]);
	}	
	double r0 = 0;
	if(nrhs>3){
		r0= *mxGetPr(prhs[3]);
	}
	int direction = 0;
	if(nrhs>4){
		direction= *mxGetPr(prhs[4]);
	}
	int maxpow = 32;
	if(nrhs>5){
		maxpow = *mxGetPr(prhs[5]);
	}

	//double *X = mxGetPr(prhs[0]);
	//double *u = mxGetPr(prhs[1]);

	Eigen::MatrixXd X = Eigen::Map<Eigen::MatrixXd>(mxGetPr(prhs[0]),3,6);
	Eigen::MatrixXd u = Eigen::Map<Eigen::MatrixXd>(mxGetPr(prhs[1]),2,6);
	
	
	int nsols;


	double *w_out, *t_out, *C_out, *v_out;

	RSDoublelinCameraPoseVector results;


	R6P2Lin(X, u, direction, r0, &results);

	nsols = 0;

	for(auto res : results){
		nsols++;
	}
	
	plhs[0] = mxCreateDoubleMatrix(3, nsols, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(3, nsols, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(3, nsols, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(3, nsols, mxREAL);
	

	v_out = mxGetPr(plhs[0]);
	w_out = mxGetPr(plhs[1]);
	C_out = mxGetPr(plhs[2]);	
	t_out = mxGetPr(plhs[3]);

	for(auto res : results){
		std::memcpy(w_out, res.w.data(), 3 * sizeof(double));
		std::memcpy(t_out, res.t.data(), 3 * sizeof(double));
		std::memcpy(C_out, res.C.data(), 3 * sizeof(double));
		std::memcpy(v_out, res.v.data(), 3 * sizeof(double));
		w_out += 3;
		t_out += 3;
		C_out += 3;
		v_out += 3;
					
	}
	
	
	

}
	
	

	
