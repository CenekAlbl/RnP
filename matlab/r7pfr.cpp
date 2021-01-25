#include <cstring>
#include <mex.h>
#include <rnp.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs<2||nlhs<6){
		mexErrMsgTxt("Minimum call: [R,w,C,t,f,r] = R7Pfr(X,u)\n X: 3x6 matrix of 3D points \n u: 2x6 matrix of 2D projections\n");
		return;
	}
	const mwSize	*dims = mxGetDimensions(prhs[0]);
	if(dims[0]!=3 || dims[1]!=7){
		mexErrMsgTxt("First argument must be a 3x6 matrix of 3D points");
	}
	dims = mxGetDimensions(prhs[1]);
	if(dims[0]!=2 || dims[1]!=7){
		mexErrMsgTxt("Second argument must be a 2x6 matrix of 2D projections");
	}
	double r0 = 0;
	if(nrhs>2){
		r0= *mxGetPr(prhs[2]);
	}
	int direction = 0;
	if(nrhs>3){
		direction= *mxGetPr(prhs[3]);
	}
	int maxiter = 32;
	if(nrhs>4){
		maxiter = *mxGetPr(prhs[4]);
	}

	//double *X = mxGetPr(prhs[0]);
	//double *u = mxGetPr(prhs[1]);

	Eigen::MatrixXd X = Eigen::Map<Eigen::MatrixXd>(mxGetPr(prhs[0]),3,7);
	Eigen::MatrixXd u = Eigen::Map<Eigen::MatrixXd>(mxGetPr(prhs[1]),2,7);
	
	
	int nsols;


	double *w_out, *t_out, *C_out, *v_out, *f_out, *r_out;

	RSDoublelinCameraPose res;

	Eigen::Vector3d v0;
	v0 << 0,0,0;

	int err =  iterativeRnP<RSDoublelinCameraPose, R7PfIter>(X, u, v0, 7, r0, direction, maxiter, res);
	
	if(!err){
		plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
		plhs[1] = mxCreateDoubleMatrix(3, 1, mxREAL);
		plhs[2] = mxCreateDoubleMatrix(3, 1, mxREAL);
		plhs[3] = mxCreateDoubleMatrix(3, 1, mxREAL);
		plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
		plhs[5] = mxCreateDoubleMatrix(1, 1, mxREAL);

		v_out = mxGetPr(plhs[0]);
		w_out = mxGetPr(plhs[1]);
		C_out = mxGetPr(plhs[2]);
		t_out = mxGetPr(plhs[3]);
		f_out = mxGetPr(plhs[4]);
		r_out = mxGetPr(plhs[5]);

		std::memcpy(w_out, res.w.data(), 3 * sizeof(double));
		std::memcpy(t_out, res.t.data(), 3 * sizeof(double));
		std::memcpy(C_out, res.C.data(), 3 * sizeof(double));
		std::memcpy(v_out, res.v.data(), 3 * sizeof(double));
		std::memcpy(f_out, &res.f, sizeof(double));
		std::memcpy(r_out, &res.rd, sizeof(double));

	}else{

		plhs[0] = mxCreateDoubleMatrix(3, 0, mxREAL);
		plhs[1] = mxCreateDoubleMatrix(3, 0, mxREAL);
		plhs[2] = mxCreateDoubleMatrix(3, 0, mxREAL);
		plhs[3] = mxCreateDoubleMatrix(3, 0, mxREAL);
		plhs[4] = mxCreateDoubleMatrix(1, 0, mxREAL);
		plhs[5] = mxCreateDoubleMatrix(1, 0, mxREAL);

	}
	

}
	
	

	
