
#include "prepareDoubleLin.h"
#include "utils.h"

int r6pDoubleLin(double * X, double * u, double r0, double * C, double * t, double *v, double * w, int direction){
	Eigen::MatrixXd H(12, 22);
	double X_1, X_2, X_3, r, c, c0;
	bool planar = false;
	if (direction!=0){
		c0 = r0;
	}
	

	for (int i = 0; i < 6; i++)
	{
		X_1 = X[i * 3];
		X_2 = X[i * 3 + 1];
		X_3 = X[i * 3 + 2];
		r = u[i * 2];
		c = u[i * 2 + 1];

		if (direction == 0){
			H.row(i * 2) << 0, -1, c, 0, r0 - r, c*r - c*r0, X_2*r - X_2*r0 - X_3*c*r + X_3*c*r0, 0, 0, X_1*r0 - X_1*r, X_3*c*r0 - X_3*c*r, X_3*r0 - X_3*r, X_1*c*r - X_1*c*r0, X_2*c*r - X_2*c*r0, X_2*r - X_2*r0, X_3 + X_2*c, -X_1*c, -X_1, X_3*r - X_3*r0 + X_2*c*r - X_2*c*r0, X_1*c*r0 - X_1*c*r, X_1*r0 - X_1*r, X_3*c - X_2;
			H.row(i * 2 + 1) << 1, 0, -r, r - r0, 0, -r*r + r0*r, X_3*r*r - X_3*r0*r, X_2*r - X_2*r0, X_3*r - X_3*r0, 0, X_1*r0 - X_1*r + X_3*r*r - X_3*r*r0, 0, -X_1*r*r + X_1*r0*r, -X_2*r*r + X_2*r0*r, X_1*r0 - X_1*r, -X_2*r, X_3 + X_1*r, -X_2, -X_2*r*r + X_2*r0*r, X_3*r - X_3*r0 + X_1*r*r - X_1*r*r0, X_2*r0 - X_2*r, X_1 - X_3*r;
		}
		else{
			H.row(i * 2) << 0, -1, c, 0, c0 - c, c*c - c0*c, X_2*c - X_2*c0 - X_3*c*c + X_3*c*c0, 0, 0, X_1*c0 - X_1*c, -X_3*c*c + X_3*c0*c, X_3*c0 - X_3*c, X_1*c*c - X_1*c0*c, X_2*c*c - X_2*c0*c, X_2*c - X_2*c0, X_3 + X_2*c, -X_1*c, -X_1, X_3*c - X_3*c0 + X_2*c*c - X_2*c*c0, -X_1*c*c + X_1*c0*c, X_1*c0 - X_1*c, X_3*c - X_2;
			H.row(i * 2 + 1) << -c, r, 0, -c*c + c0*c, c*r - c0*r, 0, X_2*c0*r - X_2*c*r, -X_2*c*c + X_2*c0*c, -X_3*c*c+ X_3*c0*c, X_1*c*r - X_1*c0*r, X_1*c*c - X_1*c0*c, X_3*c*r - X_3*c0*r, 0, 0, X_1*c*c - X_1*c*c0 - X_2*c*r + X_2*c0*r, -X_3*r, -X_3*c, X_2*c + X_1*r, X_3*c0*r - X_3*c*r, -X_3*c*c + X_3*c0*c, X_2*c*c - X_2*c*c0 + X_1*c*r - X_1*c0*r, X_2*r - X_1*c;
		}

	}
	
	Eigen::MatrixXd Helim = H.transpose();
	std::list<int> b;
	colEchelonForm(Helim, b);

	Eigen::MatrixXd Helimtr = Helim.transpose();
	
	Eigen::MatrixXd g; 
	Eigen::VectorXd p;
	std::vector<double> v1, v2, v3, w1, w2, w3;
	
    g = Helimtr.block(6, 12, 6, 10).transpose();
    p = Eigen::Map<Eigen::VectorXd>(g.data(), 6 * 10);
    solveVW3var(p.data(), v1, v2, v3, w1, w2, w3);
	
	Eigen::VectorXd x(16,1);
	Eigen::MatrixXd A = Helimtr.block(0,6,6,16);
	Eigen::VectorXd Ct(6);
	Eigen::Vector3d Cres;

	for (int i = 0; i < v1.size(); i++)
	{
		x << v1[i]*w1[i], v1[i]*w2[i], v1[i]*w3[i], v2[i]*w1[i], v2[i]*w2[i], v2[i]*w3[i], v3[i]*w1[i], v3[i]*w2[i], v3[i]*w3[i], v1[i], v2[i], v3[i], w1[i], w2[i], w3[i], 1;
		Ct = -A*x;
		memcpy(C+i*3,Ct.segment(0, 3).data(),3*sizeof(double));
		memcpy(t + i * 3, Ct.segment(3, 3).data(), 3 * sizeof(double));
		Eigen::Vector3d vr;
		vr << v1[i], v2[i], v3[i];
		memcpy(v + i * 3, vr.data(), 3 * sizeof(double));
		w[i * 3] = w1[i];
		w[i * 3 + 1] = w2[i];
		w[i * 3 + 2] = w3[i];
	}

	int nresults = v1.size();
    return nresults;
}