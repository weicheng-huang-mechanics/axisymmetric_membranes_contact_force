#ifndef hyperElasticMM_H
#define hyperElasticMM_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class hyperElasticMM
{
public:
	hyperElasticMM(elasticPlate &m_plate, timeStepper &m_stepper, double m_para1, double m_para2);
	~hyperElasticMM();
	void computeFs();
	void computeJs();
    
    void setFirstJacobian();

    VectorXd totalForce;

    double totalEnergyS;
    
private:
	elasticPlate *plate;
    timeStepper *stepper;
	
    int ndof;

    double param1;
    double param2;
    double thickness;
    
    double xk, yk, xkp1, ykp1, l_k;
    double rBar1, rBar2, r_k;

    VectorXd flocal;
    MatrixXd Jss;

    int ind1, ind2;

    Vector2d p, p1;

    VectorXi localDOF;

    VectorXd computeStretchingForce(double xa, double ya, double xb, double yb, double lBar, double rBar);
    MatrixXd computeStretchingJacobian(double xa, double ya, double xb, double yb, double lBar, double rBar);
    VectorXd ListVec(double a1, double a2, double a3, double a4);
    MatrixXd ListMat(VectorXd a1, VectorXd a2, VectorXd a3, VectorXd a4);
};

#endif
