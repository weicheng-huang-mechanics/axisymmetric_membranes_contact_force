#ifndef EXTERNALCONATCTFORCE_H
#define EXTERNALCONATCTFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class externalContactForce
{
public:
	externalContactForce(elasticPlate &m_plate, timeStepper &m_stepper, 
        double m_stiffness, double m_dBar);
	~externalContactForce();

	void computeFc();
	void computeJc();

    VectorXd ForceVec;

    double totalContactForce;

    double upperBound;

private:
	elasticPlate *plate;
    timeStepper *stepper;

    double stiffness;
    double dBar;

    int ind;

    double dEnergydD;
    double d2EnergydD2;

    Vector2d f;

    double dt;

    double fVelocity;
    Vector2d tK;

    Vector2d friction;
    Matrix2d frictionJacobian;
    Vector2d dfVelocity;
    Matrix2d dtK;

    Matrix2d Id3;
    Matrix2d IdG;
};

#endif
