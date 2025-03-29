#include "externalContactForce.h"

externalContactForce::externalContactForce(elasticPlate &m_plate, timeStepper &m_stepper, 
        double m_stiffness, double m_dBar)
{
	plate = &m_plate;
	stepper = &m_stepper;

    stiffness = m_stiffness;
    dBar = m_dBar;

    ForceVec = VectorXd::Zero(plate->ndof);

    Id3<<1,0,
         0,1;

    IdG<<1,0,
    	 0,0;

    upperBound = plate->radius + dBar;
}

externalContactForce::~externalContactForce()
{
	;
}

void externalContactForce::computeFc()
{
	ForceVec = VectorXd::Zero(plate->ndof);

	totalContactForce = 0.0;

	for(int i = 0; i < plate->nv; i++)
	{
		Vector2d xCurrent = plate->getVertex(i);

		double d = upperBound - xCurrent(1);

		if (d <= dBar)
		{
			//cout << "height :" << d << endl;

			dEnergydD = - 2 * (d - dBar) * log(d / dBar) - (d - dBar) * (d - dBar) / d;

			stepper->addForce(2 * i + 1, - stiffness * dEnergydD);

			totalContactForce = totalContactForce + stiffness * dEnergydD;

			d2EnergydD2 = - 2 * log(d / dBar) - 2 * (d - dBar) / d - 2 * (d - dBar) / d + (d - dBar) * (d - dBar) / (d * d);

			stepper->addJacobian(2 * i + 1, 2 * i + 1, stiffness * d2EnergydD2);
		}

	}
}

void externalContactForce::computeJc()
{
	;
}
