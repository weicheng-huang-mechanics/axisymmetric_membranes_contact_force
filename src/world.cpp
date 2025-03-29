#include "world.h"

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				
	saveData = m_inputData.GetBoolOpt("saveData");		

	deltaTime = m_inputData.GetScalarOpt("deltaTime");     
	totalTime = m_inputData.GetScalarOpt("totalTime");    
	YoungM = m_inputData.GetScalarOpt("YoungM");
	density = m_inputData.GetScalarOpt("density");

	thickness = m_inputData.GetScalarOpt("thickness");
	radius = m_inputData.GetScalarOpt("radius");
	Possion = m_inputData.GetScalarOpt("Possion");

	stol = m_inputData.GetScalarOpt("stol");
	forceTol = m_inputData.GetScalarOpt("forceTol");
	scaleRendering = m_inputData.GetScalarOpt("scaleRendering");
	maxIter = m_inputData.GetIntOpt("maxIter");
	gVector = m_inputData.GetVecOpt("gVector");
	viscosity = m_inputData.GetScalarOpt("viscosity");
	nv = m_inputData.GetIntOpt("nv");

	stiffness = m_inputData.GetScalarOpt("stiffness");
	dBar = m_inputData.GetScalarOpt("dBar");

	shearM = YoungM / (2 * (1 + Possion));
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

void world::OpenFile(ofstream &outfile)
{
	if (saveData==false) return;
	
	int systemRet = system("mkdir datafiles"); //make the directory
	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	// Open an input file named after the current time
	ostringstream name;
	name.precision(6);
	name << fixed;

    name << "datafiles/simDER";
    name << ".txt";

	outfile.open(name.str().c_str());
	outfile.precision(10);	
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}

	if ( timeStep % 10 != 0)
	{
		return;
	}

	outfile << m_externalPressureForce->pressure << " " << abs(m_externalContactForce->totalContactForce) << endl;
}

void world::setPlateStepper()
{
	// Create the plate 
	plate = new elasticPlate(YoungM, density, thickness, Possion, deltaTime, nv, radius);

	plateBoundaryCondition();

	plate->setup();

	stepper = new timeStepper(*plate);

	// set up force
	m_inertialForce = new inertialForce(*plate, *stepper);
	m_gravityForce = new externalGravityForce(*plate, *stepper, gVector);
	m_dampingForce = new dampingForce(*plate, *stepper, viscosity);
	m_externalPressureForce = new externalPressureForce(*plate, *stepper);
	m_hyperElasticMM = new hyperElasticMM(*plate, *stepper, 0.4375*shearM, 0.0625*shearM);
	m_externalContactForce = new externalContactForce(*plate, *stepper, stiffness, dBar);
	
	plate->updateTimeStep();

	// set up first jacobian
	m_inertialForce->setFirstJacobian();
	m_dampingForce->setFirstJacobian();
	m_hyperElasticMM->setFirstJacobian();

	stepper->first_time_PARDISO_setup();

	// time step 
	Nstep = totalTime / deltaTime;
	timeStep = 0;
	currentTime = 0.0;

	m_externalPressureForce->pressure = 0.0;
}

void world::plateBoundaryCondition()
{
	Vector2d xStart = plate->getVertex(0);
	plate->setVertexBoundaryCondition(xStart, 0);

	Vector2d xEnd = plate->getVertex(plate->nv - 1);
	plate->setOneBoundaryCondition(xEnd(0), plate->nv - 1, 0);	
}

void world::updateTimeStep()
{
	bool goodSolved = false;

	while (goodSolved == false)
	{
		// Start with a trial solution for our solution x
		plate->updateGuess(); // x = x0 + u * dt

		updateEachStep();

		goodSolved = true;
	}

	plate->updateTimeStep();

	if (render) 
	{
		cout << "time: " << currentTime << endl;
	}

	currentTime += deltaTime;
		
	timeStep++;
}

void world::updateEachStep()
{
	double normf = forceTol * 10.0;
	double normf0 = 0;
	
	bool solved = false;
	
	int iter = 0;

	m_externalPressureForce->pressure = m_externalPressureForce->pressure + 1.0 * deltaTime;

	while (solved == false)
	{
		plate->prepareForIteration();

		stepper->setZero();

		m_inertialForce->computeFi();
		//m_gravityForce->computeFg();
		m_dampingForce->computeFd();
		m_hyperElasticMM->computeFs();
		m_externalPressureForce->computeFp();
		m_externalContactForce->computeFc();
	
		normf = stepper->GlobalForceVec.norm();

		if (iter == 0) 
		{
			normf0 = normf;
		}
		
		if (normf <= forceTol)
		{
			solved = true;
		}
		else if(iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}

		//cout << normf << endl;

		normf = 0.0;
		
		if (solved == false)
		{
			m_inertialForce->computeJi();
			//m_gravityForce->computeJg();
			m_dampingForce->computeJd();
			m_hyperElasticMM->computeJs();
			m_externalPressureForce->computeJp();;
			
			stepper->integrator(); // Solve equations of motion
			plate->updateNewtonMethod(stepper->GlobalMotionVec); // new q = old q + Delta q
			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			timeStep = Nstep;

			break;
		}
	}

	if (render)
	{
		cout << "iter " << iter << endl;
	}
}

int world::simulationRunning()
{
	if (timeStep < Nstep)
	{
		return 1;
	}
	else 
	{
		return -1;
	}
}

Vector2d world::getScaledCoordinate(int i, int j)
{
	Vector2d xCurrent;
	
	if (j == 0)
	{
		xCurrent = plate->v_edgeElement[i].x_1 * scaleRendering;
	}
	if (j == 1)
	{
		xCurrent = plate->v_edgeElement[i].x_2 * scaleRendering;
	}

	return xCurrent;
}

Vector2d world::getScaledPos(int i)
{
	Vector2d xCurrent;

	xCurrent = plate->getVertex(i) * scaleRendering;

	return xCurrent;
}

int world::getNv()
{
	return plate->nv;
}

int world::numStretchingPair()
{
	return plate->edgeNum;
}

void world::computeReactionForce()
{
	//reactionForce = VectorXd::Zero(plate->ndof);
}

double world::getBound()
{
	return m_externalContactForce->upperBound * scaleRendering;
}