#include "elasticPlate.h"

elasticPlate::elasticPlate(double m_YoungM, double m_density, double m_thickness, 
		double m_Possion, double m_dt, int m_nv, double m_radius)
{
	YoungM = m_YoungM;
	density = m_density;
	thickness = m_thickness;
	Possion = m_Possion;
	dt = m_dt;
	nv = m_nv;
	radius = m_radius;

	setupGeometry();

	ndof = 2 * nv;
	x = VectorXd::Zero(ndof);
	x0 = VectorXd::Zero(ndof);
	u = VectorXd::Zero(ndof);

	for (int i = 0; i < nv; i++)
	{
		x(2 * i + 0) = v_nodes[i](0);
		x(2 * i + 1) = v_nodes[i](1);
	}
	x0 = x;
	x00 = x;

	computeEdge();

	setupMass();

	//set up constraint map
	isConstrained = new int[ndof];
    for (int i=0; i < ndof; i++)
    {
		isConstrained[i] = 0;
    }
}

elasticPlate::~elasticPlate()
{
	delete isConstrained;
	delete unconstrainedMap;
	delete fullToUnconsMap;
}

void elasticPlate::setup()
{
	ncons = 0;
    for (int i=0; i < ndof; i++)
    {
		if (isConstrained[i] > 0)
		{
			ncons++;
		}
	}
	uncons = ndof - ncons;

	unconstrainedMap = new int[uncons]; // maps xUncons to x
	fullToUnconsMap = new int[ndof];
	setupMap();
}

void elasticPlate::setupMap()
{
	int c = 0;
	for (int i=0; i < ndof; i++)
	{
		if (isConstrained[i] == 0)
		{
			unconstrainedMap[c] = i;
			fullToUnconsMap[i] = c;
			c++;
		}
	}
}

void elasticPlate::setupMass()
{
	massArray = VectorXd::Zero(ndof);

	double deltaMass;

	int index1;
	int index2;

	double localRadius;

	for (int i = 0; i < edgeNum; i++)
	{
		localRadius = (v_edgeElement[i].x_1_start(0) + v_edgeElement[i].x_2_start(0)) / 2;

		deltaMass = 2 * M_PI * localRadius * v_edgeElement[i].refLength * thickness * density / 2;

		index1 = v_edgeElement[i].nv_1;
		index2 = v_edgeElement[i].nv_2;

		massArray(2 * index1 + 0) = massArray(2 * index1 + 0) + deltaMass;
		massArray(2 * index1 + 1) = massArray(2 * index1 + 1) + deltaMass;
	
		massArray(2 * index2 + 0) = massArray(2 * index2 + 0) + deltaMass;
		massArray(2 * index2 + 1) = massArray(2 * index2 + 1) + deltaMass;
	}
}

int elasticPlate::getIfConstrained(int k)
{
	return isConstrained[k];
}

void elasticPlate::setVertexBoundaryCondition(Vector2d position, int k)
{
	isConstrained[2 * k + 0] = 1;
	isConstrained[2 * k + 1] = 1;
	
	// Store in the constrained dof vector
	x(2 * k + 0) = position(0);
	x(2 * k + 1) = position(1);
}

void elasticPlate::setOneBoundaryCondition(double position, int i, int k)
{
	isConstrained[2 * i + k] = 1;

	x(2 * i + k) = position;
}

Vector2d elasticPlate::getVertex(int i)
{
	Vector2d xCurrent;

	xCurrent(0) = x(2 * i + 0);
	xCurrent(1) = x(2 * i + 1);
	
	return xCurrent;
}

Vector2d elasticPlate::getVertexOld(int i)
{
	Vector2d xCurrent;

	xCurrent(0) = x0(2 * i + 0);
	xCurrent(1) = x0(2 * i + 1);

	return xCurrent;
}

Vector2d elasticPlate::getVertexInitial(int i)
{
	Vector2d xCurrent;

	xCurrent(0) = x00(2 * i + 0);
	xCurrent(1) = x00(2 * i + 1);

	return xCurrent;
}

Vector2d elasticPlate::getVelocity(int i)
{
	Vector2d uCurrent;

	uCurrent(0) = ( x(2 * i + 0) - x0(2 * i + 0) ) / dt;
	uCurrent(1) = ( x(2 * i + 1) - x0(2 * i + 1) ) / dt;
	
	return uCurrent;
}


Vector2d elasticPlate::getVelocityOld(int i)
{
	Vector2d uCurrent;

	uCurrent(0) = u(2 * i + 0) ;
	uCurrent(1) = u(2 * i + 1) ;
	
	return uCurrent;
}

void elasticPlate::updateEdgePair()
{
	for (int i = 0; i < edgeNum; i++)
	{
		v_edgeElement[i].x_1 = getVertex(v_edgeElement[i].nv_1);
		v_edgeElement[i].x_2 = getVertex(v_edgeElement[i].nv_2);
		v_edgeElement[i].edgeLength = (v_edgeElement[i].x_1 - v_edgeElement[i].x_2).norm();
	}
}

void elasticPlate::prepareForIteration()
{
	updateEdgePair();
}

void elasticPlate::updateTimeStep()
{
	prepareForIteration();

	// compute velocity
	u = (x - x0) / dt;

	// update x
	x0 = x;
}

void elasticPlate::updateGuess()
{
	for (int c=0; c < uncons; c++)
	{
		x[unconstrainedMap[c]] = x[unconstrainedMap[c]] + u[unconstrainedMap[c]] * dt;
	}
}

void elasticPlate::updateNewtonMethod(VectorXd m_motion)
{
	for (int c=0; c < uncons; c++)
	{
		x[unconstrainedMap[c]] -= m_motion[c];
	}
}

void elasticPlate::setupGeometry()
{
	v_nodes.clear();
    edge.clear();

    double thetaSpace = (M_PI / 2.00001);

    double deltaTheta = thetaSpace / ( nv - 1);

    double rLocal = radius;

    for (int i = 0; i < nv; i++)
    {
    	double localTheta = deltaTheta * i;

    	Vector2d xCurrent;

    	xCurrent(0) = rLocal * cos(localTheta);
    	xCurrent(1) = rLocal * sin(localTheta);

    	v_nodes.push_back(xCurrent);
    }
    nv = v_nodes.size();
    

    for (int i = 0; i < nv - 1; i++)
    {
    	Vector2i edgeCurrent;

    	edgeCurrent(0) = i;
    	edgeCurrent(1) = i + 1;

    	edge.push_back(edgeCurrent);
    }

}

void elasticPlate::computeEdge()
{
	edgeNum = 0;
	v_edgeElement.clear();

	for (int i = 0; i < edge.size(); i++)
	{
		Vector2i edgeCurrent = edge[i];

		edgeElement m_edgeElement;

		m_edgeElement.nv_1 = edgeCurrent(0);
		m_edgeElement.nv_2 = edgeCurrent(1);

		m_edgeElement.x_1_start = getVertex(m_edgeElement.nv_1);
		m_edgeElement.x_2_start = getVertex(m_edgeElement.nv_2);

		m_edgeElement.x_1 = getVertex(m_edgeElement.nv_1);
		m_edgeElement.x_2 = getVertex(m_edgeElement.nv_2);

		m_edgeElement.refLength = (m_edgeElement.x_2- m_edgeElement.x_1).norm();
		m_edgeElement.edgeLength = m_edgeElement.refLength;

		Vector2d tangent = (m_edgeElement.x_2 - m_edgeElement.x_1) / (m_edgeElement.x_2 - m_edgeElement.x_1).norm();

		m_edgeElement.phiBar = atan( tangent(1) / tangent(0) );

		m_edgeElement.arrayNum = VectorXi::Zero(4);

		m_edgeElement.arrayNum(0) = 2 * m_edgeElement.nv_1 + 0;
		m_edgeElement.arrayNum(1) = 2 * m_edgeElement.nv_1 + 1;
		
		m_edgeElement.arrayNum(2) = 2 * m_edgeElement.nv_2 + 0;
		m_edgeElement.arrayNum(3) = 2 * m_edgeElement.nv_2 + 1;
		
		v_edgeElement.push_back(m_edgeElement);

		edgeNum = edgeNum + 1;
	}
}