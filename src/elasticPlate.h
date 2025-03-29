#ifndef ELASTICPLATE_H
#define ELASTICPLATE_H

#include "eigenIncludes.h"
#include <fstream>

struct edgeElement
{
	int nv_1;
	int nv_2;

	Vector2d x_1;
	Vector2d x_2;

	double phiBar;

	Vector2d x_1_start;
	Vector2d x_2_start;

	double refLength;
	double edgeLength;

	VectorXi arrayNum;
};

class elasticPlate
{
	public:
	elasticPlate(double m_YoungM, double m_density, double m_thickness, 
		double m_Possion, double m_dt, int m_nv, double m_radius);
	~elasticPlate();

	double YoungM;
	double thickness;
	double Possion;
	double dt;
	double density;
	double radius;

	Vector2d getVertex(int i);
	Vector2d getVertexOld(int i);
	Vector2d getVertexInitial(int i);
	Vector2d getVelocity(int i);
	Vector2d getVelocityOld(int i);

	VectorXd x;
	VectorXd x0;
	VectorXd x00;
	VectorXd u;

	std::vector<Vector2d> v_nodes;
    std::vector<Vector2i> edge;

	std::vector<edgeElement> v_edgeElement;
	
	int temp;

	int nv;
	int edgeNum;

	int ndof;
	int uncons;
	int ncons;

	void setupGeometry();

	void setVertexBoundaryCondition(Vector2d position, int k);
	void setOneBoundaryCondition(double position, int i, int k);

	void computeEdge();
	
	void updateEdgePair();
	
	// boundary conditions
	int* isConstrained;
	int getIfConstrained(int k);
	int* unconstrainedMap;
	int* fullToUnconsMap;
	void setup();
	void setupMap();

	void updateTimeStep();
	void updateGuess();
	void updateNewtonMethod(VectorXd m_motion);
	void prepareForIteration();

	VectorXd massArray;
	void setupMass();

	private:
};

#endif
