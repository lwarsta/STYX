#include "CellGeomBase.h"

CellGeomBase::CellGeomBase()
{
	id = -1;
	geomType = -1;
	gridConnect = -1;
}

CellGeomBase::CellGeomBase(int idNew, char geomTypeNew, 
	std::vector<int> vertIndicesNew, std::vector<Vertex*> vertPointersNew)
{
	id = idNew;
	geomType = geomTypeNew;
	vertIndices = vertIndicesNew;
	vertPointers = vertPointersNew;
	//computeCentrePoint();
}
/*
void CellGeomBase::setNeighIndices(std::vector<int> neighIndicesNew)
{
	neighIndices = neighIndicesNew;
}
*/

int CellGeomBase::setGeomType(int geomTypeNew)
{
	// Check for illegal geometry types and return error.
	geomType = geomTypeNew;
	return 0;
}

void CellGeomBase::computeCentrePoint()
{
	centrePoint.x = 0.0;
	centrePoint.y = 0.0;
	centrePoint.z = 0.0;

	for (size_t i = 0; i < vertPointers.size(); i++) {
		centrePoint.x += vertPointers.at(i)->x;
		centrePoint.y += vertPointers.at(i)->y;
		centrePoint.z += vertPointers.at(i)->z;
	}

	if (vertPointers.size() > 0) {
		centrePoint.x /= vertPointers.size();
		centrePoint.y /= vertPointers.size();
		centrePoint.z /= vertPointers.size();
	}
}

Vertex CellGeomBase::createVector(Vertex point0, Vertex point1)
{
	Vertex vector;
	vector.x = point0.x - point1.x;
	vector.y = point0.y - point1.y;
	vector.z = point0.z - point1.z;

	return vector;
}

double CellGeomBase::computeVectorLength(Vertex vector)
{
	return sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
}

double CellGeomBase::computeVectorDotProduct(Vertex point0, Vertex point1)
{
	return point0.x * point1.x + point0.y * point1.y + point0.z * point1.z;
}

Vertex CellGeomBase::addVectors(Vertex point0, Vertex point1)
{
	point0.x = point0.x + point1.x;
	point0.y = point0.y + point1.y;
	point0.z = point0.z + point1.z;

	return point0;
}

Vertex CellGeomBase::computeUnitVector(Vertex vector)
{
	double length = computeVectorLength(vector);

	if (length > 0.0)
	{
		vector.x = vector.x / length;
		vector.y = vector.y / length;
		vector.z = vector.z / length;
		return vector;
	}

	// When a zero length vector is sent in, return the same vector.
	else
	{
		return vector;
	}
}

Vertex CellGeomBase::computeVectorCrossProduct(Vertex vector0, 
	Vertex vector1)
{
	Vertex crossProduct;
	crossProduct.x = vector0.y * vector1.z - vector0.z * vector1.y;
	crossProduct.y = vector0.z * vector1.x - vector0.x * vector1.z;
	crossProduct.z = vector0.x * vector1.y - vector0.y * vector1.x;

	return crossProduct;
}
