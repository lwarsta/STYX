#include "LinkJuncGeomBase.h"

LinkJuncGeomBase::LinkJuncGeomBase()
{
	//id = 0;
	//material = 0;
	//init_cond_ind = 0;
	//bound_cond_ind = 0;
	id = -1;
	geomType = -1;
}

LinkJuncGeomBase::LinkJuncGeomBase(int idNew, char geomTypeNew,
	std::vector<int> vertIndicesNew, std::vector<Vertex*> vertPointersNew)
{
	id = idNew;
	geomType = geomTypeNew;
	vertIndices = vertIndicesNew;
	vertPointers = vertPointersNew;
	//computeCentrePoint();
}

int LinkJuncGeomBase::setGeomType(int geomTypeNew)
{
	// Check for illegal geometry types and return error.
	geomType = geomTypeNew;
	return 0;
}

void  LinkJuncGeomBase::compCenterPoint()
{
	Algorithms algorithms;
	std::vector<Vertex> points;

	for (size_t i = 0; i < vertPointers.size(); i++) {
		points.push_back( *vertPointers.at(i));
	}

	centrePoint = algorithms.compute_centre_point(points);
}