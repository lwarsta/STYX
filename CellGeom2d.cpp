#include "CellGeom2d.h"

// Superclass constructor is also called below.
CellGeom2d::CellGeom2d(int idNew, char geomTypeNew, 
	std::vector<int> vertIndicesNew, std::vector<Vertex*> vertPointersNew)
	: CellGeomBase(idNew, geomTypeNew, vertIndicesNew, vertPointersNew)
{
	// WHY ARE THESE COMMENTED OUT?
	//assignSideVertIndices();
	//computeSideCentrePoints();
	//computeSideLengths();
	//computeArea();
}

void CellGeom2d::setNeighIndices(std::vector<int> neighIndicesNew)
{
	neighIndices = neighIndicesNew;
}

void CellGeom2d::setNeighPointers(std::vector<CellGeom2d*> neighPointersNew)
{
	neighPointers = neighPointersNew;
	//computeDistancesBetweenCells();
}

void CellGeom2d::assignSideVertIndices()
{
	// Triangle.
	if (geomType == 5)
	{
		std::vector<int> side_01{ 2,0 };
		sideVertIndices.push_back(side_01);
		std::vector<int> side_02{ 0,1 };
		sideVertIndices.push_back(side_02);
		std::vector<int> side_03{ 1,2 };
		sideVertIndices.push_back(side_03);
	}
	// Quadrilateral.
	else if (geomType == 9) {
		std::vector<int> sideBack{ 2,3 }; // 0,1
		sideVertIndices.push_back(sideBack);
		std::vector<int> sideLeft{ 3,0 }; //1,2
		sideVertIndices.push_back(sideLeft);
		std::vector<int> sideFront{ 0,1 }; //2,3
		sideVertIndices.push_back(sideFront);
		std::vector<int> sideRight{ 1,2 }; //3,0
		sideVertIndices.push_back(sideRight);
	}
}

void CellGeom2d::computeSideCentrePoints()
{
	for (size_t i = 0; i < sideVertIndices.size(); i++) {
		Vertex sideCentrePoint;
		sideCentrePoint.x = 0.0;
		sideCentrePoint.y = 0.0;
		sideCentrePoint.z = 0.0;
		for (size_t j = 0; j < sideVertIndices.at(i).size(); j++) {
			sideCentrePoint.x += vertPointers.at(sideVertIndices.at(i).at(j))->x;
			sideCentrePoint.y += vertPointers.at(sideVertIndices.at(i).at(j))->y;
			sideCentrePoint.z += vertPointers.at(sideVertIndices.at(i).at(j))->z;
		}
		if (sideVertIndices.at(i).size() > 0) {
			sideCentrePoint.x /= sideVertIndices.at(i).size();
			sideCentrePoint.y /= sideVertIndices.at(i).size();
			sideCentrePoint.z /= sideVertIndices.at(i).size();
		}
		sideCentrePoints.push_back(sideCentrePoint);
	}
}

void CellGeom2d::computeDistancesBetweenCells()
{
	for (size_t i = 0; i < sideCentrePoints.size(); i++) {
		if (neighPointers.at(i) != 0) {
			Vertex centrePointNeigh = neighPointers.at(i)->getCentrePoint();
			Vertex vector0 = createVector(centrePoint, sideCentrePoints.at(i));
			double distance0 = computeVectorLength(vector0);
			Vertex vector1 = createVector(sideCentrePoints.at(i), centrePointNeigh);
			double distance1 = computeVectorLength(vector1);
			distances.push_back(distance0 + distance1);
		}
		else {
			distances.push_back(0.0);
		}
	}
}

void CellGeom2d::computeSideLengths()
{
	for (size_t i = 0; i < sideVertIndices.size(); i++) {
		Vertex vector = createVector(*vertPointers.at(sideVertIndices.at(i).at(0)),
			*vertPointers.at(sideVertIndices.at(i).at(1)));
		sideLengths.push_back(computeVectorLength(vector));
	}
}

void CellGeom2d::computeArea()
{
	// Should this be projected area to a flat plane or true area?
	// Triangle.
	if (geomType == 5)
	{
		Vertex vector1 = createVector(*vertPointers.at(0), *vertPointers.at(2));
		Vertex vector2 = createVector(*vertPointers.at(0), *vertPointers.at(1));
		Vertex crossProduct = computeVectorCrossProduct(vector1, vector2);
		area = 0.5 * computeVectorLength(crossProduct);
	}
	// Quadrilateral.
	else if (geomType == 9) {
		Vertex vector1 = createVector(*vertPointers.at(0), *vertPointers.at(3));
		Vertex vector2 = createVector(*vertPointers.at(0), *vertPointers.at(1));
		Vertex crossProduct = computeVectorCrossProduct(vector1, vector2);
		area = 0.5 * computeVectorLength(crossProduct);
		vector1 = createVector(*vertPointers.at(2), *vertPointers.at(1));
		vector2 = createVector(*vertPointers.at(2), *vertPointers.at(3));
		crossProduct = computeVectorCrossProduct(vector1, vector2);
		area += 0.5 * computeVectorLength(crossProduct);
	}
}
