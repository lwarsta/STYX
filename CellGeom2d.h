#ifndef _CELLGEOM2D_H
#define _CELLGEOM2D_H
#include "CellGeomBase.h"

class CellGeom2d : public CellGeomBase
{
protected:
	double area;
	std::vector<int> neighIndices;
	std::vector<CellGeom2d*> neighPointers;
	std::vector<std::vector<int> > sideVertIndices;
	std::vector<Vertex> sideCentrePoints;
	std::vector<double> sideLengths;
public:
	CellGeom2d(int idNew, char geomTypeNew, std::vector<int> vertIndicesNew, 
		std::vector<Vertex*> vertPointersNew);
	void setNeighIndices(std::vector<int> neighIndicesNew);
	void setNeighPointers(std::vector<CellGeom2d*> neighPointersNew);
	size_t getNumOfNeigh() { return neighPointers.size(); }
	void assignSideVertIndices();
	void computeSideCentrePoints();
	void computeDistancesBetweenCells();
	void computeSideLengths();
	void computeArea();
	double getArea(){return area;}
	std::vector< std::vector<int> > getSideVertIndices() { 
		return sideVertIndices; }
	std::vector<int> getNeighIndices() { return neighIndices; }
	double getSideLength(int index) { return sideLengths.at(index); }
};

#endif
