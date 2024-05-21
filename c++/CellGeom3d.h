#ifndef _CELLGEOM3D_H
#define _CELLGEOM3D_H
#include "CellGeomBase.h"

class CellGeom3d : public CellGeomBase
{
protected:
	double volume;
	std::vector<int> neighIndices;
	std::vector<CellGeom3d*> neighPointers;
	std::vector< std::vector<int> > faceVertIndices;
	std::vector<Vertex> faceCentrePoints;
	std::vector<double> faceAreas;
public:
	CellGeom3d(int idNew, char geomTypeNew, std::vector<int> vertIndicesNew, std::vector<Vertex*> vertPointersNew);
	void setNeighIndices(std::vector<int> neighIndicesNew);
	//std::vector<int> * getNeighIndices(){return &neighIndices;} // NOT USED
	void setNeighPointers(std::vector<CellGeom3d*> neighPointersNew);
	size_t getNumOfNeigh(){return neighPointers.size();}
	void assignFaceVertIndices();
	void computeFaceCentrePoints();
	void computeDistancesBetweenCells();
	void computeFaceAreas();
	void computeVolume();
	std::vector<std::vector<int> > getFaceVertIndices() { return faceVertIndices; }
	std::vector< Vertex > getFaceCentrePoints() { return faceCentrePoints; }
	std::vector<int> getNeighIndices() { return neighIndices; }
	double getFaceArea(int index){return faceAreas.at(index);}
	double getVolume(){return volume;}
};

#endif
