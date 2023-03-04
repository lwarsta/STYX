#ifndef _CELLGEOMBASE_H
#define _CELLGEOMBASE_H
#include <vector>
#include <math.h>
#include "CellBase.h"
#include "Vertex.h"

class CellGeomBase : public CellBase
{
protected:
	int geomType;
	Vertex centrePoint;
	std::vector<int> vertIndices;
	std::vector<Vertex*> vertPointers;
	std::vector<double> distances;
	int gridConnect;
public:
	CellGeomBase();
	CellGeomBase(int idNew, char geomTypeNew, std::vector<int> vertIndicesNew, 
		std::vector<Vertex*> vertPointersNew);
	int getGeomType() { return geomType; }
	int setGeomType(int geomTypeNew);
	std::vector<int> getVertIndices() { return vertIndices; }
	//void setNeighIndices(std::vector<int> neighIndicesNew);
	double getDistance(int index){return distances.at(index);}
	Vertex getCentrePoint() { return centrePoint; }
	void computeCentrePoint();
	Vertex createVector(Vertex point0, Vertex point1);
	double computeVectorLength(Vertex vector);
	double computeVectorDotProduct(Vertex point0, Vertex point1);
	Vertex addVectors(Vertex point0, Vertex point1);
	Vertex computeUnitVector(Vertex vector);
	Vertex computeVectorCrossProduct(Vertex point0, Vertex point1);
	void setGridConnection(int gridConnectNew){gridConnect = gridConnectNew;}
	int getGridConnection(){return gridConnect;}
};

#endif
