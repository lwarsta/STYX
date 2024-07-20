#ifndef _LINKJUNCBASE_H
#define _LINKJUNCBASE_H
#include <iostream>
#include <vector>
#include <string>
#include "Algorithms.h"
#include "CellBase.h"
#include "Vertex.h"

class LinkJuncGeomBase : public CellBase
{
protected:
	//int id;
	//int material;
	//int init_cond_ind;
	//int bound_cond_ind;
	int geomType;
	Vertex centrePoint;
	std::vector<int> vertIndices;
	std::vector<Vertex*> vertPointers;
public:
	LinkJuncGeomBase();
	//int get_id() {return id;}
	//void set_id(int id_new) {id = id_new;}
	//int get_material() {return material;}
	//void set_material(int material_new) {material = material_new;}
	//int get_init_cond_ind() {return init_cond_ind;}
	//void set_init_cond_ind(int init_cond_ind_new) {init_cond_ind = init_cond_ind_new;}
	//int get_bound_cond_ind() {return bound_cond_ind; }
	//void set_bound_cond_ind(int bound_cond_ind_new) {bound_cond_ind = bound_cond_ind_new;}
	LinkJuncGeomBase(int idNew, char geomTypeNew, std::vector<int> vertIndicesNew,
		std::vector<Vertex*> vertPointersNew);
	int getGeomType() { return geomType; }
	int setGeomType(int geomTypeNew);
	std::vector<int> getVertIndices() { return vertIndices; }
	std::vector<Vertex*> getVertPointers() { return vertPointers; }
	void compCenterPoint();
	Vertex getCentrePoint() { return centrePoint; }
};

#endif
