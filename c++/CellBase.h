#ifndef _CELLBASE_H
#define _CELLBASE_H
#include <iostream>
#include <string>

class CellBase
{
protected:
	int id;
	int material;
	int initCondInd;
	int boundCondInd;
	int gridConnect;
public:
    CellBase();
	int getId() {return id;}
	void setId(int idNew) {id = idNew;}
	int getMaterial() {return material;}
	void setMaterial(int materialNew) {material = materialNew;}
	int getInitCondInd() {return initCondInd;}
	void setInitCondInd(int initCondIndNew) {initCondInd = initCondIndNew;}
	int getBoundCondInd() {return boundCondInd; }
	void setBoundCondInd(int boundCondIndNew) {boundCondInd = boundCondIndNew;}
	void setGridConnection(int gridConnectNew) { gridConnect = gridConnectNew; }
	int getGridConnection() { return gridConnect; }
};

#endif
