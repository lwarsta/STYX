#ifndef _LINKJUNCBASE_H
#define _LINKJUNCBASE_H
#include <iostream>
#include <string>

class LinkJuncBase
{
protected:
	int id;
	int material;
	int init_cond_ind;
	int bound_cond_ind;
public:
	LinkJuncBase();
	int get_id() {return id;}
	void set_id(int id_new) {id = id_new;}
	int get_material() {return material;}
	void set_material(int material_new) {material = material_new;}
	int get_init_cond_ind() {return init_cond_ind;}
	void set_init_cond_ind(int init_cond_ind_new) {init_cond_ind = init_cond_ind_new;}
	int get_bound_cond_ind() {return bound_cond_ind; }
	void set_bound_cond_ind(int bound_cond_ind_new) {bound_cond_ind = bound_cond_ind_new;}
};

#endif
