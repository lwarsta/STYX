#pragma once
class Vertex
{
private:
    int id;
public:
    Vertex();
    Vertex(int idNew, double xNew, double yNew, double zNew);
    double x;
    double y;
    double z;
};