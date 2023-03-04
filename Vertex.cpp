#include "Vertex.h"

Vertex::Vertex()
{
    id = -1;
    x = 0.0;
    y = 0.0;
    z = 0.0;
}

Vertex::Vertex(int idNew, double xNew, double yNew, double zNew)
{
    id = idNew;
    x = xNew;
    y = yNew;
    z = zNew;
}