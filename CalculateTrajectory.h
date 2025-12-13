#ifndef CALCULATETRAJECTORY_H
#define CALCULATETRAJECTORY_H

#include "Grid.h"

//явный метод Эйлера
//в терминах скорость-ускорение
void explicitEulerMethodVA(std::vector<Node> &inNodes,double dt);

void explicitEulerMethodVA(std::vector<Node> &inNodes, double dt)
{
    const mathVector zero=Zero;
    for(size_t i=0; i<inNodes.size(); ++i)
    {
        Node &n=inNodes[i];
        n.v+=(n.a-n.v*ALPHA)*dt;
        n.r+=n.v*dt;
        n.a=zero;
    }
}
#endif // CALCULATETRAJECTORY_H
