#include "Grid.h"
#include <ostream>

#define ArrayBorder_ij Array[interactionTargets[j]].Borders[interactionTargets[j]]

Grid PhiStabilization(Grid X){
    Grid next = X.copy();
    //Iteration
    #pragma omp parallel for
    for(int i = 0; i < X.GetSize(); i++){
        double Rho = X.Array[i].Data.rho;
        double ci = X.Array[i].Data.c;

        //Electric Potential (Phi) iteration start
        std::vector<unsigned int> interactionTargets = X.Array[i].NeigbourArray;
        double PhiChangeNumerator = 0;
        double PhiChangeDenumerator = 0;
        for(int j = 0; j < interactionTargets.size(); j++){
                double kij = X.Array[i].Data.k;

                double S = L / sqrt3;
                double Distance = L;
                double V = sqrt3 * L / 2;

                PhiChangeNumerator += X.ArrayBorder_ij.sigma * S / Distance * X.Array[interactionTargets[j]].Data.Phi;
                PhiChangeDenumerator += X.ArrayBorder_ij.sigma * S / Distance;
            }
        next.Array[i].Data.Phi = PhiChangeNumerator/PhiChangeDenumerator;
    }
    //std::cout<<X.differencePhi(next)<<"Difference"<<std::endl;
    return next;
}

Grid NextStepForTemp(Grid current, double dt = 0){
    Grid next = current.copy();
    
    double dtactmin = dt;

    //Iteration
    #pragma omp parallel for
    for(int i = 0; i < current.GetSize(); i++){
        double Rho = current.Array[i].Data.rho;
        double ci = current.Array[i].Data.c;

        //Temperature (T) iteration start
        std::vector<unsigned int> interactionTargets = next.Array[i].NeigbourArray;
        const double T = current.Array[i].Data.T;
        double TemperaturechangeFromNeighbours = 0;

        next.Array[i].Data.Q = 0;
        for(int j = 0; j < interactionTargets.size(); j++){
            mathVector e = (current.Array[interactionTargets[j]].r - current.Array[i].r).normalise();
            double kij = current.Array[i].Data.k;

            double Distance = (current.Array[interactionTargets[j]].r - current.Array[i].r)();
            double S = L / sqrt3;
            double V = sqrt3 * L / 2;

            TemperaturechangeFromNeighbours += kij/Rho/ci * (current.Array[interactionTargets[j]].Data.T - T) * S / Distance / V ;

            //Current change
            next.ArrayBorder_ij.E = e * (current.Array[interactionTargets[j]].Data.Phi - current.Array[i].Data.Phi)/Distance;
            next.ArrayBorder_ij.I = next.ArrayBorder_ij.E * next.ArrayBorder_ij.sigma;
            next.ArrayBorder_ij.Q = next.ArrayBorder_ij.I * next.ArrayBorder_ij.E;
            
            next.Array[i].Data.Q += next.ArrayBorder_ij.Q/2;
            //std::cout<<"For i = "<<i<<", neigbour "<< interactionTargets[j] << " Q = "<< next.Array[i].Data.Q/ Rho / ci <<std::endl;
        }
        next.Array[i].Data.T += current.Array[i].Data.Q * dt / Rho / ci + TemperaturechangeFromNeighbours* dt;
        //Temperature (T)  iteration end

        //CFL
        double kij = current.Array[i].Data.k;

        double S = L / sqrt3;
        double Distance = L;
        double V = sqrt3 * L / 2;

        double dtact = 1. / (kij/Rho/ci * S / Distance / V) / double(interactionTargets.size());

        if(dtact < dt) dtactmin = fmin(dtact,dtactmin);
    }
    if(dtactmin < dt )std::cout<< "Use recommended: " <<dtactmin<<std::endl;
    return next;
}

inline double ForceLine(double inL, Node a)
{
    if(inL<L)
    {
        //ветка отталкивания от 0 до L0
        return a.Borders[0].EYoung*L*(inL-L);//repulsiveK*L0*(inL/L0-1); repulsiveK*L0=repulsiveE
    }
    else
    {
        //ветка притяжения от L0 до LD
        if(inL<L)
        {
            if(inL<LP)
            {
                //участок роста от L0 до LP
                return a.Borders[0].EYoung*L*(inL-L);//attractiveK*L0*(inL/L0-1); attractiveK*L0=attractiveE
            }
            else
            {
                //участок спада от LP до LD
                // Плавный спад (кубический)
                //                               double t = (inL - LP) / (LD_T - LP);
                //                               return attractiveK * (LP - L0) * (1 - 3*t*t + 2*t*t*t); ПОДУМАТЬ!!!

                return a.Borders[0].EYoung*(LP-L)*(1-(inL-LP)/(LD-LP));//return attractiveK*L0*(LP/L0-1)*(1-(inL/LP-1)/(LD/LP-1)) attractiveK*L0=attractiveE
            }
        }
        else
        {
            return 0;
        }
    }
    return 0;
}

Grid Movement(Grid current, double dt){
    Grid next = current.copy();
    
    //#pragma omp parallel for
    for(int i = 0; i < current.GetSize(); i++){
        std::vector<unsigned int> interactionTargets = next.Array[i].NeigbourArray;
        
        mathVector Force = Zero;
        for(auto NeighbourAddress : interactionTargets){
            Node interTarget = current.Array[NeighbourAddress];
            mathVector e = (interTarget.r - current.Array[i].r).normalise();
            
            double Distance = (interTarget.r - current.Array[i].r)();
            double S = Distance / sqrt3;

            if(L != 0) Force = Force + e*(-ForceLine(Distance, current.Array[i])/L);
        }
        double m = 1.0;
        mathVector a = Force/m;
        //std::cout<<"Force of "<< current.Array[i].r.y <<" = "<<Force.y<<std::endl;
        next.Array[i].v  = current.Array[i].v + (a-current.Array[i].v*ALPHA)*dt;
        //std::cout<<"Speed "<<next.Array[i].v.y<<std::endl;
        next.Array[i].r  = current.Array[i].r + current.Array[i].v*dt;
        //std::cout<<"Speed "<<next.Array[i].r.y<<std::endl;
    }   
    return next;
}