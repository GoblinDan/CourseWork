#include <cmath>
#include <iostream>
#include <vector>

#include "Simulation.h"

#define OutStep 1
#define Ustart 1.0
// todo: дырка в центре!
// Подсчёт потенциалов, токов и нагрев с sigma = 6e7

/*      Length, x, i
        <-------------->
w, y   ^ --------------
i, j   | |            |
d      | |            |
t      | |            |
h      v --------------
*/



std::vector<unsigned int> BorderLeft;
std::vector<unsigned int> BorderRight;
std::vector<unsigned int> BorderUp;
std::vector<unsigned int> BorderDown;

void outputInFile(Grid target, int num = 0){
    FILE *outFile;
    char fileName[512];
    sprintf(fileName,"CSV/Snapshot%04i.csv",num);
    outFile=fopen(fileName,"w+");
    fprintf(outFile,"x;y;z;T;phi;Q;Ex;Ey;Ez\n");

    for(int i = 0; i < target.GetSize(); i++){
        fprintf(outFile,"%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;\n",
            target.Array[i].r.x,target.Array[i].r.y,target.Array[i].r.z,
            target.Array[i].Data.T,target.Array[i].Data.Phi,target.Array[i].Data.Q,
            0, 0, 0);
    }
}

int main(){
    int Size = 100;
    std::cout<<"Beginning Birth!"<<std::endl;
    Grid test(Size*1.16, Size);

    //Initial state
    for(int i = 0; i < test.GetSize(); i++){
        if(test.Array[i].r.x < interactionRange) BorderLeft.push_back(i);
        if(test.Array[i].r.x > (Size*L - 2*interactionRange)) BorderRight.push_back(i);
        if(test.Array[i].r.y < interactionRange) BorderUp.push_back(i);
        if(test.Array[i].r.y > (Size*L*1.16 - 2*interactionRange)) BorderDown.push_back(i);
    }

    //Electric potential approximation
    for(int i = 0; i < test.GetSize(); i++){
        test.Array[i].Data.Phi = test.Array[i].r.x/(0.01) * (Ustart/2) - Ustart/2;;
    }
    
    std::cout<<"Finished step: "<<0<<std::endl;
    outputInFile(test, 0);
    std::cout<<"Finished outputing step: "<<0<<std::endl;


    Grid next = test.copy();
    int iter = 0;
    do{
        test = next.copy();
        for(int i = 0; i < BorderLeft.size(); i++){
            test.Array[BorderLeft[i]].Data.Phi = -Ustart/2;
        }
        for(int i = 0; i < BorderRight.size(); i++){
            test.Array[BorderRight[i]].Data.Phi = Ustart/2;
        }
        next = PhiStabilization(test);
        iter++;
        if(iter%100 == 0) std::cout<<std::abs(test.differencePhi(next))<<std::endl;
    //} while(std::abs(test.differencePhi(next)) > 6e-7);//;
    } while(iter < 1e3);
    std::cout<<"Finished stabilizing Phi"<<std::endl;

    mathVector vUp = mathVector(0,1e-3);
    mathVector vDown = mathVector(0,-1e-3);

    for (int i = 1; i < 21; i++){
        test = NextStepForTemp(test, 1e-2);
        for(int j = 0; j < 10000; j++) {
            test = Movement(test, 1e-6);
        //std::cout<<"Finished step: "<<i<<std::endl;
            
            for(auto i : BorderUp){
                test.Array[i].v = vUp;
            }
            for(auto i: BorderDown){
                test.Array[i].v = vDown;
            }
        }
        if(i%OutStep == 0)outputInFile(test, i);
        if(i%20 == 0)std::cout<<"Finished step: "<<i<<std::endl;
    }
    return 0;
}