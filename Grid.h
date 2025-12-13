#include "Node.h"
#include <cmath>
#include <omp.h>

const double L = 0.1e-3; //Distance between nodes, twice the radius
const double interactionRange = L*1.001;

#define sqrt3 1.73205
#define sqrt2 1.41421

class Grid{
    private:
        unsigned int Width = 0; unsigned int Length = 0;
    public:
        std::vector<Node> Array;

        Grid copy(){
            Grid NewOne = Grid(Width,Length,false);
            NewOne.Array = Array;
            return NewOne;
        }

        Grid(unsigned int _Width, unsigned int _Length, bool fill = true){
            Width = _Width; Length = _Length;
             //create array

            if(fill){
                for(int i = 0; i < Length; i++){
                    for(int j = 0; j < Width; j++){
                        double x = i * L + (0.5 * L * double((j % 2)==0));// x is just two radiuses
                        double y = j * L * sqrt3 / 2;
                        if((mathVector(x,y) - mathVector(Length*L/2,Width*L/2/1.16)) > Width*L/4)
                            Array.push_back(Node(x, y));
                    }
                }
                std::cout<<"Neigbours"<<std::endl;
                //Check Neighbours
                for(int i = 0; i < Array.size(); i++){
                    Array[i].NeigbourArray = getNeigbours(i);
                    for(int j = 0; j < Array[i].NeigbourArray.size(); j++){
                        Array[i].Borders[Array[i].NeigbourArray[j]] = BorderData();
                    }
                }
            }
        }

        std::vector<unsigned int> getNeigbours(unsigned int Index){
            std::vector<unsigned int> output;
            mathVector coords = Array[Index].r;
            
            for(int i = 0; i < Array.size(); i++){
                    unsigned int ind = i;
                    if((coords - Array[ind].r) < interactionRange && ind != Index)
                        output.push_back(ind);
            }

            return output;
        }

        unsigned int GetSize(){return Array.size();}
        double differencePhi(Grid prev){
            double result = 0;
            for(int i = 0; i < Array.size(); i++){
                result += (Array[i].Data.Phi - prev.Array[i].Data.Phi)/Width/Length;
            }
            return result;
        }
};
