#include <cmath>
#include <unordered_map>
#include <vector>
#include <iostream>

#define Zero mathVector()

struct mathVector{
    double x, y, z;
    mathVector(){x = 0; y = 0; z = 0;}
    mathVector(double _x, double _y, double _z){x = _x; y = _y; z = _z;}
    mathVector(double _x, double _y){x = _x; y = _y; z = 0;}
    mathVector operator-(const mathVector& b){
        return mathVector(x - b.x, y - b.y, z - b.z);
    }
    mathVector operator+(const mathVector& b){
        return mathVector(x + b.x, y + b.y, z + b.z);
    }
    mathVector operator*(const double d){
        return mathVector(x*d, y*d, z*d);
    }
    double operator*(const mathVector d){
        return x*d.x + y*d.y + z*d.z;
    }
    mathVector operator/(const double d){
        return mathVector(x/d, y/d, z/d);
    }

    mathVector& operator+=(const mathVector& rhs){
        this->x += rhs.x;
        this->y += rhs.y;
        this->z += rhs.z;
        return *this;
    }
    bool operator<(const double l){
        return (x*x + y*y + z*z) < (l*l);
    }
    bool operator>(const double l){
        return (x*x + y*y + z*z) > (l*l);
    }
    double operator()(){
        return std::sqrt(x*x + y*y + z*z);
    }
    mathVector normalise(){
        mathVector me = mathVector(this->x,this->y,this->z);
        return me / me();
    }
};


struct PhysicsData{
    double Phi, T, rho, k, c, Q, density;

    PhysicsData(){
        Phi = 0;
        Q   = 0; 
        T   = 300;  
        //standart is Cu
        rho = 8960;
        k = 400;
        c = 385;
        density = 8.96; //g*cm^3
    }
    PhysicsData(double _Phi, double _T){
        Phi = _Phi;
        T   = _T;
    }
};

struct BorderData{
    double sigma, Q, EYoung;
    mathVector  I, E;
    BorderData(){
        I = Zero;
        sigma = 6e7; // standart is Cu
        EYoung = 1e9;
        E = Zero;
        Q = 0;
    }
};

struct Node{
    PhysicsData Data;
    mathVector r;
    mathVector v;
    std::vector<unsigned int> NeigbourArray;
    std::unordered_map<unsigned int, BorderData> Borders;
    Node(float _x, float  _y){
        r.x = _x; r.y = _y; r.z = 0;
        v.x = 0; v.y = 0; v.z = 0;
        Data = PhysicsData();
        NeigbourArray = std::vector<unsigned int>(0);
    }
};