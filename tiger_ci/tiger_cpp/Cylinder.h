/*
 * File:   Cylinder.h
 * Author: Francis Ricci
 *
 */

#include <vector>
#include <string>

#ifndef CYLINDER_H
#define	CYLINDER_H

class Cylinder{
public:
    // create cylinder
    Cylinder(int a, int b);
    
    // return string representation of cylinder
    std::string toString();
    
    //return first endpoint
    int firstEndpoint();
    
    //return second endpoint
    int secondEndpoint();
    
private:
    // endpoint values;
    int _a;
    int _b;
};

class Chain{
public:
    // creates an empty chain
    Chain();
    
    // creates new empty chain on a given orbital
    Chain(int orb);
    
    // add a cylinder to the chain
    void addCylinder(Cylinder cyl);
    
    // return number of cylinders in the chain
    int numCylinders();
    
    // return cylinders in the chain
    std::vector<Cylinder> getCylinders();
    
    // return the chain's orbital
    int getOrb();
    
    // get cylinders in the chain
    void print();
    
private:
    // set of cylinders in the chain
    std::vector<Cylinder> _cylinders;
    int _num_cylinders;
    int _orb_num;
};

#endif	/* CYLINDER_H */


