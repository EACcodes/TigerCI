/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/*
 * File:   Cylinder.cpp
 * Author: Francis Ricci
 *
 */

#include <sstream>
#include <string>
#include <iostream>
#include <tiger_cpp/Cylinder.h>
#include <tiger_cpp/OrbitalParsers.h>

// create cylinder
Cylinder::Cylinder(int a, int b){
    _a = a;
    _b = b;
}

// print cylinder
std::string Cylinder::toString(){
    ostringstream oss;
    oss << _a << " " << _b;
    return oss.str();
}

//return first endpoint
int Cylinder::firstEndpoint(){
    return _a;
}

//return second endpoint
int Cylinder::secondEndpoint(){
    return _b;
}

// creates an empty chain
Chain::Chain(){
    Chain(-1);
}

// creates new empty chain
Chain::Chain(int orb){
    _num_cylinders = 0;
    _orb_num = orb;
}

// add a cylinder to the chain
void Chain::addCylinder(Cylinder cyl){
    _cylinders.push_back(cyl);
    _num_cylinders++;
}

// get cylinders in the chain
void Chain::print(){
    bool first = true;
    std::cout << _orb_num << ": ";
    for (Cylinder cyl : _cylinders){
        if (!first)
            std::cout << ", ";
        first = false;
        std::cout << cyl.toString();
    }
    std::cout << std::endl;
}

// return number of cylinders in the chain
int Chain::numCylinders(){
    return _num_cylinders;
}

// return the chain's orbital
int Chain::getOrb(){
    return _orb_num;
}
// return cylinders in the chain
std::vector<Cylinder> Chain::getCylinders(){
    return _cylinders;
}
