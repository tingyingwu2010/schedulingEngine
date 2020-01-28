#ifndef _INDIVIDUAL_
#define _INDIVIDUAL_

#include<vector>
#include<iostream>

class Individual{
	public : 
		double fitness; // Fitness of the chromosome
		unsigned int fctObj; // Objective function
		std::vector<unsigned int> jobs; // Jobs list
		std::vector<unsigned int> machines; // Machines list

	public :
		Individual(int n); //Constructor 
		Individual(unsigned int obj, int n, std::vector<std::vector<int> > corresp); // Constructor with initialization
};

bool operator<(Individual const& i1, Individual const& i2); // Comparator
bool operator>(Individual const& i1, Individual const& i2); // Comparator
std::ostream& operator<<(std::ostream &out, Individual const& i1); // Displayer

#endif
