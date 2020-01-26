#ifndef _INDIVIDUAL_
#define _INDIVIDUAL_

#include<vector>
#include<iostream>

class Individual{
	public : 
		double fitness;
		unsigned int fctObj;	
		std::vector<unsigned int> jobs;
		std::vector<unsigned int> machines;

	public :
		Individual(int n); //Constructor 
		Individual(unsigned int obj, int n, std::vector<std::vector<int> > corresp);
		
};

bool operator<(Individual const& i1, Individual const& i2);
bool operator>(Individual const& i1, Individual const& i2);
std::ostream& operator<<(std::ostream &out, Individual const& i1);

#endif
