#include <cinttypes>
#include <iostream>
#include <bits/stdc++.h>
#include <Eigen/Dense>
#include "Simplex.h"

using namespace std;
using namespace Eigen;

static int courant = -999;

#define LESS 0
#define MORE 1

MatrixXd addConstraint(int mode,int index, double val, MatrixXd constraints){
	MatrixXd result(constraints.rows()+1, constraints.cols());
	MatrixXd v(1, constraints.cols());
	v = MatrixXd::Zero(1, constraints.cols());
	if(mode == LESS){
		v(0, index-1) = 1;
		v(0, constraints.cols()-1) = val;
	}
	else{
		v(0, index-1) = -1;
		v(0, constraints.cols()-1) = -val;
	}
	result << constraints, v;
	return result;
}

bool bound(Simplex solver){

	bool truth = true;
	if(!solver.hasSolution)	truth = false;
	if(solver.optimum < courant){
	
		truth = false;
	
	} else {
	
		courant = solver.optimum; 
	
	}
	return truth;

}

vector<MatrixXd> branch(MatrixXd constraints, pair<int, double>index){
	vector<MatrixXd> result;
	result.push_back(addConstraint(LESS, index.first, floor(index.second), constraints));
	result.push_back(addConstraint(MORE, index.first, ceil(index.second), constraints));
	return result;
}

Simplex BB(VectorXd objectiveFunction, MatrixXd constraints){
	Simplex solver(objectiveFunction, constraints);
	if(solver.hasSolution){
		if(solver.intSolution) return solver;
		else{
			vector<pair<int, double>>::iterator it;
			for(it = solver.solution.begin(); it < solver.solution.end(); it++){
				if(ceil(it->second) != floor(it->second))	break;
			}
			vector<MatrixXd> branches = branch(constraints, *it);
			Simplex solver1 = BB(objectiveFunction, branches.at(0));
			Simplex solver2 = BB(objectiveFunction, branches.at(1));
			if(bound(solver1) && bound(solver2)){
				if(solver1.optimum > solver2.optimum)	return solver1;
				else return solver2;
			} else {
				if(bound(solver1))	return solver1;
				if(bound(solver2))	return solver2;
				solver.optimum = -999;
				solver.intSolution = false;
				return solver;
			}
		}
	} else {
		solver.optimum = -999;
		return solver;
	}

}

int main()
{
    int nbConstraints = 2;
    int nbVariables = 2;  
    MatrixXd constraints(nbConstraints, nbVariables + 1);
    VectorXd objectiveFunction(nbVariables);
    
    constraints << 10, 2, 17,
	       	  10, 3, 47;
    objectiveFunction << 5, 2;

	Simplex solver;
	solver = BB(objectiveFunction, constraints);
	if(solver.hasSolution && solver.intSolution){
		cout << "Solution :" << endl;
		for(int j = 1; j <= nbVariables; j++){
			vector<pair<int, double>>::iterator i;
			for(i = solver.solution.begin(); i < solver.solution.end(); i++)	if(i->first == j) break;

			if(i->first == j)	cout << "X" << i->first << " = " << i->second << endl;
			else cout << "X" << j << " = " << 0 << endl;
		}
		cout << "Optimum :" << endl;
		cout << "Z* = " << solver.optimum <<endl;

	}

	return 0;
}

