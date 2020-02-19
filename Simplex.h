#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <iostream>
#include <bits/stdc++.h>
#include <math.h>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int p = 999999;

class Simplex{
public:
    double **coeffMatrix = 0;
    double **coeffMatrixPhase1 = 0;
    int m = 0,n = 0, a = 0, ukn = 0;

    bool phase1flag1 = false;
    bool basisContainsArtificial = false;
    bool hasSolution = false;
    vector<pair<int, double>> basisVector ;
    vector<pair<int, double>> basisVectorPhase1 ;//int = x_i th var from 1 but matrix has zero so take care
    vector<pair<int, double>> solution;
    double optimum = -p;

    bool intSolution = true;	
    
    MatrixXd convertToStandars(int nbVariables, int nbConstraints, VectorXd objectiveFunction, MatrixXd constraints){

        int n = nbVariables + nbConstraints + 1;
        int m = nbConstraints + 1;
        int a = 0;
        vector<int> vindex;
        for(int i = 0; i < nbConstraints; i++){
            if(constraints(i, nbVariables) < 0){
                a++;
                n++;
                vindex.push_back(i);
            }
        }
        vector<int>::iterator it;
        for(it = vindex.begin(); it < vindex.end(); it++){
            for(int i = 0; i < nbVariables+1; i++) constraints(*it, i) *= -1;
        }
        MatrixXd result(m, n);
        result = MatrixXd::Zero(m, n);
        for(int i = 0; i < m; i++){
            for(int j = 0; j < nbVariables; j++){
                if(i == 0){
                    result(i, j) = -objectiveFunction(j);
                } else {
                    result(i, j) = constraints(i-1, j);
                }                
            }
        }
        for(int i = 0; i < nbConstraints; i++){
            if(find(vindex.begin(), vindex.end(), i) < vindex.end()) result(1+i, nbVariables+i) = -1;
            else result(1+i, nbVariables+i) = 1;
        }
        int i = nbVariables + nbConstraints;
        for(it = vindex.begin(); it < vindex.end(); it++){
            result((*it)+1, i) = 1;
            i++;
        }
        for(int j = 0; j < m; j++) if(j > 0) result(j, n-1) = constraints(j-1, nbVariables);

        for(it = vindex.begin(); it < vindex.end(); it++){
            for(int i = 0; i < n ; i++) result(0, i) -= result((*it)+1, i);
        }

        for(int i = nbConstraints + nbVariables; i < n-1 ; i++) result(0, i) = 1;
        return result;
    }

    Simplex(){}

    Simplex(VectorXd objectiveFunction, MatrixXd constraints)
    {
        int nbVariables = objectiveFunction.size();
        int nbConstraints = constraints.rows();
        int i = 0,j = 0;
        n = nbVariables + nbConstraints + 1;
        m = nbConstraints + 1;
        for(int i = 0; i < nbConstraints; i++){
            if(constraints(i, nbVariables) < 0){
                n++;
            }
        }
        a = n - (nbVariables + nbConstraints + 1);
        ukn = nbVariables;
        MatrixXd matrix = convertToStandars(nbVariables, nbConstraints, objectiveFunction, constraints);

        // dynamically allocated coefficient matrix for simplex tableau
        //solutionVector = new double* [n];
        coeffMatrix = new double* [m];
        coeffMatrixPhase1 = new double* [m];

        for(i = 0; i < m; i++)
        {
            coeffMatrix[i] = new double [n];
            coeffMatrixPhase1[i] = new double [n];
        }


        for(int i = 0 ; i < m ; i++ )
        {
            for( int j = 0 ; j < n ; j++ )
            {
                coeffMatrix[i][j] = 0;
                coeffMatrixPhase1[i][j] = 0;
            }
        }


        for( int i = 0 ; i < m ; i++ )
        {
            for( int j = 0 ; j < n ; j++ )
            {

            coeffMatrix[i][j] = matrix(i, j);

            }

        }
        //copy  coeff matrix to phase1 matrix
        for( int i = 0 ; i < m ; i++ )
        {
            for( int j = 0 ; j < n ; j++ )
            {
                if((0 == i)&&(j < (n-a-1)))
                {
                    //all vars except artificial are 0
                }
                else{
                coeffMatrixPhase1[i][j] = coeffMatrix[i][j];
                }
            }
        }

        findIdentityMatrixPhase1();
        simplexPhase1();
        

        //postphase1 procesing and case handling :
        if(phase1flag1)
        {
            if(coeffMatrixPhase1[0][n-1] == 0)//Z* = 0
            {
                //check if basis has artificial var or not
                for(int i = 0; i < basisVectorPhase1.size(); i++)
                {
                    if((basisVectorPhase1[i].first - 1 >= (n-a-1) ))//xi starts from 1
                    {
                        basisContainsArtificial = true;
                    }
                    //
                }
                simplexPhase2();
                solution =  basisVectorPhase1;
                vector<pair<int, double>>::iterator it;
                for(it = solution.begin(); it < solution.end(); it++) optimum += -matrix(0, it->first - 1)*it->second;
                optimum += p; 
                hasSolution = true;
            }
        }
        if(hasSolution){
            vector<pair<int, double>>::iterator it;
            for(it = solution.begin(); it < solution.end(); it++){
                if(ceil(it->second) != floor(it->second)){ intSolution = false; break;}
            }
        }
        else intSolution = false;

    }

    //calculates min reduced cost, to decide the index of entering column to basis, returns index of min reduced cost column
    //user shall check if its positive or negative to decide optimal solution is reached or not.
    int searchMinz()
    {
        int min_index = 0;
        int min_val = coeffMatrix[0][0];
        int j;
        for(j = 0; j < n -1; j++)
        {
            if(coeffMatrix[0][j] < 0)
            {
                break;//at least one negative z element
            }
        }

        if(j == n - 1)
        {
            return -1;
        }

        //else there is atleast one negative ele

        for(int i = 1; i < n -1; i++)
        {
            if(min_val > coeffMatrix[0][i])
            {
                min_val = coeffMatrix[0][i];
                min_index = i;
            }
        }

        if(min_val < 0)
        {
            int i = 0;
            //check if atleast one column value is positive , if not solution is unbounded
            for(i = 1; i < m; i++)
            {
                if(coeffMatrix[i][min_index] > 0)
                    break;
            }

            if(i == m)
            {
            //all elements in column are negative
            return -1;
            }
        }
        return min_index;

    }

    //calculates leaving basis vector index on the basis of min-ratio
    //sanity checks to ensure min ratio is positive only.
    int calMinRatio_LeavingBasisVector(int enteringVectorIndex)
    {
        double min_ratio_val = 999999, new_ratio;//default large min ratio
        int min_ratio_index = -1;//default row index;

        for(int i = 1; i < m; i++)//iterate all constraint rows
        {
            if(coeffMatrix[i][enteringVectorIndex] != 0)
            {
                new_ratio = coeffMatrix[i][n-1]/coeffMatrix[i][enteringVectorIndex];
                if(new_ratio < 0)
                {
                    continue;
                }
                else if(((min_ratio_val) > new_ratio ))
                {
                    min_ratio_val = new_ratio;
                    min_ratio_index = i;
                }

            }
        }

        if(min_ratio_val >= 0)
        {
            return min_ratio_index;//return leaving  basis row;
        }
        else
        {
            return -1;// all min ratios are negative
        }
    }

    //update the solution vector at every iteration
    void updateSolutionvector(int leavingIndex, int enteringIndex)
    {
        basisVector[leavingIndex - 1].first = enteringIndex+1;//stores x1, x2 from 0th index of vector
        basisVector[leavingIndex - 1].second = coeffMatrix[leavingIndex][n-1];
        for(int i = 0; i < basisVector.size(); i++)
        {
            if(basisVector[i].first != enteringIndex+1)
            {
                basisVector[i].second = coeffMatrix[i+1][n-1];
            }
        }

    }

    int findIdentityMatrixPhase1()
    {
    //needs #unknows = m-1
    int countNonzero = 0;
    int Numslacks = n- a -ukn -1;
        for(int i = ukn; i < n-1 ; i++)//each col except last as it has nothing for identity chkin
        {
            for(int j = 1; j < m; j++) //evry row inc z
            {
                if(i < (ukn + Numslacks))//for slack/surplus cols
                {
                    if(coeffMatrixPhase1[j][i] == 1)
                    {
                        basisVectorPhase1.push_back(pair<int, double> (i + 1, coeffMatrixPhase1[j][n-1]));
                        break;
                    }
                }
                else{//for artificial vars
                        if(coeffMatrixPhase1[j][i] == 1)
                        {
                            basisVectorPhase1.push_back(pair<int, double> (i + 1, coeffMatrixPhase1[j][n-1]));

                            //make tableau a proper tableau
                            for(int k = 0; k < n; k++)//each column k of rowj effect row 0
                            {
                                coeffMatrixPhase1[0][k] = coeffMatrixPhase1[0][k] - (coeffMatrixPhase1[j][k]);
                            }


                        }

                }
            }

        }

        for(int i = 0; i < basisVectorPhase1.size(); i++)
        {
            //cout <<"x"<<basisVectorPhase1[i].first<<" = "<<basisVectorPhase1[i].second<<endl;

            basisVector.push_back(pair<int, double> (basisVectorPhase1[i].first, basisVectorPhase1[i].second));
        }

        for(int i = 0; i < basisVectorPhase1.size(); i++)
        {
            //cout <<"x"<<basisVectorPhase1[i].first<<" = "<<basisVectorPhase1[i].second<<endl;
            int colIndex = basisVector[i].first - 1;
            for(int j = 1; j <m; j++)
            {
                if(coeffMatrixPhase1[j][colIndex] == 1)
                {
                    basisVectorPhase1[j-1].first = basisVector[i].first;
                    basisVectorPhase1[j-1].second = basisVector[i].second;
                    break;
                }
            }
        }


    }

    void simplexPhase1()
    {
        while(1)
        {
            int entering_index =  searchMinzPhase1();
            if((entering_index == -1)||(coeffMatrixPhase1[0][entering_index] > 0))
            {
                break;//optimal or unbounded solution
            }
            int leaving_index = calMinRatio_LeavingBasisVectorPhase1(entering_index);

            if(leaving_index == -1)
            {
                break;
            }

            double pivot_element = coeffMatrixPhase1[leaving_index][entering_index];

            //make pivot element unity along with row operation by dividing whole row by pivot element
            if(pivot_element != 0){
                for(int j = 0; j < n; j++)
                {
                    coeffMatrixPhase1[leaving_index][j] /= pivot_element;
                }
            }

            pivot_element = coeffMatrixPhase1[leaving_index][entering_index];//resetting the pivot element to the value after the roe operation ie 1
            
            //make all elements of pivot column 0 expect pivot ele, Ri = Ri - Rp*Ri
            for( int i = 0 ; i < m ; i++ )
            {
                double multiply_factor = coeffMatrixPhase1[i][entering_index];
                if(leaving_index != i)
                {
                    for(int j = 0; j < n; j++)
                    {
                        coeffMatrixPhase1[i][j] = coeffMatrixPhase1[i][j] - ((multiply_factor) * (coeffMatrixPhase1[leaving_index][j]));
                        if(fabs(coeffMatrixPhase1[i][j]) < 0.000001)
                        {
                            coeffMatrixPhase1[i][j] = 0;

                        }

                    }
                }
            }

            updateSolutionvectorPhase1(leaving_index, entering_index);

        }


    }

    void updateSolutionvectorPhase1(int leavingIndex, int enteringIndex)
    {


        basisVectorPhase1[leavingIndex - 1].first = enteringIndex+1;//stores x1, x2 from 0th index of vector
        basisVectorPhase1[leavingIndex - 1].second = coeffMatrixPhase1[leavingIndex][n-1];
        for(int i = 0; i < basisVectorPhase1.size(); i++)
        {
            if(basisVectorPhase1[i].first != enteringIndex+1)
            {
                basisVectorPhase1[i].second = coeffMatrixPhase1[i+1][n-1];
            }
        }
    }

    void updateSolutionvectorPhase2(int leavingIndex, int enteringIndex)
    {
    ///basis vector of phase 1 is used but same coeffmatrix is used in phase2

        basisVectorPhase1[leavingIndex - 1].first = enteringIndex+1;//stores x1, x2 from 0th index of vector
        basisVectorPhase1[leavingIndex - 1].second = coeffMatrix[leavingIndex][n-1];
        for(int i = 0; i < basisVectorPhase1.size(); i++)
        {
            if(basisVectorPhase1[i].first != enteringIndex+1)
            {
                basisVectorPhase1[i].second = coeffMatrix[i+1][n-1];
            }
        }
    }

    int searchMinzPhase1()
    {
        int min_index = 0;
        int min_val = coeffMatrixPhase1[0][0];
        int j;
        for(j = 0; j < n-1; j++)
        {
            if(coeffMatrixPhase1[0][j] < 0)
            {
                break;//at least one negative z element
            }
        }

        if(j == n-1)
        {
            phase1flag1 = true;
            return -1;
        }

        //else there is atleast one negative ele


        for(int i = 1; i < n-1; i++)//n-1 so as to not to compare with value of z as last col is z
        {
            if(min_val > coeffMatrixPhase1[0][i])
            {
                min_val = coeffMatrixPhase1[0][i];
                min_index = i;
            }
        }

        if(min_val < 0)
        {
            int i = 0;
            //check if atleast one column value is positive , if not solution is unbounded
            for(i = 1; i < m; i++)
            {
                if(coeffMatrixPhase1[i][min_index] > 0)
                    break;
            }

            if(i == m)
            {
            return -1;
            }
        }
        return min_index;

    }

    int calMinRatio_LeavingBasisVectorPhase1(int enteringVectorIndex)
    {
        double min_ratio_val = 999999, new_ratio;//default large min ratio
        int min_ratio_index = -1;//default row index;

        for(int i = 1; i < m; i++)//iterate all constraint rows
        {
            if(coeffMatrixPhase1[i][enteringVectorIndex] != 0)
            {
                new_ratio = coeffMatrixPhase1[i][n-1]/coeffMatrixPhase1[i][enteringVectorIndex];
                if(new_ratio < 0)
                {
                    continue;
                }
                else if(((min_ratio_val) > new_ratio ))
                {
                    min_ratio_val = new_ratio;
                    min_ratio_index = i;
                }
            }
        }

        if(min_ratio_val >= 0)
        {
            return min_ratio_index;//return leaving  basis row;
        }
        else
        {
            return -1;// all min ratios are negative
        }

    }

    void simplexPhase2()
    {
    if(!basisContainsArtificial)
    {
        //copy orignal z to coeffmatrix
        for(int j = 0 ;j < n; j++)
        {
            coeffMatrixPhase1[0][j] = coeffMatrix[0][j];
        }
        //remove artificial vars columns
        //copy n-1th col to first artificial var col and reduce n
        //first artificial var col n-a-1
        for(int i = 0; i < m; i++)
        {
            coeffMatrixPhase1[i][n-a-1] = coeffMatrixPhase1[i][n-1];
        }
        //now set new n
        n = n-a;

            //to work with coeffmatrix
            for( int i = 0 ; i < m ; i++ )
            {
                for( int j = 0 ; j < n ; j++ )
                {

                    coeffMatrix[i][j] = coeffMatrixPhase1[i][j] ;

                }
            }
    }
    else//z* = 0 but contains a artificial var with val 0
    {
        for(int j = 0 ;j < n; j++)
        {
            coeffMatrixPhase1[0][j] = coeffMatrix[0][j];
        }

        //remove artificial vars columns which is not in basis
        //set whole column to 0
        for(int j = (n-a-1) ;j < n-1; j++)
        {
            if(!artColinBasis(j+1))//var indexed from 1 is sent
            {
            //set all rows of this col = 0
                for(int i = 0 ;i < m; i++)
                {
                    coeffMatrixPhase1[i][j] = 0;
                }
            }

        }

        //to work with coeffmatrix
        for( int i = 0 ; i < m ; i++ )
        {
            for( int j = 0 ; j < n ; j++ )
            {

                coeffMatrix[i][j] = coeffMatrixPhase1[i][j] ;

            }
        }


    }

    //make tableu legitimiate
        for(int i = 0; i < basisVectorPhase1.size(); i++)
        {
            //cout <<"x"<<basisVectorPhase1[i].first<<" = "<<basisVectorPhase1[i].second<<endl;
            int colofBasis = basisVectorPhase1[i].first - 1;
            if(coeffMatrix[i+1][colofBasis] == 1)
            {
                //make tableau a proper tableau
                int multplier = coeffMatrix[0][colofBasis];
                for(int k = 0; k < n; k++)//each column k of rowj effect row 0
                {
                    coeffMatrix[0][k] = coeffMatrix[0][k] - ((multplier)*coeffMatrix[i+1][k]);
                }
            }
        }

    //#if 0
        //now solve  within loop iterations
        while(1)
        {
            int entering_index =  searchMinz();
            if((entering_index == -1)||(coeffMatrix[0][entering_index] > 0))
            {
                break;//optimal or unbounded solution
            }
            int leaving_index = calMinRatio_LeavingBasisVector(entering_index);

            if(leaving_index == -1)
            {
                break;
            }

            double pivot_element = coeffMatrix[leaving_index][entering_index];

            //make pivot element unity along with row operation by dividing whole row by pivot element
            if(pivot_element != 0){
                for(int j = 0; j < n; j++)
                {
                    coeffMatrix[leaving_index][j] /= pivot_element;
                }
            }

            pivot_element = coeffMatrix[leaving_index][entering_index];//resetting the pivot element to the value after the roe operation ie 1

            //make all elements of pivot column 0 expect pivot ele, Ri = Ri - Rp*Ri
            for( int i = 0 ; i < m ; i++ )
            {
                double multiply_factor = coeffMatrix[i][entering_index];
                if(leaving_index != i)
                {
                    for(int j = 0; j < n; j++)
                    {
                        coeffMatrix[i][j] = coeffMatrix[i][j] - ((multiply_factor) * (coeffMatrix[leaving_index][j]));
                        if(fabs(coeffMatrix[i][j]) < 0.000001)
                        {
                            coeffMatrix[i][j] = 0;

                        }

                    }

                }
            }

            updateSolutionvectorPhase2(leaving_index, entering_index);
            
        }

    //#endif
    }

    bool artColinBasis(int colIndex)//col index starts from 1
    {
        for(int i = 0; i < basisVectorPhase1.size(); i++)
        {
            if(basisVectorPhase1[i].first == colIndex) //belongs to artificial vars col
            {
                return true;
            }
        }
        return false;
    }

};
#endif