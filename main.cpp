/* 
 * File:   main.cpp
 * Author: User
 *
 * Created on February 2, 2016, 1:21 PM
 */
#define PI 3.1415926535897
#include <iostream>
# include <stdlib.h>
# include <math.h>
# include <fstream>
#include <cmath>
#include <time.h>
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>
#include <fenv.h>
#include <stdio.h>
#include <limits.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <omp.h>
#include <time.h>
#include <cctype>
#include <cstring>

///Global Variable Declaration
    double cfl;
    double rho_inf; 
    double v_inf;
    double M_inf;
    double g;              //gamma
    double alpha;          //Angle of attack
    double s_inf;
    
    int ndimn;
    int nelem;
    int npoin;
    int nface;
    int neqns=4;
    
using namespace std;

string doubletostr (double n)
{
    stringstream ss;
    ss<<n;
    return ss.str();
}


///Function to get the path of current directory
string getexepath() {
   char cwd[1024];
   stringstream ss2;
   string dir;
   if (getcwd(cwd, sizeof(cwd)) != NULL)
   {
       ss2 << cwd;
       ss2 >> dir;
   }
   else
       perror("getcwd() error");
   return dir;
}

///Function to select the case being simulated
string dirpath(int grid)
{
    string in;
    if(grid == 0)
    {
        in="/test";
    }
    if(grid == 1)
    {
        in="/NACA0012";
    }
    else if (grid == 2)
    {
        in="/Cylinder-medium";
    }
    else if (grid == 3)
    {
        in="/Cylinder-fine";
    }
    else if (grid == 4)
    {
        in="/Cylinder-veryfine";
    }
    else if (grid == 5)
    {
        in="/Cylinder-coarse";
    }
    else if (grid == 6)
    {
        in="/Channel";
    }
    
    return in;
}

///Header files required for the program
#include "gridinput.h"
#include "controlinput.h"
#include "preprocessing.h"
#include "init_fields.h"
#include "output_gen.h"
#include "rhs_vector.h"



int main()
{
    ///Catches mathematical exceptions
    feenableexcept(FE_INVALID | FE_OVERFLOW |FE_DIVBYZERO);
    
    ///Environment variable to select number of processors
    omp_set_num_threads(max(atoi(getenv("OMP_NUM_THREADS")),1));     
    
    /**********Declaration of variables****************************************/
    vector< vector<double> > coord(1,vector<double> (1));
    vector< vector<int> > inpoel(1,vector<int> (1));
    vector< vector<int> > bface(1,vector<int> (1));
    int stages;
    int totalsteps;
    int fluxscheme;
    int grid;
    int method;
    int limiter=0;
    int gradmethod=1;
    
    /*********************Read Control File************************************/    
    control(stages, totalsteps,fluxscheme,grid,method,limiter);
    /**************************************************************************/
    
    /**Print Primitive variables and Residuals in a file for each iteration*/
    string in=dirpath(grid);
    ofstream out;
    out.open(getexepath() +in +"/conservation.txt",ios::trunc);
    out.close();
    
    out.open(getexepath() +in +"/residuals.txt",ios::trunc);
    out.close();
    /**************************************************************************/
     
    /*********************Read Grid file***************************************/
    gridcall(coord,inpoel,bface,grid);    
    /**************************************************************************/
        
    /****************************Preprocessing*********************************/
    vector<vector<double> > geoel(nelem,vector<double> (3,0));
    vector< vector <int> > esuel (nelem, vector<int> (3,0));
    vector<fdata> intfac;
    int nbfac,nafac;
    preprocess(coord,inpoel,bface, geoel, intfac,esuel,nbfac,nafac);
    
    /*Find the Centroids of Ghost Cells for higher order methods*/
    vector<ghostdata> geogh;
    ghostcentroid(geogh,intfac,nbfac,coord,geoel);
    
    /*Find the Least Square matrix for higher order reconstruction*/
    vector< vector<double> > lsquare(nelem, vector<double>(3,0));
    lsquare_matrix(esuel, geoel, lsquare, geogh);
    
    /*Store the vectors connecting Element centroid to Gauss point*/
    vector< vector< vector<double> > > vec_to_faces(nelem, vector< vector<double> > (3, vector<double> (2,0)));
    elem_vectofaces(vec_to_faces, intfac, coord, esuel, geoel);
    /**************************************************************************/ 
    
    
    /*************************Declare vectors for variables********************/
    vector <vector<double> > W(nelem,vector<double> (neqns,0));  //Primitive Variables
    vector <vector<double> > U(nelem,vector<double> (neqns,0));  //Conservative Variables
    /**************************************************************************/
    
    /*****************************Initialize vectors***************************/
    initialize(W, U);
    /**************************************************************************/
    
    /*************************Assign Boundary Conditions***********************/
    boundary(intfac,bface,U);
    /**************************************************************************/
    
    /**************************Initialize RHS vector***************************/
    vector< vector<double> > R(nelem, vector<double>(neqns,0));    
    /**************************************************************************/
    
    /****************Initialize Gradient vector- for higher order**************/
    vector< vector< vector<double> > > gradW(nelem, vector< vector<double> > (neqns, vector<double> (2,0)));
    /**************************************************************************/
    
    bool exitflag=false;
    int t=0;
    vector<double> ires(neqns,0);
    
    int totaliter=0;
    
    time_t t1,t2;
    t1=time(0);
    
    /************************Iteration Loop************************************/
    while (t < totalsteps)
    {
        double deltat;
        vector<double> dt(nelem,0);
        time_step(&U, &geoel, &intfac, &esuel, deltat, &dt);
        if(method == 1)
        {
            #pragma omp parallel for schedule (dynamic)
            for(int i=0; i<nelem; i++)
            {
                dt[i]=deltat;
            }
        }
        
        vector< vector<double> > Uini(nelem, vector<double>(neqns,0));
        #pragma omp parallel for schedule (dynamic)
        for(int i=0; i<nelem; i++)
        {
            for(int j=0; j<neqns; j++)
            {
                Uini[i][j]=U[i][j];
            }
        }
        /*********************1-Stage Runge-Kutta******************************/
        if (stages == 1)
        {
            if(fluxscheme != 1)
            {
                cellcentergrad(&lsquare, &esuel, &geogh, &geoel, &U, &gradW);
            }
        rhsvector(&R, &U, &intfac, &esuel, fluxscheme,&gradW,limiter,&vec_to_faces,&geoel);
        euler (&U, &geoel, &R, &dt);
        bound_update(&bface, &intfac, &U);
        }
        
        
        /*********************2-Stage Runge-Kutta******************************/
        if (stages ==2)
        {
            vector< vector<double> > Utemp(nelem, vector<double> (neqns,0));
            for (int i=0;i<nelem;i++)
            {
                for (int k=0;k<neqns;k++)
                {
                    Utemp[i][k]=U[i][k];
                }
            }
            
            if(fluxscheme != 1)
            {
                cellcentergrad(&lsquare, &esuel, &geogh, &geoel, &U, &gradW);
            }
            rhsvector(&R, &U, &intfac, &esuel, fluxscheme,&gradW,limiter,&vec_to_faces,&geoel);
            euler (&U, &geoel, &R, &dt);
            bound_update(&bface, &intfac, &U);
            
            
            if(fluxscheme != 1)
            {
                cellcentergrad(&lsquare, &esuel, &geogh, &geoel, &U, &gradW);
            }
            rhsvector(&R, &U, &intfac, &esuel, fluxscheme,&gradW,limiter,&vec_to_faces,&geoel);
            euler (&U, &geoel, &R, &dt);
            
            #pragma omp parallel for schedule (dynamic)
            for (int i=0;i<nelem;i++)
            {
                for (int k=0;k<neqns;k++)
                {
                    U[i][k]= 0.5*Utemp[i][k]+0.5*U[i][k];
                }
            }
            bound_update(&bface, &intfac, &U);
            
        }
        /*********************3-Stage Runge-Kutta******************************/
        else if (stages == 3)
        {
            vector< vector<double> > Utemp(nelem, vector<double> (neqns,0));
            
            for (int i=0;i<nelem;i++)
            {
                for (int k=0;k<neqns;k++)
                {
                    Utemp[i][k]=U[i][k];
                }
            }
            
            if(fluxscheme != 1)
            {
                cellcentergrad(&lsquare, &esuel, &geogh, &geoel, &U, &gradW);
            }
            rhsvector(&R, &U, &intfac, &esuel, fluxscheme,&gradW,limiter,&vec_to_faces,&geoel);
            euler (&U, &geoel, &R, &dt);
            bound_update(&bface, &intfac, &U);
            
            
            if(fluxscheme != 1)
            {
                cellcentergrad(&lsquare, &esuel, &geogh, &geoel, &U, &gradW);
            }
            rhsvector(&R, &U, &intfac, &esuel, fluxscheme,&gradW,limiter,&vec_to_faces,&geoel);
            euler (&U, &geoel, &R, &dt);
            bound_update(&bface, &intfac, &U);
            
            #pragma omp parallel for schedule (dynamic)
            for (int i=0;i<nelem;i++)
            {
                for (int k=0;k<neqns;k++)
                {
                    U[i][k]= (3.0/4.0)*Utemp[i][k]+(1.0/4.0)*U[i][k];
                }
            }
            bound_update(&bface, &intfac, &U);
            
            
            if(fluxscheme != 1)
            {
                cellcentergrad(&lsquare, &esuel, &geogh, &geoel, &U, &gradW);
            }
            rhsvector(&R, &U, &intfac, &esuel, fluxscheme,&gradW,limiter,&vec_to_faces,&geoel);
            euler (&U, &geoel, &R, &dt);
            
            #pragma omp parallel for schedule (dynamic)
            for (int i=0;i<nelem;i++)
            {
                for (int k=0;k<neqns;k++)
                {
                    U[i][k]= (1.0/3.0)*Utemp[i][k]+(2.0/3.0)*U[i][k];
                }
            }
            bound_update(&bface, &intfac, &U);
        }
        t++;
        monitor_cons(&U, &geoel, t, exitflag,grid);
        
        monitor_res(&Uini,&U,&geoel,t, exitflag,ires,grid);
        if(exitflag == true)
        {
            totaliter=t;
            cout<<"Solution Converged!! Hurrrray!!!"<<endl;
            break;
        }
        
    }
    
        
    t2=time(0);
    
    double seconds= difftime(t2,t1);
    
    cout<<"Total run time: "<<seconds<<" secs"<<endl;

    /****************Generate Tecplot Output file******************************/
    output(geoel,coord,inpoel,U,grid,totaliter,seconds);
    
} 