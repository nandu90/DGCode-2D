/* 
 * File:   controlinput.h
 * Author: User
 *
 * Created on February 5, 2016, 12:57 PM
 */

#ifndef CONTROLINPUT_H
#define	CONTROLINPUT_H
#include <string>
#include <sstream>
#include <iostream>

//Reads the control.txt file and the inputs specified therein
void control(int &stages, int &totalsteps, int &fluxscheme, int &grid, int &method, int &limiter)
{
    ifstream controlfile;
    controlfile.open("control.txt");
    string line;
    while(getline(controlfile,line))
    {
        string junk;
        int junk2=-1;
        istringstream ss(line);
        ss>>junk;   
        if(junk.compare("rho_inf") == 0)
        {
            ss >> rho_inf;
        }
        else if(junk.compare("v_inf") == 0)
        {
            ss >> v_inf;
        }
        else if(junk.compare("M_inf") == 0)
        {
            ss >> M_inf;
        }
        else if(junk.compare("gamma") == 0)
        {
            ss >> g;
        }
        else if(junk.compare("angle_of_attack") == 0)
        {
            ss >> alpha;
        }
        else if(junk.compare("CFL") == 0)
        {
            ss >> cfl;
        }
        else if(junk.compare("Stages") == 0)
        {
            ss >> stages;
        }
        else if(junk.compare("Totalsteps") == 0)
        {
            ss >> totalsteps;
        }
        else if(junk.compare("Van-Leer") == 0)
        {
            fluxscheme=1;
        }
        else if(junk.compare("Reconstruct") == 0)
        {
            fluxscheme=2;
        }
        else if(junk.compare("Barth-Jasperson") == 0)
        {
            limiter=1;
        }
        else if(junk.compare("Venkatakrishnan") == 0)
        {
            limiter=2;
        }
        else if(junk.compare("NACA0012") == 0)
        {
            grid=1;
        }
        else if(junk.compare("Cylinder-medium") == 0)
        {
            grid=2;
        }
        else if(junk.compare("Cylinder-fine") == 0)
        {
            grid=3;
        }
        else if(junk.compare("Cylinder-veryfine") == 0)
        {
            grid=4;
        }
        else if(junk.compare("Cylinder-coarse") == 0)
        {
            grid=5;
        }
        else if(junk.compare("Channel") == 0)
        {
            grid=6;
        }
        else if(junk.compare("test") == 0)
        {
            grid=0;
        }
        else if(junk.compare("Global-Step") == 0)
        {
            method=1;
        }
        else if(junk.compare("Local-Step") == 0)
        {
            method=2;
        }
    }
    controlfile.close();
}
#endif	/* CONTROLINPUT_H */

