/* 
 * File:   rhs_vector.h
 * Author: User
 *
 * Created on March 31, 2016, 4:40 PM
 */

#ifndef RHS_VECTOR_H
#define	RHS_VECTOR_H

struct elemvalues
{
    int no;
    double rho;
    double u;
    double v;
    double p;
    double c;
};

#include "higherOrderReconstruct.h"
/*********************Upwind Van-Leer Flux routine*****************************/
void fluxes(vector<double> &flux, int fluxscheme, elemvalues elemi, elemvalues elemj, vector< vector<double> >gradWi, vector< vector<double> >gradWj, int limiter, vector< vector< vector<double> > > Usurr, vector< vector< vector<double> > > local_vectofaces, vector< vector<double> > vectoface, double nx, double ny, vector<double> area)
{
    vector<double> Fl(neqns,0);
    vector<double> Fr(neqns,0);
    
    double uL,vL,rhoL,pL,cL;
    double uR,vR,rhoR,pR,cR;
    
    if(fluxscheme == 1) //1st Order - Face value corresponds to element center values
    {
        rhoL=elemi.rho;
        uL=elemi.u;
        vL=elemi.v;
        pL=elemi.p;
        
        rhoR=elemj.rho;
        uR=elemj.u;
        vR=elemj.v;
        pR=elemj.p;
    }  
    else if (fluxscheme == 2)
    {
        vector<double> WL(neqns,0);
        vector<double> WR(neqns,0);
        higherorder(elemi, elemj, WL, WR, gradWi, gradWj,limiter,Usurr,local_vectofaces,vectoface,area);
        rhoL=WL[0];
        uL=WL[1];
        vL=WL[2];
        pL=WL[3];
        
        rhoR=WR[0];
        uR=WR[1];
        vR=WR[2];
        pR=WR[3];
        
    }
    
    cL=sqrt(fabs(g*pL/rhoL));
    cR=sqrt(fabs(g*pR/rhoR));
    /***Calculate the normal component of velocities w.r.t. face***/
        /**Dot Product**/
    double vnL=uL*nx+vL*ny;
    double vnR=uR*nx+vR*ny;
    
    /***Calculate the respective Mach Numbers***/
    double MnL=vnL/cL;
    double MnR=vnR/cR;
        
    if(MnL >= 1)
    {
        Fl[0]=rhoL*vnL;
        Fl[1]=rhoL*vnL*uL + pL*nx;
        Fl[2]=rhoL*vnL*vL + pL*ny;
        double rhoE=pL/(g-1) + 0.5*rhoL*(pow(uL,2) + pow(vL,2));
        Fl[3]=vnL*(rhoE + pL);
    }
        
    else if (MnL <= -1)
    {
        Fl[0]=0;
        Fl[1]=0;
        Fl[2]=0;
        Fl[3]=0;
    }
    else 
    {
        double mflux=rhoL*cL*pow(MnL + 1,2)/4;
        Fl[0]=mflux;
        Fl[1]=mflux*(uL + (nx*(-vnL + 2*cL))/g);
        Fl[2]=mflux*(vL + (ny*(-vnL + 2*cL))/g);
        Fl[3]=mflux*((((pow(uL,2)+pow(vL,2)) - pow(vnL,2))/2.0) + ((pow((g-1)*vnL + 2*cL,2))/(2*(pow(g,2)-1))));
    }
        
        
    if(MnR >= 1)
    {
        Fr[0]=0;
        Fr[1]=0;
        Fr[2]=0;
        Fr[3]=0;
    }
    else if(MnR <= -1)
    {
        Fr[0]=rhoR*vnR;
        Fr[1]=rhoR*vnR*uR + pR*nx;
        Fr[2]=rhoR*vnR*vR + pR*ny;
        double rhoE=pR/(g-1) + 0.5*rhoR*(pow(uR,2) + pow(vR,2));
        Fr[3]=vnR*(rhoE + pR);
    }
    else
    {
        double mflux=-rhoR*cR*pow(MnR - 1,2)/4;
        Fr[0]=mflux;
        Fr[1]=mflux*(uR + (nx*(-vnR - 2*cR))/g);
        Fr[2]=mflux*(vR + (ny*(-vnR - 2*cR))/g);
        Fr[3]=mflux *((((pow(uR,2)+pow(vR,2)) - pow(vnR,2))/2.0) + ((pow((g-1)*vnR - 2*cR,2))/(2*(pow(g,2)-1))));
    }
        
    for(int i=0; i<neqns; i++)
    {
        flux[i]=Fl[i]+Fr[i];
    }
}

void rhsvector(vector< vector<double> > *R, vector< vector<double> > *U, vector<fdata> *intfac, vector< vector<int> > *esuel, int fluxscheme, vector< vector< vector<double> > > *gradW, int limiter, vector< vector< vector<double> > > *vec_to_faces, vector< vector<double> > *geoel)
{
    
    
    /*Reinitialize R vector*/
    for(int i=0; i<nelem; i++)
    {
        for (int j=0; j<neqns; j++)
        {
            (*R)[i][j]=0;
        }
    }
    #pragma omp parallel for schedule (dynamic)
    for(int i=0; i< (*intfac).size(); i++)
    {
        
        vector<double> faceflux(neqns,0);
        /***Note this is only Finite Volume, therefore variables will have the center value at the face***/        
        elemvalues elemi;
        cons_to_prim((*U)[(*intfac)[i].ei -1],elemi.rho,elemi.u,elemi.v,elemi.p);
        elemi.c=sqrt(fabs(g*elemi.p/elemi.rho));
        elemi.no=(*intfac)[i].ei;
        
        elemvalues elemj;
        cons_to_prim((*U)[(*intfac)[i].ej -1],elemj.rho,elemj.u,elemj.v,elemj.p);
        elemj.c=sqrt(fabs(g*elemj.p/elemj.rho));
        elemj.no=(*intfac)[i].ej;
        
        double nx=(*intfac)[i].nx;
        double ny=(*intfac)[i].ny;
        double length=(*intfac)[i].len;
        
        /*********************For Reconstruction*************************/
        /**Find the Gradients at cell centers***/
        vector< vector<double> > gradWi(neqns, vector<double> (2,0));
        vector< vector<double> > gradWj(neqns, vector<double> (2,0));
    
        for(int i=0; i< neqns; i++)
        {
           for(int j=0; j<2; j++)
           {
               gradWi[i][j]=(*gradW)[elemi.no-1][i][j];
               if(elemj.no <= nelem)
               {
                   gradWj[i][j]=(*gradW)[elemj.no-1][i][j];
               }
           }
        }
        /****vector from cell center of i and j to face center*****/ 
        vector< vector<double> > vectoface(2,vector<double>(2,0));
        for(int j=0;j<3;j++)
        {
            if((*esuel)[elemi.no-1][j] == elemj.no)
            {
                vectoface[0][0]=(*vec_to_faces)[elemi.no-1][j][0];
                vectoface[0][1]=(*vec_to_faces)[elemi.no-1][j][1];
                break;
            }
        }
        if(elemj.no <= nelem)
        {
            for(int j=0; j<3; j++)
            {
                if((*esuel)[elemj.no-1][j] == elemi.no)
                {
                    vectoface[1][0]=(*vec_to_faces)[elemj.no-1][j][0];
                    vectoface[1][1]=(*vec_to_faces)[elemj.no-1][j][1];
                    break;
                }
            }
        }
        /******************************************************************/
        
        /********************For Limiter application***************************/
        vector< vector< vector<double> > > Usurr(2, vector< vector<double> > (3,vector<double> (neqns,0)));
        for(int j=0; j<3; j++)
        {
            for(int k=0; k<neqns; k++)
            {
                Usurr[0][j][k] = (*U)[(*esuel)[elemi.no-1][j]-1][k];
            }
        }
        
        for(int j=0; j<3; j++)
        {
            for(int k=0; k<neqns; k++)
            {
                if(elemj.no <= nelem)
                {
                    Usurr[1][j][k] = (*U)[(*esuel)[elemj.no-1][j]-1][k];
                }
            }
        }
        
        
        vector< vector< vector<double> > > local_vectofaces(2, vector< vector<double> > (3,vector<double> (2,0)));
        for(int j=0; j<3; j++)
        {
            local_vectofaces[0][j][0]=(*vec_to_faces)[elemi.no-1][j][0];
            local_vectofaces[0][j][1]=(*vec_to_faces)[elemi.no-1][j][1];
        }
        
        for(int j=0; j<3; j++)
        {
            if(elemj.no <= nelem)
            {
                local_vectofaces[1][j][0]=(*vec_to_faces)[elemj.no-1][j][0];
                local_vectofaces[1][j][1]=(*vec_to_faces)[elemj.no-1][j][1];
            }
        }
        
        vector<double> area(2,0);
        area[0]= (*geoel)[elemi.no-1][0];
        if(elemj.no <= nelem)
        {
            area[1]=(*geoel)[elemj.no-1][0];
        }
        /****************************************************************/
        
        
        
        
        fluxes(faceflux,fluxscheme, elemi, elemj,gradWi,gradWj,limiter, Usurr, local_vectofaces, vectoface,nx,ny,area);
        
        for(int j=0; j<neqns; j++)
        {
            faceflux[j]=faceflux[j]*length;
        }
        
        /***Add this flux to left element and subtract from right element***/

        for(int j=0; j<neqns; j++)
        {
            # pragma omp critical
            {
            (*R)[(*intfac)[i].ei-1][j]=(*R)[(*intfac)[i].ei-1][j] + faceflux[j];
            if((*intfac)[i].ej <= nelem)
            {
                (*R)[(*intfac)[i].ej-1][j]=(*R)[(*intfac)[i].ej-1][j] - faceflux[j];
            }
            }
        }
        
    }
    
    for(int i=0; i<nelem; i++)
    {
        for (int j=0; j<neqns; j++)
        {
            (*R)[i][j]=-(*R)[i][j];
        }
    }
}

#endif	/* RHS_VECTOR_H */

