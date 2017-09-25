/* 
 * File:   rhs_vector.h
 * Author: User
 *
 * Created on March 31, 2016, 4:40 PM
 */

#ifndef HIGHERORDERRECONSTRUCT_H
#define	HIGHERORDERRECONSTRUCT_H

#include <algorithm>

/*********************Calculates Cell center gradients*************************/
void cellcentergrad(vector< vector<double> > *lsquare, vector< vector<int> > *esuel, vector<ghostdata> *geogh, vector< vector<double> > *geoel, vector< vector<double> > *U, vector< vector< vector<double> > > *gradW)
{
    /*Reintialize*/
    #pragma omp parallel for schedule (dynamic)
    for(int i=0; i<nelem; i++)
    {
        for(int j=0; j<neqns; j++)
        {
            (*gradW)[i][j][0]=0;
            (*gradW)[i][j][1]=0;
        }
    }
    /************/
    #pragma omp parallel for schedule (dynamic)
    for(int i=0; i<nelem; i++)
    {
       double a=(*lsquare)[i][0];
       double b=(*lsquare)[i][1];
       double c=(*lsquare)[i][2];
       vector<double> D(neqns,0);
       vector<double> E(neqns,0);
    
       double Disc=a*c-b*b;
    
       double xi=(*geoel)[i][1];
       double yi=(*geoel)[i][2];
       vector<double> tempwi(neqns,0);
       cons_to_prim((*U)[i],tempwi[0],tempwi[1],tempwi[2],tempwi[3]);
    
       for(int j=0; j<3; j++)
       {
           vector<double> tempwj(neqns,0);
           double xj,yj;
           int selem;
           if((*esuel)[i][j] > nelem)//ghost cell
           {
               for(int k=0; k<(*geogh).size(); k++)
               {
                  if((*esuel)[i][j] == (*geogh)[k].elem)
                  {
                     xj=(*geogh)[k].cenx;
                     yj=(*geogh)[k].ceny;
                     selem=(*geogh)[k].elem;
                     break;
                  }
                }
            }
            else
            {
                xj=(*geoel)[(*esuel)[i][j]-1][1];
                yj=(*geoel)[(*esuel)[i][j]-1][2];
                selem=(*esuel)[i][j];
            }
            cons_to_prim((*U)[selem-1],tempwj[0],tempwj[1],tempwj[2],tempwj[3]);
        
            double w=1/sqrt(pow(xj-xi,2) + pow(yj-yi,2));
        
            for(int k=0; k< neqns; k++)
            {
                D[k]=D[k]+pow(w,2)*(xj-xi)*(tempwj[k]-tempwi[k]);
                E[k]=E[k]+pow(w,2)*(yj-yi)*(tempwj[k]-tempwi[k]);
            }       
       }
    
    /***Solving Matrix System by Cramer's Rule***/
       for(int k=0; k<neqns; k++)
       {
            double DX=c*D[k] - b*E[k];
            double DY=a*E[k] - b*D[k];
            (*gradW)[i][k][0]=DX/Disc;
            (*gradW)[i][k][1]=DY/Disc;
       }
    }
}

/********************Implementation of Barth-Jasperson Limiter*****************/
void Barth_Jasperson(vector< vector<double> > Usurr, vector<double> W, vector< vector<double> > local_vectofaces, vector< vector<double> >gradW, vector<double> &phi)
{
    /***Convert the Conservative variables to Primitve for surrounding neighbours*/
    vector< vector<double> > Wsurr(3, vector<double> (neqns,0));
    for(int i=0; i<3; i++)
    {
        cons_to_prim(Usurr[i], Wsurr[i][0], Wsurr[i][1], Wsurr[i][2], Wsurr[i][3]);
    }
    for (int i=0; i<neqns; i++)
    {
        vector <double> deltaW(3,0);
        for (int j=0; j<3; j++)
        {
            deltaW[j]=fabs(Wsurr[j][i] - W[i]);
            //cout<<deltaW[j]<<" ";
        }
        //cout<<endl;
        double dWmax= *max_element(deltaW.begin(), deltaW.end());
        double dWmin= *min_element(deltaW.begin(), deltaW.end());
        
        /***Compute the Unconstrained Value at each Gauss Point***/
        vector<double> phigauss(3,0);
        vector<double> Wgauss(3,0);
        for (int j=0; j<3; j++)
        {
            Wgauss[j] = gradW[i][0]*local_vectofaces[j][0] + gradW[i][1]*local_vectofaces[j][1];
            
            double Wgaussplus=0;
            if(Wgauss[j] > 0)
            {
                Wgaussplus=dWmax;
            }
            else if (Wgauss[j] < 0)
            {
                Wgaussplus=dWmin;
            }
            
            if(Wgauss[j] == 0)
            {
                phigauss[j]=1;
            }
            else
            {
                double one=1;
                phigauss[j]= min(one, fabs(Wgaussplus/Wgauss[j]));
            }
        }
        
        phi[i] = *min_element(phigauss.begin(), phigauss.end());
    }
}


/*****************Function used by Venkatakrishnan limiter routine*************/
double limiter(double a , double b, double area)
{
    double k=0.3;
    double eps=pow(k*sqrt(area),1.5);
    
    double phi;
    if(b == 0)
    {
        phi=1;
    }
    else
    {
        phi=(1/b)*((a*a + eps*eps)*b + 2*b*b*a)/(a*a + 2*b*b + a*b + eps*eps);
    }
    return phi;
}

/********************Implementation of Venkatakrishnan Limiter*****************/
void venkatakrishnan(vector< vector<double> > Usurr, vector<double> W, vector< vector<double> > local_vectofaces, vector< vector<double> >gradW, vector<double> &phi,double area)
{
    /***Convert the Conservative variables to Primitve for surrounding neighbours*/
    vector< vector<double> > Wsurr(3, vector<double> (neqns,0));
    for(int i=0; i<3; i++)
    {
        cons_to_prim(Usurr[i], Wsurr[i][0], Wsurr[i][1], Wsurr[i][2], Wsurr[i][3]);
    }
    for (int i=0; i<neqns; i++)
    {
        vector <double> tempW(4,0);
        tempW[0]=W[i];
        for (int j=1; j<4; j++)
        {
            tempW[j]=(Wsurr[j-1][i]);
        }
        double Wmax= *max_element(tempW.begin(), tempW.end());
        double Wmin= *min_element(tempW.begin(), tempW.end());
        
        /***Compute the Unconstrained Value at each Gauss Point***/
        vector<double> phig(3,0);
        for (int j=0; j<3; j++)
        {
            double Wgminus = gradW[i][0]*local_vectofaces[j][0] + gradW[i][1]*local_vectofaces[j][1];
            
            double Wgplus=0;
            if(Wgminus > 0)
            {
                Wgplus=Wmax-W[i];
            }
            else if (Wgminus < 0)
            {
                Wgplus=Wmin-W[i];
            }
            phig[j]=limiter(Wgplus,Wgminus,area);
        }
        
        phi[i] = *min_element(phig.begin(), phig.end());
    }
}

/***********************Higher Order Reconstruction Routine********************/
void higherorder(elemvalues elemi, elemvalues elemj, vector<double> &WL, vector<double> &WR, vector< vector<double> >gradWi, vector< vector<double> >gradWj, int limiter, vector< vector< vector<double> > > Usurr, vector< vector< vector<double> > > local_vectofaces, vector< vector<double> > vectoface,vector<double> area)
{
    /*Convenient vector for loops*/
    vector<double> Wi(neqns,0);
    Wi[0]=elemi.rho;
    Wi[1]=elemi.u;
    Wi[2]=elemi.v;
    Wi[3]=elemi.p;
    vector<double> Wj(neqns,0);
    Wj[0]=elemj.rho;
    Wj[1]=elemj.u;
    Wj[2]=elemj.v;
    Wj[3]=elemj.p;
    
    
    vector<double> phiL(neqns,0);
    vector<double> phiR(neqns,0);
    if(limiter == 0)
    {
        for(int i=0; i<neqns; i++)
        {
            phiL[i]=1;
            phiR[i]=1;
        }
    }
    else if (limiter == 2)  //Venkatakrishnan Limiter
    {
        venkatakrishnan(Usurr[0], Wi, local_vectofaces[0], gradWi, phiL,area[0]);
        if(elemj.no <= nelem)
        {
            vector< vector<double> > tempgradWj (neqns, vector<double> (2,0));
            for(int i=0; i<neqns; i++)
            {
                tempgradWj[i][0]=gradWj[i][0];
                tempgradWj[i][1]=gradWj[i][1];
            }
            venkatakrishnan(Usurr[1], Wj, local_vectofaces[1], tempgradWj, phiR,area[1]);
        }
    }
    
    else if (limiter == 1)  //Barth-Jasperson Limiter
    {
        Barth_Jasperson(Usurr[0], Wi, local_vectofaces[0], gradWi, phiL);
        if(elemj.no <= nelem)
        {
            vector< vector<double> > tempgradWj (neqns, vector<double> (2,0));
            for(int i=0; i<neqns; i++)
            {
                tempgradWj[i][0]=gradWj[i][0];
                tempgradWj[i][1]=gradWj[i][1];
            }
            Barth_Jasperson(Usurr[1], Wj, local_vectofaces[1], tempgradWj, phiR);
        }
    }
    
    
    
    
    
    for(int i=0; i<neqns ;i++)
    {
        WL[i]= Wi[i] + phiL[i]*(gradWi[i][0]*vectoface[0][0] + gradWi[i][1]*vectoface[0][1]);
    }
    if(elemj.no <= nelem)
    {
        for(int i=0; i< neqns; i++)
        {
            WR[i]= Wj[i] + phiR[i]*(gradWj[i][0]*vectoface[1][0] + gradWj[i][1]*vectoface[1][1]);
        }
    }
    else
    {
        WR[0]=elemj.rho;
        WR[1]=elemj.u;
        WR[2]=elemj.v;
        WR[3]=elemj.p;
    }
    

}

#endif	/* HIGHERORDERRECONSTRUCT_H */

