/* 
 * File:   init_fields.h
 * Author: nash
 *
 * Created on March 29, 2016, 2:52 PM
 */

#ifndef INIT_FIELDS_H
#define	INIT_FIELDS_H
#include "functions.h"

/**************Function to initialize flow Field*******************************/
void initialize(vector< vector<double> > &W, vector< vector<double> > &U)
{
    for (int i=0;i<nelem;i++)
    {
        W[i][0]=rho_inf;
        W[i][1]=v_inf*cos(alpha);
        W[i][2]=v_inf*sin(alpha);
        W[i][3]=pow(v_inf,2)*rho_inf/(pow(M_inf,2)*g);
    }
    s_inf=W[0][3]/(pow(rho_inf,g));
    prim_to_cons(W,U);
}

/******************Function to update Ghost cell values************************/
void bound_update(vector< vector<int> > *bface, vector<fdata> *intfac, vector< vector<double> > *U)
{
    #pragma omp parallel for schedule (dynamic)
    for(int i=0; i<(*bface).size(); i++)
    {
        for(int j=nelem; j< (*U).size(); j++)
        {
            if(j+1 == (*bface)[i][0])
            {
                for(int k=0; k<(*bface).size(); k++)
                {
                    if((*bface)[i][1] == (*intfac)[k].pi && (*bface)[i][2] == (*intfac)[k].pj)
                    {
                        int selem=(*intfac)[k].ei;
                        if((*bface)[i][3] == 2 ) //wall BC (mirror) (normal component of vel at boundary must be 0)
                        {
                            if(j+1 != (*intfac)[k].ej)
                            {
                                cout<<"error"<<endl;
                                exit(1);
                            }
                            (*U)[j][0]=(*U)[selem-1][0];
                            (*U)[j][3]=(*U)[selem-1][3];
                            
                            double nx=(*intfac)[k].nx;
                            double ny=(*intfac)[k].ny;
                            double vn=(*U)[selem-1][1]*nx+(*U)[selem-1][2]*ny;
                            
                            
                            (*U)[j][1]=(*U)[selem-1][1] - 2*vn*nx;
                            (*U)[j][2]=(*U)[selem-1][2] - 2*vn*ny;
                            
                        }
                        
                        else if ((*bface)[i][3] == 4) //In/Out Flow   For now, set everything to free stream values
                        {
                            (*U)[j][0]=rho_inf;
                            (*U)[j][1]=rho_inf*v_inf*cos(alpha);
                            (*U)[j][2]=rho_inf*v_inf*sin(alpha);
                            double p=pow(v_inf,2)*rho_inf/(pow(M_inf,2)*g);
                            (*U)[j][3]=p/(g-1) + 0.5*(pow((*U)[j][1],2) + pow((*U)[j][2],2))/(*U)[j][0];
                        }
                        break;
                    }
                }
                break;
            }
        }
    }
}

/******************Assigns element numbers to Ghost cells**********************/
void boundary(vector<fdata> intfac, vector< vector<int> > &bface, vector< vector<double> > &U)
{
    int nbfac=bface.size();   
    
    /***Now resize the U matrix to accomodate ghost elements in the end****/
    U.resize(U.size()+nbfac);  //No of boundary elements = number of boundary faces
    for(int i=nelem; i<U.size(); i++)
    {
        U[i].resize(4,0);
    }
    
    
    for(int i=nelem; i<U.size(); i++) //start loop from end of internal element index i.e. nelem
    {
        for(int j=0; j<nbfac; j++)
        {
            if(i+1 == intfac[j].ej) //intfac is structured so that always the second element is boundary element
            {
                for(int k=0; k<nbfac; k++)
                {
                    if((bface[k][1] == intfac[j].pi && bface[k][2] == intfac[j].pj))
                    {
                        bface[k][0] =i+1;  //overwrite the element numbers in bface with new element numbers for boundary faces
                        break;
                    }
                }
                break;
            }
        }
    }
    
    bound_update(&bface, &intfac, &U);   
    
}

#endif	/* INIT_FIELDS_H */

