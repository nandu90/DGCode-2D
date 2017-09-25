/* 
 * File:   functions.h
 * Author: nash
 *
 * Created on March 29, 2016, 4:31 PM
 */

#ifndef FUNCTIONS_H
#define	FUNCTIONS_H

/*Converts primitive variables to conservative*/
void prim_to_cons(vector< vector<double> > W, vector <vector<double> > &U)
{
    for (int i=0;i<nelem;i++)
    {
        U[i][0]=W[i][0];
        U[i][1]=W[i][1]*W[i][0];
        U[i][2]=W[i][2]*W[i][0];
        U[i][3]=W[i][3]/(g-1) + 0.5*W[i][0]*(pow(W[i][1],2) + pow(W[i][2],2));
    }
}
/*Converts Conservative variables to primitive*/
void cons_to_prim(vector<double>  U, double &rho, double &u, double &v, double &p)
{
    rho=U[0];
    u=U[1]/U[0];
    v=U[2]/U[0];
    p=(g-1)*(U[3]-0.5*(pow(U[1],2) + pow(U[2],2))/U[0]);
}

/*Calculates unit normal vectors to faces and face lengths*/
void face_normal(int e1, int e2, vector<fdata> *intfac, double &nx, double &ny, double &len)
{
    for(int i=0; i<(*intfac).size();i++)
    {
        if((*intfac)[i].ei == e1 && (*intfac)[i].ej == e2)
        {
            nx=(*intfac)[i].nx;
            ny=(*intfac)[i].ny;
            len=(*intfac)[i].len;
            break;
        }
        else if((*intfac)[i].ei == e2 && (*intfac)[i].ej == e1)
        {
            nx=-(*intfac)[i].nx;
            ny=-(*intfac)[i].ny;
            len=(*intfac)[i].len;
            break;
        }
    }
}

/*Calculates global and local time step based on CFL*/
void time_step(vector <vector<double> > *U, vector <vector<double> > *geoel, vector<fdata> *intfac, vector< vector <int> > *esuel, double &deltat, vector<double> *dt)
{
       
    #pragma omp parallel for schedule (dynamic)
    for(int i=0;i<nelem;i++)
    {
        double rhoi,ui,vi,pi;
        cons_to_prim((*U)[i],rhoi,ui,vi,pi);       
        double ci=sqrt(fabs(g*pi/rhoi));
        double surfint=0;
        
        for(int j=0; j<3; j++)
        {
            int selem=(*esuel)[i][j];
            double nx,ny,length;
            face_normal(i+1,selem,intfac,nx,ny,length);
            
            double rhoj,uj,vj,pj;
            cons_to_prim((*U)[selem-1],rhoj,uj,vj,pj);            
            double cj=sqrt(fabs(g*pj/rhoj));
            
            surfint=surfint + length*0.5*(ci+cj + fabs((ui+uj)*nx + (vi+vj)*ny));
            
        }
        
        (*dt)[i]=cfl*(*geoel)[i][0]/surfint;
    }
    
    deltat=(*dt)[0];
    for(int i=1; i<nelem; i++)
    {
        if((*dt)[i] <= deltat)
        {
            deltat=(*dt)[i];
        }
    }
}

/*Prints variables in a file for each iteration*/
void monitor_cons(vector< vector<double> > *U, vector< vector<double> > *geoel, int t, bool &exitflag, int grid)
{
    vector<double> totalU (neqns,0);
    double area=0;
    for(int i=0; i< nelem; i++ )
    {
        for(int j=0; j<neqns; j++)
        {
           totalU[j] = totalU[j] + (*U)[i][j]*(*geoel)[i][0];
        }
        area=area+(*geoel)[i][0];
    }
    string in=dirpath(grid);
    ofstream out;
    out.open(getexepath() +in + "/conservation.txt",ios::app);
    out<<t<<" "<<totalU[0]/area<<" "<<totalU[1]/area<<" "<<totalU[2]/area<<" "<<totalU[3]/area<<endl;
    out.close();
}

/*Prints variable residuals in a file for each iteration*/
void monitor_res (vector< vector<double> > *Ui, vector< vector<double> > *U, vector< vector<double> > *geoel, int t, bool &exitflag, vector <double> &ires, int grid)
{
    vector<double> res (neqns,0);
    for(int i=0; i< nelem; i++ )
    {
        for(int j=0; j<neqns; j++)
        {
            res[j]=res[j] + pow((*Ui)[i][j] - (*U)[i][j],2)*(*geoel)[i][0];
        }
    }
    
    for(int j=0; j<neqns; j++)
    {
        res[j]=sqrt(res[j]);
    }
    
    
    
    if(t==1)
    {
        for(int j=0; j<neqns; j++)
        {
            ires[j]=res[j];
        }
    }
    else
    {
        string in=dirpath(grid);
        ofstream out;
        out.open(getexepath() + in + "/residuals.txt",ios::app);
        out<<t<<" "<<res[0]/ires[0]<<" "<<res[1]/ires[1]<<" "<<res[2]/ires[2]<<" "<<res[3]/ires[3]<<endl;
        out.close();
        cout<<"Step: "<<t<<" Mass: "<<res[0]/ires[0]<<" xMom: "<<res[1]/ires[1]<<" yMom: "<<res[2]/ires[2]<<" Energy: "<<res[3]/ires[3]<<endl;
        
        if(res[0]/ires[0] < 1e-4 && res[1]/ires[1] < 1e-4 && res[2]/ires[2] < 1e-4 && res[3]/ires[3] < 1e-3)
       {
           exitflag=true;
       }
    }
}


/***Function for Time integration by forward Euler Method***/
void euler (vector< vector<double> > *U, vector< vector<double> > *geoel, vector< vector<double> > *R, vector<double>* dt)
{
    #pragma omp parallel for schedule (dynamic)
    for ( int i=0;i<nelem;i++)
    {
        for (int j=0;j<neqns;j++)
        {
                (*U)[i][j] = (*U)[i][j] + (*dt)[i]*(*R)[i][j]/(*geoel)[i][0];
        
        }
    }
}


#endif	/* FUNCTIONS_H */

