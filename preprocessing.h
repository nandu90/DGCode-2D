/* 
 * File:   preprocessing.h
 * Author: User
 *
 * Created on March 25, 2016, 6:16 PM
 */

#ifndef PREPROCESSING_H
#define	PREPROCESSING_H


struct fdata
{
    int ei;
    int ej;
    int pi;
    int pj;
    double nx;
    double ny;
    double len;
};

/***********Calculates elements surrounding point******************************/
void find_elsupo(vector< vector<int> > inpoel, vector<int> &esup1, vector<int> &esup2)
{
    /***Initialize esup2 vector****/    
    esup2.resize(npoin+1);
    for (int i=0;i<npoin+1;i++)
    {
        esup2[i]=0;
    }
    
    /****Find how many elements contain each point****/
    for (int i=0;i<nelem;i++)
    {
        for (int j=0; j<3; j++)  // for triangular elements only
        {
            esup2[inpoel[i][j]]= esup2[inpoel[i][j]]+1;
        }
    }
    /**************************************************/
    
    /****Update Storage Counter*****/
    for(int i=1;i < npoin+1; i++)
    {
        esup2[i]=esup2[i]+esup2[i-1];
        //cout<<esup2[i]<<endl;
    }
    esup2[0]=0;
    /*******************************/
    
    /******Initialize esup1*****/
    esup1.resize(esup2[npoin]+10);
    for(int i=0;i<esup1.size();i++)
    {
        esup1[i]=0;
    }
    
    /****Fill esup1 vector with surrounding element numbers****/
    for (int i=0;i<nelem;i++)
    {
        for(int j=0;j<3;j++)
        {
            int ipoin=inpoel[i][j];
            int istor=esup2[ipoin-1];
            
            esup1[istor]=i+1;
            esup2[ipoin-1]=istor+1;            
        }
    }
    /**********************************************************/
    /***************Resize esup1***********************/
    int size=0;
    for(int i=0;i<esup1.size();i++)
    {
        if(esup1[i] == 0)
        {
            size=i;
            break;
        }
    }
    esup1.resize(size);
    /**********************************************************/
       
    /****Restructure esup2*****/
    for(int i=npoin;i>=1;i--)
    {
        esup2[i]=esup2[i-1];
    }
    esup2[0]=0;
}

/******************Find points surrounding points******************************/
void find_psupo(vector< vector<int> > inpoel, vector<int> esup1, vector<int> esup2, vector<int> &psup1, vector<int> &psup2)
{
    /***Initialize psup2 vector*****/
    psup2.resize(npoin+1);
    for (int i=0;i<psup2.size(); i++)
    {
        psup2[i]=0;
    }
    /********************************/
    
    
    
    int istor=0;
    
    for(int i=0; i<npoin; i++)
    {
        
        /*******Initialize help array*****/
        vector<int> poin1(1,i+1);
        /*Loop over surrounding elements*/
        for(int j=esup2[i]; j<esup2[i+1]; j++)
        {
            int ielem=esup1[j];
            /*Loop over the nodes of the element*/
            for (int k=0;k<3;k++)
            {
                int jpoin=inpoel[ielem-1][k];
                bool match=false;
                for(int l=0;l<poin1.size();l++)
                {
                    if(jpoin == poin1[l])
                    {
                        match=true;
                    }
                }
                if(match == false)
                {
                        psup1.resize(psup1.size()+1);
                        psup1[istor]=jpoin;
                        istor++;
                        
                        poin1.resize(poin1.size()+1);
                        poin1[poin1.size()-1]=jpoin;
                
                }
            }
        }
        
        psup2[i+1]=istor;
    }
    
}


void localfaces(vector< vector<int> > inpoel, vector< vector<int> > &help, int ielem)
{
    for(int j=0; j<3; j++)
        {
            if(j==2)
            {
               help[j][0]=inpoel[ielem][j];
               help[j][1]=inpoel[ielem][0];
            }
            else
            {
               help[j][0]=inpoel[ielem][j];
               help[j][1]=inpoel[ielem][j+1]; 
            }
        }
}

/******************Find elements surrounding elements**************************/
void find_esuel(vector< vector<int> > inpoel, vector<int> esup1, vector<int> esup2, vector<vector <int> > &esuel)
{
    for(int i=0;i<nelem;i++)
    {
        /***Construct Local faces***/
        vector< vector<int> > help(3, vector<int> (2,0));        
        localfaces(inpoel,help,i);
        /**********************/
        
        /**Loop over element faces**/
        for(int j=0;j<3;j++)
        {
            int ip1=help[j][0]-1;
            int ip2=help[j][1]-1;
            
            /*Loop over the elements surrounding ip1*/
            for(int istor=esup2[ip1]; istor<esup2[ip1+1]; istor++)
            {
                int jelem=esup1[istor];
                
                if(jelem != i+1)
                {
                    /*Loop over elements surrounding ip2*/
                    for(int istor2=esup2[ip2]; istor2<esup2[ip2+1]; istor2++)
                    {
                        int jelem2=esup1[istor2];
                        if (jelem == jelem2)
                        {
                            esuel[i][j]=jelem2;
                            goto next;
                        }
                    }
                }
            }
            next:;
        }
    }
}

/************************Find and store face data******************************/
void find_faces(vector< vector<int> > inpoel, vector< vector<int> > &esuel, vector<fdata> &intfac, int &nbfac, int &nafac, vector< vector<double> > coord)
{
    nbfac=0;
    
    /*****Get the boundary faces****/
    for(int i=0; i<nelem; i++)
    {
        vector< vector<int> > localfac(3, vector<int> (2,0));        
        localfaces(inpoel,localfac,i);
        for(int j=0; j<3; j++)
        {
            if(esuel[i][j] == 0)
            {
                intfac.resize(intfac.size()+1);
                intfac[nbfac].ei=i+1;
                intfac[nbfac].ej=nelem+nbfac+1;
                intfac[nbfac].pi=localfac[j][0];
                intfac[nbfac].pj=localfac[j][1];
                esuel[i][j]=nelem+nbfac+1;
                
                intfac[nbfac].nx=coord[intfac[nbfac].pj-1][1] - coord[intfac[nbfac].pi-1][1];
                intfac[nbfac].ny=-(coord[intfac[nbfac].pj-1][0] - coord[intfac[nbfac].pi-1][0]);
                intfac[nbfac].len=sqrt(pow(intfac[nbfac].nx,2)+ pow(intfac[nbfac].ny,2));
                
                intfac[nbfac].nx=intfac[nbfac].nx/intfac[nbfac].len;
                intfac[nbfac].ny=intfac[nbfac].ny/intfac[nbfac].len;
                
                nbfac++;
            }
        }
    }
    /****************************/
    nafac=nbfac;
    /*****Get the internal faces*****/
    for(int i=0;i<nelem;i++)
    {
        vector< vector<int> > localfac(3, vector<int> (2,0));        
        localfaces(inpoel,localfac,i);
        for(int j=0;j<3;j++)
        {
            int je=esuel[i][j];
            
            if(je > i+1 && je < nelem+1)
            {
                intfac.resize(intfac.size()+1);
                intfac[nafac].ei=i+1;
                intfac[nafac].ej=je;
                intfac[nafac].pi=localfac[j][0];
                intfac[nafac].pj=localfac[j][1];
                
                
                intfac[nafac].nx=coord[intfac[nafac].pj-1][1] - coord[intfac[nafac].pi-1][1];
                intfac[nafac].ny=-(coord[intfac[nafac].pj-1][0] - coord[intfac[nafac].pi-1][0]);
                intfac[nafac].len=sqrt(pow(intfac[nafac].nx,2)+ pow(intfac[nafac].ny,2));
                
                intfac[nafac].nx=intfac[nafac].nx/intfac[nafac].len;
                intfac[nafac].ny=intfac[nafac].ny/intfac[nafac].len;
                
                nafac++;
            }
        }
    }
}

void preprocess(vector< vector<double> > coord, vector< vector<int> > inpoel, vector< vector<int> > bface, vector<vector<double> > &geoel, vector<fdata> &intfac, vector< vector <int> > &esuel, int &nbfac, int &nafac)
{
    
    /***Find elements surrounding points****/
    vector<int> esup1;
    vector<int> esup2;
    find_elsupo(inpoel, esup1,esup2);
    /***************************************/
    
    /***Find Points surrounding points*******/
    vector<int> psup1;
    vector<int> psup2;
    find_psupo(inpoel, esup1, esup2, psup1, psup2);
    /***************************************/
    
    /****Find Elements surrounding elements****/
    find_esuel(inpoel,esup1,esup2,esuel);
    /******************************************/
    
    /***Get the face information***/
    find_faces(inpoel,esuel,intfac,nbfac,nafac,coord);
    
    
        
    
    /****Get the area of the elements - Mass Matrix*/
    for(int i=0; i < nelem; i++)
    {
        double a1,a2,b1,b2;        
        a1=coord[inpoel[i][1]-1][1]-coord[inpoel[i][2]-1][1];
        a2=coord[inpoel[i][2]-1][1]-coord[inpoel[i][0]-1][1];
        
        b1=coord[inpoel[i][2]-1][0]-coord[inpoel[i][1]-1][0];
        b2=coord[inpoel[i][0]-1][0]-coord[inpoel[i][2]-1][0];
        
        double D=fabs(a1*b2-a2*b1);
        
        geoel[i][0]=D/2;
        
        /*********Also store the coordinates of centroid of elements*****/
        geoel[i][1]=(coord[inpoel[i][0]-1][0]+coord[inpoel[i][1]-1][0]+coord[inpoel[i][2]-1][0])/3;  //xc
        geoel[i][2]=(coord[inpoel[i][0]-1][1]+coord[inpoel[i][1]-1][1]+coord[inpoel[i][2]-1][1])/3;  //yc
    }
    
}


struct ghostdata
{
    int elem;
    double cenx;
    double ceny;
};


void ghostcentroid(vector<ghostdata> &geogh, vector<fdata> intfac, int nbfac, vector< vector<double> > coord, vector< vector<double> > geoel)
{
    geogh.resize(nbfac);  //No of ghost cells = number of boundary faces
    
    for(int i=0; i<nbfac; i++)
    {
        geogh[i].elem=intfac[i].ej; //Boundary face data is stored first in intfac
        
        /*Find the shortest distance of centroid of inner cell from face*/
        double x1=coord[intfac[i].pi-1][0];
        double y1=coord[intfac[i].pi-1][1];
        
        double x2=coord[intfac[i].pj-1][0];
        double y2=coord[intfac[i].pj-1][1];
        
        double xc=geoel[intfac[i].ei-1][1];
        double yc=geoel[intfac[i].ei-1][2];
        
        double d=fabs((x2-x1)*(y1-yc)-(x1-xc)*(y2-y1))/sqrt(pow(x2-x1,2)+pow(y2-y1,2));
        
        geogh[i].cenx=xc+2*d*intfac[i].nx;
        geogh[i].ceny=yc+2*d*intfac[i].ny;
        
    }
    
}

void lsquare_matrix(vector< vector <int> > esuel, vector< vector<double> > geoel, vector< vector<double> > &lsquare, vector<ghostdata> geogh)
{
    #pragma omp parallel for schedule (dynamic)
    for (int i=0; i<nelem; i++)
    {
        double xi=geoel[i][1]; //Centroid coordinates
        double yi=geoel[i][2];
        double a=0;
        double b=0;
        double c=0;
        
        for(int j=0; j<3; j++)
        {
            double xj,yj;
            if(esuel[i][j] > nelem)
            {
                for(int k=0; k<geogh.size(); k++)
                {
                    if(esuel[i][j] == geogh[k].elem)
                    {
                        xj=geogh[k].cenx;  //Centroid coordinates
                        yj=geogh[k].ceny;
                        break;
                    }
                }
            }
            else
            {
                xj=geoel[esuel[i][j]-1][1]; //Centroid coordinates
                yj=geoel[esuel[i][j]-1][2];
            }
            
            double w=1/sqrt(pow(xj-xi,2) + pow(yj-yi,2));
            a=a+pow(w,2)*pow(xj-xi,2);
            b=b+pow(w,2)*(xj-xi)*(yj-yi);
            c=c+pow(w,2)*pow(yj-yi,2);
        }
        lsquare[i][0]=a;
        lsquare[i][1]=b;
        lsquare[i][2]=c;
    }
}



void elem_vectofaces(vector< vector< vector<double> > > &vec_to_faces, vector<fdata> intfac, vector< vector<double> > coord, vector< vector<int> > esuel, vector< vector<double> > geoel)
{
    #pragma omp parallel for schedule (dynamic)
    for (int i=0; i<nelem; i++)
    {
        for (int j=0; j<3; j++)
        {
            for(int k=0; k<intfac.size(); k++)
            {
                if(i+1 == intfac[k].ei && esuel[i][j] == intfac[k].ej || esuel[i][j] == intfac[k].ei && i+1 == intfac[k].ej)
                {
                    double posx=(coord[intfac[k].pi-1][0]+coord[intfac[k].pj-1][0])*0.5;
                    double posy=(coord[intfac[k].pi-1][1]+coord[intfac[k].pj-1][1])*0.5;
                    //cout<<i+1<<" "<<j+1<<" "<<posx<<" "<<posy<<" "<<geoel[i][1]<<" "<<geoel[i][2]<<endl;
                    vec_to_faces[i][j][0]=posx-geoel[i][1];
                    vec_to_faces[i][j][1]=posy-geoel[i][2];
                    
                    break;
                }
            }
        }
    }
}
#endif	/* PREPROCESSING_H */

