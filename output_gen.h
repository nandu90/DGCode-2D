/* 
 * File:   ouput_gen.h
 * Author: nash
 *
 * Created on March 29, 2016, 4:40 PM
 */

#ifndef OUTPUT_GEN_H
#define	OUTPUT_GEN_H


void output(vector< vector<double> > geoel, vector <vector<double> >coord, vector< vector<int> > inpoel, vector < vector<double> > U, int grid, int iter, double secs)
{
    vector< vector<double> > W(nelem, vector<double>(neqns,0));
    vector<double> M(nelem,0);
    double error=0;
    
    /****Convert Conservative variables to primitive****/
    double sum=0;
    for(int i=0; i<nelem; i++)
    {
        cons_to_prim(U[i], W[i][0], W[i][1], W[i][2], W[i][3]);
        double c= sqrt(fabs(g*W[i][3]/W[i][0]));
        M[i] = sqrt(pow(W[i][1],2) + pow(W[i][2],2))/c;
        
        double sgen;
        sgen=((W[i][3]/(pow(W[i][0],g))) - s_inf)/s_inf; 
        sum=sum+pow(sgen,2)*geoel[i][0];
    }
    error=sqrt(sum);
    
    //********Map data from cell centers to nodes*************//
    vector <vector<double> > Wpoin(npoin, vector<double>(neqns,0));
    vector<double> Mpoin(npoin,0);
    vector<double> we(npoin,0);
    //Distance function
    for(int i=0;i<npoin;i++)
    {
        for(int j=0;j<nelem;j++)
        {
            for(int k=0;k<3;k++)
            {
                if(inpoel[j][k]-1 == i)
                {
                    double length=sqrt(pow(geoel[j][1]-coord[inpoel[j][k]-1][0],2)+pow(geoel[j][2]-coord[inpoel[j][k]-1][1],2));
                    for(int l=0;l<neqns;l++)
                    {
                        Wpoin[i][l]=Wpoin[i][l]+W[j][l]*(1/length);
                    }
                    Mpoin[i]=Mpoin[i] + M[j]*(1/length);
                    we[i]=we[i]+(1/length);
                }
            }
        }
    }
    
    for(int i=0;i<npoin;i++)
    {
        for(int j=0;j<neqns;j++)
        {
            Wpoin[i][j]=Wpoin[i][j]/we[i];            
        }
        Mpoin[i]=Mpoin[i]/we[i];
    }
    
    
    string in=dirpath(grid);
    
    //***********Output Generation*********************************//
    ofstream output;
    output.open(getexepath() + in + "/tecplot.dat",ios::trunc);
    output << "title = \"result\" " << endl;
    output << "variables = \"x\", \"y\", \"rho\", \"u\", \"v\", \"p\", \"M\" " << endl;
    output << "zone n=" << npoin << ", e=" << nelem << ", f=fepoint, et=triangle" << endl;
    for(int i=0; i<npoin; i++)
    {
        output << coord[i][0] <<" "<< coord[i][1]<<" "<<Wpoin[i][0]<<" "<<Wpoin[i][1]<<" "<<Wpoin[i][2]<<" "<<Wpoin[i][3]<<" "<<Mpoin[i]<<endl;
    }
    for(int i=0; i<nelem; i++)
    {
        output << inpoel[i][0] <<" "<< inpoel[i][1] <<" "<< inpoel[i][2] <<" "<< endl;
    }
    output.close();
    
    ofstream output2;
    output2.open(getexepath() + in + "/info.txt",ios::trunc);
    output2<<"Total Iterations for Convergence = "<<iter<<endl;
    output2<<"Total time = "<<secs<<" secs"<<endl<<endl<<endl;
    
    output2<<"2-Norm of Error  = "<<error<<endl;
    output2<<"Overall Order = "<<-log10(fabs(error))<<endl;
    output2.close();
}

#endif	/* OUTPUT_GEN_H */

