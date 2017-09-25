/* 
 * File:   gridinput.h
 * Author: User
 *
 * Created on February 2, 2016, 4:43 PM
 */

#ifndef GRIDINPUT_H
#define	GRIDINPUT_H

/*Reads gridfile from the specified case folder*/
void gridcall(vector< vector<double> > &coord, vector< vector<int> > &inpoel, vector< vector<int> > &bface, int grid)
{
    string in=dirpath(grid);
    ifstream gridread;
    gridread.open(getexepath()+in+"/gridfile.txt");
    if (gridread.fail())
    {
        cout<<"Grid file not read. Please check the existence of file"<<endl;
        exit(1);
    }
    
    string line;
    while(getline(gridread,line))
    {
        string junk;
        istringstream ss(line);
        ss>>junk;
        if(junk.compare("ndimn") == 0)
        {
            getline(gridread,line);
            istringstream ss2(line);
            ss2>>ndimn;            
        }
        else if(junk.compare("nelem") == 0)
        {
            getline(gridread,line);
            istringstream ss2(line);
            ss2>>nelem>>npoin>>nface;
            
            inpoel.resize(nelem);
            for(int i=0;i<nelem;i++)
            {
                inpoel[i].resize(3);
            }
            
            coord.resize(npoin);
            for(int i=0;i<npoin;i++)
            {
                coord[i].resize(ndimn);
            }
            
            bface.resize(nface);
            for(int i=0;i<nface;i++)
            {
                bface[i].resize(4);
            }
            goto inpoelread;
        }       
    }
    
    inpoelread:
    getline(gridread,line);
    for(int i=0;i<nelem;i++)
    {
        getline(gridread,line);
        istringstream ss(line);
        int junk;
        ss>>junk>>inpoel[i][0]>>inpoel[i][1]>>inpoel[i][2];
    }
    getline(gridread,line);
    for(int i=0;i<npoin;i++)
    {
        getline(gridread,line);
        istringstream ss(line);
        int junk;
        ss>>junk>>coord[i][0]>>coord[i][1];
    }
    getline(gridread,line);
    for(int i=0;i<npoin;i++)
    {
        getline(gridread,line);
    }
    getline(gridread,line);
    for(int i=0;i<nface;i++)
    {
        getline(gridread,line);
        istringstream ss(line);
        int junk;
        ss>>bface[i][0]>>bface[i][1]>>bface[i][2]>>bface[i][3];
    }
    gridread.close();
}

#endif	/* GRIDINPUT_H */

