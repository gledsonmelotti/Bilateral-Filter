/*************************************************************************/
/* Upsampling from dense 3D LIDAR                                        */
/*                                                                       */
/* C.Premebida: June/2014                                                */
/* http://webmail.isr.uc.pt/~cpremebida/IROS14/LaserVisionFusion.html    */
/*                                                                       */
/*                                                                       */
/*                                                                       */
/*                                                                       */
/* All rights reserved.                                                  */
/*************************************************************************/



 

#include "mex.h"
#include <iostream>
#include <cmath>
#include "matrix.h"
#include <list>
#include <iostream>
#include <algorithm>
#include <vector>
 

//#if defined(__WIN32__) || defined(__WIN64__)

int roundr(double x)
{
    return x >= 0.0f ? floor(x + 0.5f) : ceil(x - 0.5f);
}
//#endif

class Point{
public:
        double x;
        double y;
        double z;
        double i;
        
        Point(){
            x=0;
            y=0;
            z=0;
            i=0;
        }
        
        Point(double x_,double y_,double z_,double i_){
            x=x_;
            y=y_;
            z=z_;
            i=i_;
        }
        
        
        bool xyEqual(const Point &pt)const {
            
            return (roundr(x*100.0)==roundr(pt.x*100.0)) && (roundr(y*100.0)==roundr(pt.y*100.0));
            
        }
        
        double distance(Point &a){
            return sqrt((x-a.x)*(x-a.x)+(y-a.y)*(y-a.y)+(z-a.z)*(z-a.z));
        }
        
        double distance2D(Point &a){
            return sqrt((x-a.x)*(x-a.x)+(y-a.y)*(y-a.y));
        }
        
//          bool operator() (Point &i,Point &j) const{
//              return (i.z<j.z);
//          }
//         bool operator< (Point &j) const {
//             return (z<j.z);
//         }
        
        static bool sorter (Point &g,Point &j)  {
            return (g.z<j.z);
        }
        
};





static bool clustersorter (std::pair<int,std::pair<int,double> > &g,std::pair<int,std::pair<int,double> > &j)  {
    
    if(g.first==j.first){
        return (g.second.second<j.second.second);
    }
    return (g.first>j.first);
}



static double clusterzmin (std::vector<Point> &vect)  {
    
    double out=99999;
    
    for(int i=0;i<vect.size();i++){
        
        
        if(vect[i].z<out){
            out=vect[i].z;
        }
    }
    
    return out;
}

static double clusternear (std::vector<Point> &vect,int u,int v)  {
    
    double out=0;
    double dist=99999;
    
    for(int i=0;i<vect.size();i++){
        
        double temp=sqrt((u-vect[i].x)*(u-vect[i].x)+(v-vect[i].y)*(v-vect[i].y));
        if(temp<dist){
            dist=temp;
            out=vect[i].z;
        }
    }
    
    return out;
}

 



std::vector<Point> operator*(const std::vector<Point>& v, const std::vector<Point>& vl)
{
    std::vector<Point>  out;
    Point pout;
    
    for(int i=0;i<v.size();i++){
        for(int g=0;g<vl.size();g++){
            
            if(v[i].xyEqual(vl[g])){
                
                pout=v[i];
                pout.z=v[i].z*vl[g].z;
                out.push_back(pout);
                continue;
            }
        }
    }
    
    
    return out;
    
}

template <typename T>
        class Mat
{
    unsigned int index;
    
        public:
            
            
            unsigned int rows;
            unsigned int columns;
            unsigned int count;
            T *matrix;
            
            Mat(unsigned int rows_, unsigned int columns_){
                rows=rows_;
                columns=columns_;
                count=rows_*columns_;
                matrix =new T [count];
            }
            
            T get(unsigned int row,unsigned int column) const{
                if((column>=columns)||(row>=rows)){
                    return 0;
                }else{
                    return matrix[column*rows+row];
                }
            }
            
            
            
            void set(unsigned int row,unsigned int column,T value) const{
                
                if((column>=columns)||(row>=rows)){
                    //  return std::vector<Point>();
                }else{
                    matrix[column*rows+row]=(value);
                }
            }
            
            ~Mat(){
                delete [] matrix;
            }
            
};

class Matrix
{
    unsigned int index;
    
public:
    
    
    unsigned int rows;
    unsigned int columns;
    unsigned int count;
    std::list<int> *matrix;
    
    Matrix(unsigned int rows_, unsigned int columns_){
        rows=rows_;
        columns=columns_;
        count=rows_*columns_;
        matrix =new std::list<int> [count];
    }
    
    std::list<int> get(unsigned int row,unsigned int column) const{
        if((column>=columns)||(row>=rows)){
            return std::list<int>();
        }else{
            return matrix[column*rows+row];
        }
    }
    
    std::list<int> getUnique(unsigned int row,unsigned int column,double *x, int sd) const{
        if((column>=columns)||(row>=rows)){
            return std::list<int>();
        }else{
            
            auto dal1=matrix[column*rows+row];
            if(dal1.size()==0){
                return std::list<int>();
            }
            
            std::vector<int> dal;
            
            for(auto it=dal1.begin();it!=dal1.end();it++)
                dal.push_back(*it);
            
            
            
            std::vector<double> da;
            da.reserve(dal.size());
            for(int i=0;i<dal.size();i++ ){
                
                da.push_back(x[dal[i]+sd]);
                
            }
            
            
            
            
            
            auto vmin=std::min_element(std::begin(da), std::end(da));
            
            int index=std::distance(std::begin(da), vmin);
            std::list<int> out;
//                 auto iter = dal.begin();
//                 std::advance(iter, index);
//                 int x = *iter;
//
            
            out.push_back(dal[index]);
            return out;
        }
    }
    
    
    
    
    
    void add(unsigned int row,unsigned int column,int value) const{
        
        if((column>=columns)||(row>=rows)){
            //  return std::vector<Point>();
        }else{
            matrix[column*rows+row].push_back(value);
        }
    }
    
    ~Matrix(){
        delete [] matrix;
    }
    
};


double sum(std::vector<Point> pt){
    
    double val=0;
    for(unsigned int i=0;i<pt.size();i++){
        
        val+=pt[i].z;
    }
    return val;
}




double distance(double x0,double y0,double z0,double x1,double y1,double z1){
    return sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0));
}

double distance(double x0,double y0,double x1,double y1){
    return sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
}
/****************************/
/* The computational routine */
//sd number of points in pointcloud
void calc_Dense(double *x, double *y, double *dim, int sd,int IDMETHOD,int grid_,double threshold, double epsstep)
{
    
    
//     mexPrintf("method: %d grid: %d\n",IDMETHOD,grid_);
//     mexEvalString("disp('--------')");
    
    mxArray *rhs[1];
    rhs[0]= mxCreateDoubleMatrix(grid_, grid_, mxREAL);
    double *patch=mxGetPr(rhs[0]);
    
    for(int i=0;i<grid_*grid_;i++){
        patch[i]=0.0;
    }
    
//     mexEvalString("disp('--------')");
    
    int globalgrid=grid_; // Size of the mask/kernel
    int grid=globalgrid;
    int maxgrid=grid_;
    
    
    Matrix matrix((int)dim[2]+grid,(int)dim[1]+grid);
    Mat<double> astore((int)dim[2],(int)dim[1]);
    Mat<double> bstore((int)dim[2],(int)dim[1]);
    
    
    int kin=0;
    double mr = 0;
    
    double Gs=0, Gr=0;
    double S=0, Y=0, d=0, WGain=0;
    
    /*******************/
    
    
    for (int k=0; k<sd; k=k+1){
        
        
        // Point pt;
        // pt.x=x[(*iterator)+0*sd];
        // pt.y=x[k+1*sd];
        // pt.z=x[k+2*sd];
        // pt.i=x[k+3*sd];
        int id1x=roundr(x[k+0*sd])-1;
        int id1y=roundr(x[k+1*sd]-(int)dim[0]);
        
        // int id1x=floor(x[k+0*sd]);
        // int id1y=floor(x[k+1*sd]-dim[0]);
        //
        // int id2x=id1x;
        // int id2y=id1y+1;
        //
        // int id3x=id1x;
        // int id3y=id1y-1;
        //
        // int id4x=id1x-1;
        // int id4y=id1y+1;
        //
        // int id5x=id1x-1;
        // int id5y=id1y;
        //
        // int id6x=id1x-1;
        // int id6y=id1y-1;
        //
        // int id7x=id1x+1;
        // int id7y=id1y+1;
        //
        // int id8x=id1x+1;
        // int id8y=id1y;
        //
        // int id9x=id1x+1;
        // int id9y=id1y-1;
        
        
        matrix.add(id1x,id1y,k);
        
        
        // if(distance(id1x,id1y,x[k],x[k+sd]-dim[0])<0.65){
        //     matrix.add(id1x,id1y,k);
        // }
        //
        // if(distance(id2x,id2y,x[k],x[k+sd]-dim[0])<0.65){
        //     matrix.add(id2x,id2y,k);
        // }
        //
        // if(distance(id3x,id3y,x[k],x[k+sd]-dim[0])<0.65){
        //     matrix.add(id3x,id3y,k);
        // }
        //
        // if(distance(id4x,id4y,x[k],x[k+sd]-dim[0])<0.65){
        //     matrix.add(id4x,id4y,k);
        // }
        //
        // if(distance(id5x,id5y,x[k],x[k+sd]-dim[0])<1.65){
        //     matrix.add(id5x,id5y,k);
        // }
        //
        // if(distance(id6x,id6y,x[k],x[k+sd]-dim[0])<0.65){
        //     matrix.add(id6x,id6y,k);
        // }
        //
        // if(distance(id7x,id7y,x[k],x[k+sd]-dim[0])<0.65){
        //     matrix.add(id7x,id7y,k);
        // }
        //
        // if(distance(id8x,id8y,x[k],x[k+sd]-dim[0])<0.65){
        //     matrix.add(id8x,id8y,k);
        // }
        //
        // if(distance(id9x,id9y,x[k],x[k+sd]-dim[0])<0.65){
        //     matrix.add(id9x,id9y,k);
        // }
        
        
        
        //   int x1=floor(x[k+0*sd]);
        //   int y1=ceil(x[k+1*sd]-dim[0]);
        //
        //   int x2=ceil(x[k+0*sd]);
        //   int y2=floor(x[k+1*sd]-dim[0]);
        //
        //   int x3=floor(x[k+0*sd]);
        //   int y3=floor(x[k+1*sd]-dim[0]);
        //
        //   int x4=ceil(x[k+0*sd]);
        //   int y4=ceil(x[k+1*sd]-dim[0]);
        //
        //
        //
        //     matrix.add(roundr(x[k+0*sd]),roundr(x[k+1*sd]-dim[0]),k);
        //     matrix.add(x2,y2,k);
        //     matrix.add(x3,y3,k);
        //     matrix.add(x4,y4,k);
        
        
        
        //matrix.add(roundr(x[k+0*sd]),roundr(x[k+1*sd]-dim[0]),k);
    }
    
    
    
    
    // // compute maximum in z???
    for (  int k=0; k<sd; k=k+1){
        
        
        if  (x[k+2*sd] > mr) {
            
            mr = x[k+2*sd];
        }
    }
    
    
    
    
    
    /*******************/
    
    int griddisx=floor(grid/2);
    int griddisy=floor(grid/2);
    
    
    
//     mexEvalString("disp('--------')");
    
    // column
    for (int u=0; u<(int)dim[2]; u=u+1)
    {
        // line
        for (int v=0; v<(int)dim[1]; v=v+1)
        {
            // for each column... up down...
            
            
            bool changed =false;
            grid=globalgrid;
            while(grid<=maxgrid){
                //griddisx=floor(grid/2);2
                griddisx=floor(grid/2);
                griddisy=floor(grid/2);
                
                S=0; Y=0; // 9999 on min
                
                std::vector<int> pointscell;
                std::vector<int> upointscell;
                std::vector<int> vpointscell;
                
                // x
                
                if(IDMETHOD==13 || IDMETHOD==14){
                    
                    for (int i=(u-griddisx);i<=(u+griddisx);i++){
                        // y
                        for (int f=(v-griddisy);f<=(v+griddisy);f++){
                            
                            std::list<int> temp =matrix.getUnique(i,f,x,2*sd);
                            //std::list<int> temp =matrix.get(i,f);
                            upointscell.push_back(i);
                            vpointscell.push_back(f);
                            
                            for(auto it=temp.begin();it!=temp.end();it++)
                                pointscell.push_back(*it);
                            
                        }
                        
                    }
                    
                    
                }else{
                    
                    for (int i=(u-griddisx);i<=(u+griddisx);i++){
                        // y
                        for (int f=(v-griddisy);f<=(v+griddisy);f++){
                            
                            // std::list<int> temp =matrix.getUnique(i,f,x,2*sd);
                            std::list<int> temp =matrix.get(i,f);
                            
                            for(auto it=temp.begin();it!=temp.end();it++)
                                pointscell.push_back(*it);
                            
                        }
                        
                    }
                    
                    
                    
                }
                
                
                //////////////////
                //////////////////
                //////////////////
                //////////////////
                //////////////////
                
                // counter
                 if(IDMETHOD==-2){
                    
                    Y=0;
                    S=1;
                    
                    for (unsigned int i=0;i<pointscell.size();i++){
                        Y=std::max(Y,x[(pointscell[i])+2*sd]);
                        changed=true;
                    }
                    
                    //mean
                }else if(IDMETHOD==-1){
                    
                    S=1;
                    
                    Y=pointscell.size();
                    changed=true;
                    
                    //mean
                }else if(IDMETHOD==0){
                    
                    Y=999;
                    S=1;
                    
                    for (unsigned int i=0;i<pointscell.size();i++){
                        Y=std::min(Y,x[(pointscell[i])+2*sd]);
                        changed=true;
                    }
                    
                    //mean
                }else if(IDMETHOD==1){
                    
                    Y=0;
                    S=0;
                    
                    for (unsigned int i=0;i<pointscell.size();i++){
                        Y=Y+x[(pointscell[i])+2*sd];
                        S=S+1;
                        changed=true;
                    }
                    if(changed){
                        
                        Y=Y/S;
                    }
                    
                    //median
                }else if(IDMETHOD==2)
                {
                    
                    if(pointscell.size()==0){
                        changed=false;
                    }else if(pointscell.size()==1){
                        changed=true;
                        Y=x[(pointscell[0])+2*sd];
                    }else if(pointscell.size()%2==0){
                        
                        //std::size_t middleIdx = size/2;
                        //RandAccessIter target = begin + middleIdx;
                        //int idx1=std::nth_element(begin, target, end);
                        
                        Y=(x[(pointscell[pointscell.size()/2-1])+2*sd] + x[(pointscell[pointscell.size()/2])+2*sd]) / 2.0;
                        changed=true;
                        
                    }else{
                        sort(pointscell.begin(), pointscell.end());
                        changed=true;
                        
                        Y=x[(pointscell[pointscell.size()/2])+2*sd];
                    }
                    
                    //IDW
                }else if(IDMETHOD==3)
                {
                    
                    for (unsigned int i=0;i<pointscell.size();i++){
                        
                        
                        d =  (u - x[(pointscell[i])])*(u - x[(pointscell[i])]) + (v+dim[0]-x[(pointscell[i])+sd])*(v+dim[0]-x[(pointscell[i])+sd]);
                        
                        WGain = 1.0/sqrt(d);
                        
                        S = S + WGain;
                        Y = Y + WGain*(x[(pointscell[i])+2*sd]);
                        changed=true;
                        
                    }
                    
                    if(changed==true){
                        
                        Y=Y/S;
                        
                    }
                    //BF
                }else if(IDMETHOD==4)
                {
                    
                    for (unsigned int i=0;i<pointscell.size();i++){
                        
                        
                        Gr = x[(pointscell[i])+2*sd]/mr;
                        Gs =  ((u - x[(pointscell[i])])*(u - x[(pointscell[i])]) + (v+dim[0]-x[(pointscell[i])+sd])*(v+dim[0]-x[(pointscell[i])+sd]) );
                        WGain = 1.0/sqrt(Gs*Gr);
                        S = S + WGain;
                        Y = Y + WGain*(x[(pointscell[i])+2*sd]);
                        changed=true;
                        
                    }
                    
                    if(changed==true){
                        
                        Y=Y/S;
                        
                    }
                    
                    //HMF
                }else if(IDMETHOD==5)
                {
                    
                    
                    Y=0;
                    S=0;
                    
                    for (unsigned int i=0;i<pointscell.size();i++){
                        Y=Y+x[(pointscell[i])+2*sd];
                        S=S+1;
                        changed=true;
                    }
                    //                      if(changed){
                    //
                    //                         Y=Y/S;
                    //                      }
                    
                    
                    if(changed){
                        
                        
                        double mean1=Y/S;
                        S=0;
                        Y=0;
                        
                        for (unsigned int i=0;i<pointscell.size();i++) {
                            
                            //mean
                            S=S+1;
                            
                            Y=Y+(x[(pointscell[i])+2*sd]-mean1)*(x[(pointscell[i])+2*sd]-mean1);
                            
                        }
                        
                        double var=Y/S;
                        // here starts the more interesting part... the backgroundr contains a higher variance... so, while
                        // it is possible to use mean to solve the problem in most cases, situations where variance is higher than X
                        // a new solution should be employed..
                        
                        
                        
                        // just to output variance
                        //Y=var;
                        if(var>0.25){
                            
                            if(pointscell.size()%2==0){
                                
                                //std::size_t middleIdx = size/2;
                                //RandAccessIter target = begin + middleIdx;
                                //int idx1=std::nth_element(begin, target, end);
                                
                                Y=(x[(pointscell[pointscell.size()/2-1])+2*sd] + x[(pointscell[pointscell.size()/2])+2*sd]) / 2;
                                changed=true;
                                
                            }else{
                                sort(pointscell.begin(), pointscell.end());
                                changed=true;
                                
                                Y=x[(pointscell[pointscell.size()/2])+2*sd];
                            }
                            
                            
                        }else{
                            Y=mean1;
                        }
                        
                        
                        
                        
                    }
                    
                }else if(IDMETHOD==6)
                {
                    //local max
                    mr=1;
                    for (unsigned int i=0;i<pointscell.size();i++){
                        mr=std::max(mr,x[(pointscell[i])+2*sd]);
                    }
                    
                    
                    for (unsigned int i=0;i<pointscell.size();i++){
                        
                        
                        Gr = x[(pointscell[i])+2*sd]/mr;
                        Gs =  ((u - x[(pointscell[i])])*(u - x[(pointscell[i])]) + (v+dim[0]-x[(pointscell[i])+sd])*(v+dim[0]-x[(pointscell[i])+sd]) );
                        WGain = 1.0/sqrt(Gs*Gr);
                        S = S + WGain;
                        Y = Y + WGain*(x[(pointscell[i])+2*sd]);
                        changed=true;
                        
                    }
                    
                    if(changed==true){
                        
                        Y=Y/S;
                        
                    }
                    
                    //HMF
                }else if(IDMETHOD==7){
                    
                    Y=0;
                    S=0;
                    
                    for (unsigned int i=0;i<pointscell.size();i++){
                        Y=Y+x[(pointscell[i])+2*sd];
                        S=S+1;
                        changed=true;
                    }
                    //                      if(changed){
                    //
                    //                         Y=Y/S;
                    //                      }
                    
                    
                    if(changed){
                        
                        
                        double mean1=Y/S;
                        S=0;
                        Y=0;
                        
                        for (unsigned int i=0;i<pointscell.size();i++) {
                            
                            //mean
                            S=S+1;
                            
                            Y=Y+(x[(pointscell[i])+2*sd]-mean1)*(x[(pointscell[i])+2*sd]-mean1);
                            
                        }
                        
                        double var=Y/S;
                        // here starts the more interesting part... the backgroundr contains a higher variance... so, while
                        // it is possible to use mean to solve the problem in most cases, situations where variance is higher than X
                        // a new solution should be employed..
                        
                        
                        
                        // just to output variance
                        Y=mean1;
                        // Y=var;
                        
                        if(var>0.40){
                            changed=false;
                            
                            //                             bool converted=false;
                            //
                            //                             for(int gi=1;gi<=(globalgrid-2);gi+=2){
                            //
                            //                                 converted=false;
                            //                                 griddisx=floor(gi/2);
                            //                                 griddisy=floor(gi/2);
                            //
                            //
                            //
                            //
                            //                                 pointscell.clear();
                            //                                 // x
                            //                                 for (int i=(u-griddisx);i<=(u+griddisx);i++){
                            //                                     // y
                            //                                     for (int f=(v-griddisy);f<=(v+griddisy);f++){
                            //
                            //                                         std::list<int> temp =matrix.get(i,f);
                            //                                         for(auto it=temp.begin();it!=temp.end();it++)
                            //                                             pointscell.push_back(*it);
                            //                                     }
                            //                                 }
                            //
                            //                                 Y=0;
                            //                                 S=0;
                            //
                            //                                 for (unsigned int i=0;i<pointscell.size();i++){
                            //                                     Y=Y+x[(pointscell[i])+2*sd];
                            //                                     S=S+1;
                            //                                     converted=true;
                            //                                 }
                            //
                            //                                 if(converted==false){
                            //                                     continue;
                            //                                 }
                            //
                            //                                 double meani=Y/S;
                            //                                 S=0;
                            //                                 Y=0;
                            //
                            //                                 for (unsigned int i=0;i<pointscell.size();i++) {
                            //                                     S=S+1;
                            //                                     Y=Y+(x[(pointscell[i])+2*sd]-meani)*(x[(pointscell[i])+2*sd]-meani);
                            //                                 }
                            //
                            //
                            //                                 double vari=Y/S;
                            //
                            //                                 if(vari<0.1){
                            //                                     Y=meani;
                            //                                     changed=true;
                            //                                     converted=true;
                            //                                     break;
                            //                                 }else{
                            //                                 converted=false;
                            //                                 }
                            //                             }
                            //
                            //                             if(!converted){
                            //                                 changed=false;
                            //                             }
                        }
                    }
                    
                    
                }else if(IDMETHOD==8){
                    
                    Y=0;
                    S=0;
                    
                    for (unsigned int i=0;i<pointscell.size();i++){
                        Y=Y+x[(pointscell[i])+2*sd];
                        S=S+1;
                        changed=true;
                    }
                    //                      if(changed){
                    //
                    //                         Y=Y/S;
                    //                      }
                    
                    
                    if(changed){
                        
                        
                        double mean1=Y/S;
                        S=0;
                        Y=0;
                        
                        for (unsigned int i=0;i<pointscell.size();i++) {
                            
                            //mean
                            S=S+1;
                            
                            Y=Y+(x[(pointscell[i])+2*sd]-mean1)*(x[(pointscell[i])+2*sd]-mean1);
                            
                        }
                        
                        double var=Y/S;
                        // here starts the more interesting part... the backgroundr contains a higher variance... so, while
                        // it is possible to use mean to solve the problem in most cases, situations where variance is higher than X
                        // a new solution should be employed..
                        
                        
                        
                        // just to output variance
                        Y=mean1;
                        // Y=var;
                        
                        if(var>0.20){
                            //create patch
                            
                            
                            // x
                            for (int ui_=(u-griddisx);ui_<=(u+griddisx);ui_++){
                                // y
                                for (int vf_=(v-griddisy);vf_<=(v+griddisy);vf_++){
                                    
                                    std::list<int> temp =matrix.get(ui_,vf_);
                                    
                                    double tf=0;
                                    int count=0;
                                    for(auto it=temp.begin();it!=temp.end();it++){
                                        tf+=x[((*it))+2*sd];
                                        count++;
                                    }
                                    if(count==0){
                                        if((ui_<0) || (vf_<0) || (ui_>= (grid*grid)) || (vf_>= (grid*grid))){
                                            
                                        }else{
                                            int indexk=(ui_-(u-griddisx))*(int)grid  + (vf_-(v-griddisy));
                                            
                                            patch[indexk]=0;
                                        }
                                    }else{
                                        int indexk=(ui_-(u-griddisx))*(int)grid  + (vf_-(v-griddisy));
                                        
                                        
                                        patch[indexk]=tf/(double)count;
                                        
                                    }
                                    
                                    //pointscell.push_back(*it);
                                    //pointscell.splice(pointscell.end(), temp);
                                    
                                }
                            }
                            
                            mxArray *lhs[1];
                            
                            mexCallMATLAB(1, lhs, 1, rhs, "Fdense");
                            
                            
                            double *pt=mxGetPr(lhs[0]);
                            
                            if(pt[60]<=0.01){
                                changed=false;
                            }else{
                                Y=pt[60];
                                changed=true;
                            }
                            // mexCallMATLAB(0, 0, 0, 0, "pause");
                            
                            
                            mxDestroyArray(lhs[0]);
                            
                            for(int i=0;i<121;i++){
                                patch[i]=0.0;
                            }
                            
                            
                        }
                    }
                    
                    
                }else if(IDMETHOD==9){
                    
                    Y=0;
                    S=0;
                    
                    for (unsigned int i=0;i<pointscell.size();i++){
                        Y=Y+x[(pointscell[i])+2*sd];
                        S=S+1;
                        changed=true;
                    }
                    //                      if(changed){
                    //
                    //                         Y=Y/S;
                    //                      }
                    
                    
                    if(changed){
                        
                        
                        double mean1=Y/S;
                        S=0;
                        Y=0;
                        
                        for (unsigned int i=0;i<pointscell.size();i++) {
                            
                            //mean
                            S=S+1;
                            
                            Y=Y+(x[(pointscell[i])+2*sd]-mean1)*(x[(pointscell[i])+2*sd]-mean1);
                            
                        }
                        
                        double var=Y/S;
                        // here starts the more interesting part... the backgroundr contains a higher variance... so, while
                        // it is possible to use mean to solve the problem in most cases, situations where variance is higher than X
                        // a new solution should be employed..
                        
                        
                        
                        // just to output variance
                        Y=mean1;
                        // Y=var;
                        
                        if(var>0.20){
                            //create patch
                            Y=-1;
                            
                            
                        }
                        
                        
                    }
                    
                    
                }else if(IDMETHOD==10){
                    
                 }else if(IDMETHOD==12){
                    //Just nearest point value (NN);
                    
                    //std::vector<Point> zlist;
                    //zlist.reserve(pointscell.size());
                    double val=999999.0;
                    Y=0;
                    
                    for (unsigned int i=0;i<pointscell.size();i++){
                        
                        
                        double local=sqrt( (u-x[(pointscell[i])])*(u-x[(pointscell[i])]) + (v+dim[0]-x[(pointscell[i])+sd])*(v+dim[0]-x[(pointscell[i])+sd])); //+x[(pointscell[i])+2*sd]*x[(pointscell[i])+2*sd]
                        if(local<val){
                            val=local;
                            Y=x[(pointscell[i])+2*sd];
                        }
                        
                    }
                    changed=true;
                    
                    
                }else if(IDMETHOD==13){
                    
                    Y=0;
                    if  (pointscell.size()==0){
                        
                    }else{
                        
                        std::vector<Point> zlist;
                        std::vector<Point> zlistcluster;
                        zlist.reserve(pointscell.size());
                        
                        for (unsigned int i=0;i<pointscell.size();i++){
                            zlist.push_back(Point(vpointscell[i],upointscell[i],x[(pointscell[i])+2*sd],pointscell[i]));
                        }
                        
                        std::sort(zlist.begin(), zlist.end(),Point::sorter);
                        
                        std::vector<std::vector<int> > clusters;
                        clusters.push_back(std::vector<int>());
                        
                        int cluster=0;
                        double jdistance=0.05;
                        double last=zlist[0].z;
                        
                        for (unsigned int i=0;i<zlist.size();i++){
                            
                            if(fabs((last-zlist[i].z)/(last+zlist[i].z))>jdistance){
                                cluster++;
                                clusters.push_back(std::vector<int>());
                                clusters[cluster].push_back(i);
                            }else{
                                clusters[cluster].push_back(i);
                            }
                            
                            last=zlist[i].z;
                        }
                        
                        if(clusters.size()==0){
                            
                        }
                        else if(clusters.size()==1){
                            zlistcluster=zlist;
                        }
                        else if(clusters.size()==2){
                            
                            double cth=(double)clusters[0].size()/(double)clusters[1].size();
                            
                            if(cth>threshold){
                                cluster=0;
                            }else{
                                cluster=1;
                            }
                            
                            for(int i=0;i<clusters[cluster].size();i++){
                                zlistcluster.push_back(zlist[clusters[cluster][i]]);
                            }
                            
                            // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                            
                            
                        }else {
                            
                            std::vector<std::pair<int,std::pair<int,double> > > ctest;
                            
                            for(int i=1;i<clusters.size();i++){
                                ctest.push_back(std::pair<int,std::pair<int,double> >(clusters[i].size(),std::pair<int,double>(i,zlist[clusters[i][0]].z)));
                            }
                            
                            std::sort(ctest.begin(), ctest.end(),clustersorter);
                            
                            auto cluster1=clusters[0];
                            auto cluster2=clusters[ctest[1].second.first];
                            
                            std::vector<std::vector<int> > clusters1;
                            clusters1.push_back(cluster1);
                            clusters1.push_back(cluster2);
                            
                            double cth=(double)clusters1[0].size()/(double)clusters1[1].size();
                            
                            if(cth>threshold){
                                cluster=0;
                            }else{
                                cluster=1;
                            }
                            
                            for(int i=0;i<clusters1[cluster].size();i++){
                                zlistcluster.push_back(zlist[clusters1[cluster][i]]);
                            }
                            
                        }
                        
                        //function q = guidedfilter(I, p, r, eps)
                        // %   GUIDEDFILTER   O(1) time implementation of guided filter.
                        // %
                        // %   - guidance image: I (should be a gray-scale/single channel image)
                        // %   - filtering input image: p (should be a gray-scale/single channel image)
                        // %   - local window radius: r
                        // %   - regularization parameter: eps
                        // [hei, wid] = size(I);
                        // N = boxfilter(ones(hei, wid), r); % the size of each local patch; N=(2r+1)^2 except for boundary pixels.
                        
                        // lest try I==p
                        // zlist=zlistcluster;
                        
                        double N=zlist.size();
                        double N_zlistcluster=zlistcluster.size();
                        //
                        // mean_I = boxfilter(I, r) ./ N;
                        double mean_I = sum(zlistcluster)/N_zlistcluster;
                        
                        // mean_p = boxfilter(p, r) ./ N;
                        double mean_p = sum(zlist)/N;
                        
                        // mean_Ip = boxfilter(I.*p, r) ./ N;
                        auto Ip=zlistcluster*zlist;
                        double mean_Ip = sum(Ip)/std::min(N,N_zlistcluster);
                        
                        // cov_Ip = mean_Ip - mean_I .* mean_p; % this is the covariance of (I, p) in each local patch.
                        
                        double cov_Ip = mean_Ip - mean_I * mean_p;
                        //
                        // mean_II = boxfilter(I.*I, r) ./ N;
                        auto II=zlistcluster*zlistcluster;
                        double mean_II=sum(II)/N_zlistcluster;
                        // var_I = mean_II - mean_I .* mean_I;
                        //
                        double var_I= mean_II - mean_I * mean_I;
                        
                        // a = cov_Ip ./ (var_I + eps); % Eqn. (5) in the paper;
                        double eps=epsstep;
                        double a=cov_Ip / (var_I + eps);
                        // b = mean_p - a .* mean_I; % Eqn. (6) in the paper;
                        
                        double b= mean_p - a * mean_I;
                        
                        
                        //
//                         // mean_a = boxfilter(a, r) ./ N;
//                         double mean_a = boxfilter(a, r) ./ N; ??
//                             // mean_b = boxfilter(b, r) ./ N;
//                             double mean_b = boxfilter(b, r) ./ N; ??
//                                 // q = mean_a .* I + mean_b; % Eqn. (8) in the paper;
//                                 Y=a * ownvalue + b;
                        // end
                        
                        
//  mexPrintf("PixelData: (%f, %f) - (%f, %f, %f) (%f, %f, %f) (%f, %f) --(%d, %d)\n",N,N_zlistcluster,mean_I,mean_p,mean_Ip,cov_Ip,mean_II,var_I,a,b,Ip.size(),II.size());
                        
                        
                        
                        
                                               Y=a*clusternear(zlistcluster,u,v)+b;

                        //N??o pode ser assim aqui!!! nem em 14!!
                               // Y=a*mean_I+b;
                        // astore.set(u,v,a);
                        // bstore.set(u,v,b);
                        
                        
                    }
                    
                    changed=true;
//
//                     https://github.com/atilimcetin/guided-filter
//                                 https://www.google-melange.com/gsoc/project/details/google/gsoc2014/vludv/5724160613416960
                }else if(IDMETHOD==14){
                    
                    Y=0;
                    if  (pointscell.size()==0){
                        
                    }else{
                        
                        std::vector<Point> zlist;
                        std::vector<Point> zlistcluster;
                        zlist.reserve(pointscell.size());
                        
                        for (unsigned int i=0;i<pointscell.size();i++){
                            zlist.push_back(Point(vpointscell[i],upointscell[i],x[(pointscell[i])+2*sd],pointscell[i]));
                        }
                        
                        //function q = guidedfilter(I, p, r, eps)
                        // %   GUIDEDFILTER   O(1) time implementation of guided filter.
                        // %
                        // %   - guidance image: I (should be a gray-scale/single channel image)
                        // %   - filtering input image: p (should be a gray-scale/single channel image)
                        // %   - local window radius: r
                        // %   - regularization parameter: eps
                        // [hei, wid] = size(I);
                        // N = boxfilter(ones(hei, wid), r); % the size of each local patch; N=(2r+1)^2 except for boundary pixels.
                        
                        // lest try I==p
                        zlistcluster=zlist;
                        
                        double N=zlist.size();
                        double N_zlistcluster=zlistcluster.size();
                        //
                        // mean_I = boxfilter(I, r) ./ N;
                        double mean_I = sum(zlistcluster)/N_zlistcluster;
                        
                        // mean_p = boxfilter(p, r) ./ N;
                        double mean_p = sum(zlist)/N;
                        
                        // mean_Ip = boxfilter(I.*p, r) ./ N;
                        auto Ip=zlistcluster*zlist;
                        double mean_Ip = sum(Ip)/std::min(N,N_zlistcluster);
                        
                        // cov_Ip = mean_Ip - mean_I .* mean_p; % this is the covariance of (I, p) in each local patch.
                        
                        double cov_Ip = mean_Ip - mean_I * mean_p;
                        //
                        // mean_II = boxfilter(I.*I, r) ./ N;
                        auto II=zlistcluster*zlistcluster;
                        double mean_II=sum(II)/N_zlistcluster;
                        // var_I = mean_II - mean_I .* mean_I;
                        //
                        double var_I= mean_II - mean_I * mean_I;
                        
                        // a = cov_Ip ./ (var_I + eps); % Eqn. (5) in the paper;
                        double eps=epsstep;//0.001;
                        double a=cov_Ip / (var_I + eps);
                        // b = mean_p - a .* mean_I; % Eqn. (6) in the paper;
                        
                        double b= mean_p - a * mean_I;
                        //
//                         // mean_a = boxfilter(a, r) ./ N;
//                         double mean_a = boxfilter(a, r) ./ N; ??
//                             // mean_b = boxfilter(b, r) ./ N;
//                             double mean_b = boxfilter(b, r) ./ N; ??
//
//
//                                 //
//                                 // q = mean_a .* I + mean_b; % Eqn. (8) in the paper;
//                                 Y=a * ownvalue + b;
                        // end
//   mexPrintf("PixelData: (%f, %f) - (%f, %f, %f) (%f, %f, %f) (%f, %f) --(%d, %d)\n",N,N_zlistcluster,mean_I,mean_p,mean_Ip,cov_Ip,mean_II,var_I,a,b,Ip.size(),II.size());
                        
                        
                       Y=a*clusternear(zlistcluster,u,v)+b;
                       // Y=a*mean_I+b;
                        
                        //  Y=mean_I;
                        //  astore.set(u,v,a);
                        //  bstore.set(u,v,b);
                    }
                    
                    changed=true;
//
//                     https://github.com/atilimcetin/guided-filter
//                                 https://www.google-melange.com/gsoc/project/details/google/gsoc2014/vludv/5724160613416960
                }else if(IDMETHOD==17){
                    
                    
                    if  (pointscell.size()==0){
                        
                    }else{
                        
                        //                        Y=0;
                        //                        changed=true;
                        
                        
                        std::vector<Point> zlist;
                        zlist.reserve(pointscell.size());
                        for (unsigned int i=0;i<pointscell.size();i++){
                            zlist.push_back(Point(x[(pointscell[i])],x[(pointscell[i])+sd],x[(pointscell[i])+2*sd],pointscell[i]));
                        }
                        
                        double nearz=clusterzmin(zlist);
                        
                        
                        
                          for(int i=0;i<zlist.size();i++){
                                
                                Gs =  ((u - zlist[i].x)*(u -zlist[i].x) + (v+dim[0]-zlist[i].y)*(v+dim[0]-zlist[i].y) );
                              //  WGain = 1.0/(0.01+sqrt(Gs));
                                WGain = (1.0/(1.0+sqrt(Gs)))*(1.0/ (1+fabs(nearz-zlist[i].z)));
                                
                               // WGain =(1.0/sqrt(2.0*M_PI))*exp(-Gs/10.0)*(1.0/ zlist[clusters[cluster][i]].z);

                                S = S + WGain;
                                Y = Y + WGain*(zlist[i].z);
                                changed=true;
                            }
                            
                            
                            Y=Y/S;
                        
                    }
                    
                }
                
                
                
                else if(IDMETHOD==20){
                
 
                    
                    
                    if  (pointscell.size()==0){
                        
                    }else{
                        
                        //                        Y=0;
                        //                        changed=true;
                        
                        
                        std::vector<Point> zlist;
                        zlist.reserve(pointscell.size());
                        for (unsigned int i=0;i<pointscell.size();i++){
                            zlist.push_back(Point(x[(pointscell[i])],x[(pointscell[i])+sd],x[(pointscell[i])+2*sd],pointscell[i]));
                        }
                        
                        double nearz=clusterzmin(zlist);

                        //sort
                        std::sort(zlist.begin(), zlist.end(),Point::sorter);
                        
                        
                        //limits
                        //                         auto maxvalue=std::max_element(zlist.begin(), zlist.end());
                        //                         int maxid=std::distance(zlist.begin(), maxvalue);
                        //                         auto minvalue=std::min_element(zlist.begin(), zlist.end());;
                        //                         int minid=std::distance(zlist.begin(), minvalue);
                        
                        
                        
                        std::vector<std::vector<int> > clusters;
                        clusters.push_back(std::vector<int>());
                        
                        
                        
                        //cluster
                        //http://home.isr.uc.pt/~cpremebida/files_cp/Segmentation%20and%20Geometric%20Primitives%20Extraction%20from%202D%20Laser%20Range%20Data%20for%20Mobile%20Robot%20Applications.pdf
                        int cluster=0;
                        double jdistance=epsstep;
                        double last=zlist[0].z;
                        
                        for (unsigned int i=0;i<zlist.size();i++){
                            
                            if(fabs((last-zlist[i].z)/(last+zlist[i].z))>jdistance){
                                cluster++;
                                clusters.push_back(std::vector<int>());
                                clusters[cluster].push_back(i);
                            }else{
                                clusters[cluster].push_back(i);
                            }
                            
                            last=zlist[i].z;
                        }
                        
                        
                        
                        S=0;
                    
                        
                        
                             changed=true;
                        
                             Y=(clusters.size()>1)?30:0;
                             
                             }
                
                     }

                else if(IDMETHOD==11){
                    
                    
                    if  (pointscell.size()==0){
                        
                    }else{
                        
                        //                        Y=0;
                        //                        changed=true;
                        
                        
                        std::vector<Point> zlist;
                        zlist.reserve(pointscell.size());
                        for (unsigned int i=0;i<pointscell.size();i++){
                            zlist.push_back(Point(x[(pointscell[i])],x[(pointscell[i])+sd],x[(pointscell[i])+2*sd],pointscell[i]));
                        }
                        
                        double nearz=clusterzmin(zlist);

                        //sort
                        std::sort(zlist.begin(), zlist.end(),Point::sorter);
                        
                        
                        //limits
                        //                         auto maxvalue=std::max_element(zlist.begin(), zlist.end());
                        //                         int maxid=std::distance(zlist.begin(), maxvalue);
                        //                         auto minvalue=std::min_element(zlist.begin(), zlist.end());;
                        //                         int minid=std::distance(zlist.begin(), minvalue);
                        
                        
                        
                        std::vector<std::vector<int> > clusters;
                        clusters.push_back(std::vector<int>());
                        
                        
                        
                        //cluster
                        //http://home.isr.uc.pt/~cpremebida/files_cp/Segmentation%20and%20Geometric%20Primitives%20Extraction%20from%202D%20Laser%20Range%20Data%20for%20Mobile%20Robot%20Applications.pdf
                        int cluster=0;
                        double jdistance=epsstep;
                        double last=zlist[0].z;
                        
                        for (unsigned int i=0;i<zlist.size();i++){
                            
                            if(fabs((last-zlist[i].z)/(last+zlist[i].z))>jdistance){
                                cluster++;
                                clusters.push_back(std::vector<int>());
                                clusters[cluster].push_back(i);
                            }else{
                                clusters[cluster].push_back(i);
                            }
                            
                            last=zlist[i].z;
                        }
                        
                        
                        
                        S=0;
                        
                        
                        //clean up
                        
//
//                         for(int i=0;i<clusters.size();i++){
//
//                             if(clusters[i].size()<=1){
//
//                                 clusters.erase(clusters.begin()+i);
//                                 i--;
//
//                             }
//
//
//                         }
                        
                        
                        
                        
//                         if(clusters.size()>1){
//
//                             Y=9;
//                         }else{
//                             Y=1;
                        Y=0;
//                         }
                        
                        if(clusters.size()==0){
                            
                            Y=0;
                            
                        }
                        else if(clusters.size()==1){
                            
                            
                            for(int i=0;i<clusters[cluster].size();i++){
                                
                                Gs =  ((u - zlist[clusters[cluster][i]].x)*(u - zlist[clusters[cluster][i]].x) + (v+dim[0]-zlist[clusters[cluster][i]].y)*(v+dim[0]-zlist[clusters[cluster][i]].y) );
                              //  WGain = 1.0/(0.01+sqrt(Gs));
                                WGain = (1.0/(1.0+sqrt(Gs)))*(1.0/ (1+fabs(nearz-zlist[clusters[cluster][i]].z)));
                                
                               // WGain =(1.0/sqrt(2.0*M_PI))*exp(-Gs/10.0)*(1.0/ zlist[clusters[cluster][i]].z);

                                S = S + WGain;
                                Y = Y + WGain*(zlist[clusters[cluster][i]].z);
                                changed=true;
                            }
                            
                            
                            Y=Y/S;
//                                Y=33;
                        }
                        else if(clusters.size()==2){
                            
                            double cth=((double)clusters[0].size())/((double)clusters[1].size());
                            
                            if(cth>threshold){
                                cluster=0;
                            }else{
                                cluster=1;
                            }
 
                            for(int i=0;i<clusters[cluster].size();i++){
                                
//                                 mexPrintf("points (%f, %f) - (%f, %f) \n",(double)u,(double)v,zlist[clusters[cluster][i]].x,-(double)dim[0]+zlist[clusters[cluster][i]].y);
                                
                                Gs =  (((double)u - zlist[clusters[cluster][i]].x)*((double)u - zlist[clusters[cluster][i]].x) + ((double)v+(double)dim[0]-zlist[clusters[cluster][i]].y)*((double)v+(double)dim[0]-zlist[clusters[cluster][i]].y) );
                                WGain = (1.0/(1.0+sqrt(Gs)))*(1.0/ (1+fabs(nearz-zlist[clusters[cluster][i]].z)));
                               // WGain =(1.0/sqrt(2.0*M_PI))*exp(-Gs/10.0)*(1.0/ zlist[clusters[cluster][i]].z);

                                S = S + WGain;
                                Y = Y + WGain*(zlist[clusters[cluster][i]].z);
                                changed=true;
                            }
                            
                            
                            Y=Y/S;
//
//                             if(cluster==0){
//                                 Y=33;
//                             }else{
//                                 Y=66;
//                             }
                            
                        }else {
                            Y=0;
//                              Y=100;
//
                            std::vector<std::pair<int,std::pair<int,double> > > ctest;
                            
                            for(int i=1;i<clusters.size();i++){
                                ctest.push_back(std::pair<int,std::pair<int,double> >(clusters[i].size(),std::pair<int,double>(i,zlist[clusters[i][0]].z)));
                            }
                            
                            
                            std::sort(ctest.begin(), ctest.end(),clustersorter);
                            
                            auto cluster1=clusters[0];
                            auto cluster2=clusters[ctest[1].second.first];
                            
                            std::vector<std::vector<int> > clusters1;
                            //
//                             if(zlist[cluster1[0]].z<zlist[cluster2[0]].z){
                            
                            clusters1.push_back(cluster1);
                            clusters1.push_back(cluster2);
                            
                            
//                             }else{
                            
//                                 clusters.push_back(cluster2);
//                                 clusters.push_back(cluster1);
//
//                             }
                            
                            
                           // double cth=(double)clusters1[0].size()/(double)clusters1[1].size();
                           double cth=((double)clusters1[0].size())/((double)clusters1[1].size());

                            
                            if(cth>threshold){
                                cluster=0;
                            }else{
                                cluster=1;
                            }
                            
                            
                            
                            for(int i=0;i<clusters1[cluster].size();i++){
                                Gs =  ((u - zlist[clusters1[cluster][i]].x)*(u - zlist[clusters1[cluster][i]].x) + (v+dim[0]-zlist[clusters1[cluster][i]].y)*(v+dim[0]-zlist[clusters1[cluster][i]].y) );
                                //WGain = 1.0/sqrt(Gs);
                               // WGain = (1.0/(0.01+sqrt(Gs)))*(1.0/ (0.01+fabs(nearz-zlist[clusters1[cluster][i]].z)));
                               WGain = (1.0/(1.0+sqrt(Gs)))*(1.0/ (1+fabs(nearz-zlist[clusters1[cluster][i]].z)));

                               //   WGain =(1.0/sqrt(2.0*M_PI))*exp(-Gs/10.0)*(1.0/ zlist[clusters1[cluster][i]].z);

                                S = S + WGain;
                                Y = Y + WGain*(zlist[clusters1[cluster][i]].z);
                                changed=true;
                            }
                            
                            
                            Y=Y/S;
                            
                            
                            
                        }
                        
                        
                        
                        
                        
                        
                        
                        
                        changed=true;
                        
                        
                        
//
//
//                           if(clusters.size()==2){
//
//
//                               mexPrintf("clusters %d : %d - [",clusters.size(),zlist.size());
//
//                               for(int i=0;i<zlist.size();i++){
//                                mexPrintf(" %f",zlist[i].z);
//                               }
//                               mexPrintf(" ]  -- ");
//
//
//                               for(int i=0;i<clusters.size();i++){
//                                   mexPrintf(" ( ");
//                                   for(int g=0;g<clusters[i].size();g++){
//                                mexPrintf(" %d",clusters[i][g]);
//                               }
//                                                                   mexPrintf(" ) ");
//
//                                   }
//                                                               mexPrintf(" \n ");
//
//
//
//
//
//                           }
//
//
                        
                        
                        
//                        //decision
                        
//                        cluster=0;
//                        if(clusters.size()==1){
                        
                        
                        
                        
//                        }else if (clusters.size()>1){
                        
                        
                        
                        
                        
//                        }else{
//                            //size == 0
                        
                        
                        
//                        }
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        //                        for (unsigned int i=0;i<pointscell.size();i++){
                        
                        
                        //                            // Gr = x[(pointscell[i])+2*sd]/mr;
                        
                        
                        
                        
                        //                            Gs =  ((u - x[(pointscell[i])])*(u - x[(pointscell[i])]) + (v+dim[0]-x[(pointscell[i])+sd])*(v+dim[0]-x[(pointscell[i])+sd]) );
                        //                            WGain = 1.0/sqrt(Gs);
                        //                            S = S + WGain;
                        //                            Y = Y + WGain*(x[(pointscell[i])+2*sd]);
                        //                            changed=true;
                        
                        //                        }
                        
//                         if(changed==true){
//
//                             Y=Y/S;
//
//                         }
                        
                        
                    }
                    
                }
                
                
                
                
                if(!changed){
                    
                    y[u*(int)dim[1]  + v]=0;
                    
                }
                
                if(changed){
                    
                    y[u*(int)dim[1]  + v] = Y;
                    grid=globalgrid;
                    break;
                    
                }
                
                grid=grid+2;
                
            }
            
            
            
        }
        
        
    }
    
    /*
     * if(IDMETHOD==13||IDMETHOD==14){
     *
     *
     *
     * //http://research.microsoft.com/en-us/um/people/kahe/publications/pami12guidedfilter.pdf
     * // column
     * for (int u=0; u<(int)dim[2]; u=u+1)
     * {
     * // line
     * for (int v=0; v<(int)dim[1]; v=v+1)
     * {
     * //                double ownvalue=0;
     * //
     * //                std::list<int> temp =matrix.getUnique(u,v,x,2*sd);
     * //
     * //                for(auto it=temp.begin();it!=temp.end();it++)
     * //                    ownvalue=x[(*it)+2*sd];
     *
     *
     * double a_line=astore.get(u,v);
     * double b_line=bstore.get(u,v);
     *
     * int ct=1;
     *
     *
     * Y=(a_line/((double)ct))*y[u*(int)dim[1]  + v]+(b_line/((double)ct));
     * y[u*(int)dim[1]  + v] = Y;
     *
     *
     * //                 for (int i=(u-griddisx);i<=(u+griddisx);i++){
     * //                     // y
     * //                     for (int f=(v-griddisy);f<=(v+griddisy);f++){
     * //
     * //                         // std::list<int> temp =matrix.getUnique(i,f,x,2*sd);
     * //                         double temp =astore.get(i,f);
     * //
     * //                         if(temp!=0.0){
     * //                             ct++;
     * //                             a_line+=temp;
     * //                             temp =bstore.get(i,f);
     * //                             b_line+=temp;
     * //                         }
     * //                     }
     * //                  }
     * //                 if(ct!=0){
     * //                     Y=(a_line/((double)ct))*y[u*(int)dim[1]  + v]+(b_line/((double)ct));
     * //                     y[u*(int)dim[1]  + v] = Y;
     * //                 }
     * }
     * }
     *
     *
     * }
     */
} /*End "calc_Dense" Function*/
/****************************/

void mexFunction(
        int          nlhs,
        mxArray      *plhs[],
        int          nrhs,
        const mxArray *prhs[]
        )
{
    
    double  *Lidar;
    double  *par;
    const mwSize  *dims;
    double *Y;
    
    /* Check for proper number of arguments */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin",
                "MEXCPP requires 2 input arguments.");
    } else if (nlhs >= 2) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
                "MEXCPP requires 1 output argument.");
    }
    
    Lidar   = mxGetPr(prhs[0]);
    par     = mxGetPr(prhs[1]);
    dims    = mxGetDimensions(prhs[0]);
    
    plhs[0]= mxCreateDoubleMatrix(par[1],par[2],mxREAL);
    Y =   mxGetPr(plhs[0]);
    
   // mexPrintf("dims %d %d---- %d , %d\n",sizeof(size_t),mxIsDouble(prhs[1]),dims[0],dims[1]);
    mexPrintf("(%d) \tINPUT: %d  \t%d  \t%d  \t%d  \t%d  \t%f  \t%f \n",(int)par[7],(int)par[0],(int)par[1],(int)par[2],(int)par[3],(int)par[4],par[5],par[6]);
    
    
   // mexEvalString("disp('--------')");
    calc_Dense(Lidar,Y,par,(int)dims[0],(int)par[3],(int)par[4],(double)par[5],(double)par[6]); // <= method goes here!!!!!!!!!!!!!!!!!!!!!!!!
  //  mexEvalString("disp('--------')");
    
}
