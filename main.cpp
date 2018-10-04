#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

#define eps 0.0001


double function1(double x, double y)
{
    return tan(x*y)-pow(x,2);
}

double function2(double x, double y)
{
    return 0.8*pow(x,2)+2*pow(y,2)-1;
}

double func11(double x, double y)
{
   // return  -y/x + 2*pow(cos(x*y),2);
    return (y*pow((1/cos(x*y)),2))-(2*x);
}

double func12(double x, double y)
{
   // return x/(x-y+x*cos(2*x*y));
    return x*pow((1/cos(x*y)),2);
}

double func21(double x, double y)
{
    //return - (2*x)/(5*y);
    return 1.6*x;
}

double func22(double x, double y)
{
   // return -(5*x)/(2*y);
    return 4*y;
}

double Jac(double x,double y)
{
    return func11(x,y)*func22(x,y) - func12(x,y)*func21(x,y);

}
float* del(float *a, float M, int N){
    for (int i = 0; i < N; i++)
        a[i] /= M;
    return a;
 }
double* Jordan_Gaus(double mass[2][2], double dx, double dy,double x[2]){
    float *buf1, buf2, el;
    int r=2;
    float **copymass = new float*[r+1];

    for(int j = 0; j < r; j++)
        copymass[j] = new float[r];
        for(int i = 0;i < r; i++ ){
        for(int j = 0; j < r; j++ )
            copymass[i][j] = mass[i][j];
        }

       for(int j = 1; j < r; j++)
       if(copymass[0][0]<copymass[j][0]){
        buf1=copymass[j];
        buf2=x[j];
        copymass[j]=copymass[0];
        x[j]=x[0];
        copymass[0]=buf1;
        x[0]=buf2;
       }

    for(int i = 0; i < r; i++ ){
            x[i] /= copymass[i][i];
            copymass[i] = del(copymass[i], copymass[i][i], r);
            if (i < r-1){


                for(int k = i+1; k < r; k++){
                float el = copymass[k][i];
                x[k]-= el * x[i];
                    for(int m = 0; m < r; m++){
                        copymass[k][m] -= el * copymass[i][m];

                    }
                }
            }
    }
        for(int i = r-1; i > -1; i--){
            for(int k = i-1 ; k > -1; k--){
                float el = copymass[k][i];
                    x[k]-= el * x[i];
                    for(int m = r-1; m > -1; m--){
                        copymass[k][m] -= el * copymass[i][m];

                    }
            }
        }
        dx=x[0];
        dy=x[1];

     double *re=new double[2];
     re[0]=dx;
     re[1]=dy;
     return re;

}

int main()
{
    double x, y,l=0,k=1;
    cout << "x = ";
    cin >> x ;
    cout << "y = ";
    cin >> y;
    double a[2][2], b[2];
    double dx, dy, *dxy;
    do{

     a[0][0]=func11(x,y);
     a[0][1]=func12(x,y);
     a[1][0]=func21(x,y);
     a[1][1]=func22(x,y);
     b[0]=-function1(x,y);
     b[1]=-function2(x,y);
     dxy=Jordan_Gaus(a,dx,dy,b);
     dx=dxy[0];
     dy=dxy[1];
    x+=dx;
    y+=dy;
    k++;

    }while((abs(dx)>eps)&&(abs(dy)>eps));
    cout << fixed << setprecision(3)<<"x = "<<x<<" y = "<<y<<endl;
    cout << fixed << setprecision(3)<<"F1 = "<<function1(x,y)<<endl<<"F2 = "<<function2(x,y);
    return 0;

}
