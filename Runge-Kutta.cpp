#include <iostream>
#include <math.h>
#include <fstream>

double tol=pow(10,-8);
double x=0, y=0;
double epsilon=pow(10,-7);
double delta=pow(10,-8);
int step=0;

using namespace std;

double fx(double t, double x, double y)
{
    double res;
    res=y+0*t+0*x;
    return res;
}

double fy(double t, double x, double y)
{
    double res;
    res=0.1*(1-x*x)*y-x+0*t;
    return res;
}

void graph(double T,double tol, double x0, double y0, ofstream &fout)
{
    //int n=0;
    double /*x0=0, y0=1.0,*/ x1=0, x2=0, y1=0, y2=0, count=0.0, h=0.1, errx=0, erry=0, err=0, facmax=1.3, facmin=0.7, fac=0.9;
    double kx1=0, kx2=0, kx3=0, kx4=0, kx5=0, kx6=0, kx7=0;
    double ky1=0, ky2=0, ky3=0, ky4=0, ky5=0, ky6=0, ky7=0;
    double c2=1.0/5, c3=3.0/10, c4=4.0/5, c5=8.0/9, c6=1.0, c7=1.0;
    double a21=1.0/5;
    double a31=3.0/40, a32=9.0/40;
    double a41=44.0/45, a42=-56.0/15, a43=32.0/9;
    double a51=19372.0/6561, a52=-25360.0/2187, a53=64448.0/6561, a54=-212.0/729;
    double a61=9017.0/3168, a62=-355.0/33, a63=46732.0/5247, a64=49.0/176, a65=-5103.0/18656;
    double a71=35.0/384, a72=0.0, a73=500.0/1113, a74=125.0/192, a75=-2187.0/6784, a76=11.0/84;
    double b1=35.0/384, b2=0.0, b3=500.0/1113, b4=125.0/192, b5=-2187.0/6784, b6=11.0/84, b7=0.0;
    double p1=5179.0/57600, p2=0.0, p3=7571.0/16695, p4=393.0/640, p5=-92097.0/339200, p6=187.0/2100, p7=1.0/40;

    while(count<T)
    {
        if(count+h>T)
        {
            h=T-count;
        }
        kx1=fx(count,x0,y0);
        ky1=fy(count,x0,y0);

        kx2=fx(count+c2*h, x0+h*a21*kx1, y0+h*a21*ky1);
        ky2=fy(count+c2*h, x0+h*a21*kx1, y0+h*a21*ky1);

        kx3=fx(count+c3*h, x0+h*(a31*kx1+a32*kx2), y0+h*(a31*ky1+a32*ky2));
        ky3=fy(count+c3*h, x0+h*(a31*kx1+a32*kx2), y0+h*(a31*ky1+a32*ky2));

        kx4=fx(count+c4*h, x0+h*(a41*kx1+a42*kx2+a43*kx3), y0+h*(a41*ky1+a42*ky2+a43*ky3));
        ky4=fy(count+c4*h, x0+h*(a41*kx1+a42*kx2+a43*kx3), y0+h*(a41*ky1+a42*ky2+a43*ky3));

        kx5=fx(count+c5*h, x0+h*(a51*kx1+a52*kx2+a53*kx3+a54*kx4), y0+h*(a51*ky1+a52*ky2+a53*ky3+a54*ky4));
        ky5=fy(count+c5*h, x0+h*(a51*kx1+a52*kx2+a53*kx3+a54*kx4), y0+h*(a51*ky1+a52*ky2+a53*ky3+a54*ky4));

        kx6=fx(count+c6*h, x0+h*(a61*kx1+a62*kx2+a63*kx3+a64*kx4+a65*kx5), y0+h*(a61*ky1+a62*ky2+a63*ky3+a64*ky4+a65*ky5));
        ky6=fy(count+c6*h, x0+h*(a61*kx1+a62*kx2+a63*kx3+a64*kx4+a65*kx5), y0+h*(a61*ky1+a62*ky2+a63*ky3+a64*ky4+a65*ky5));

        kx7=fx(count+c7*h, x0+h*(a71*kx1+a72*kx2+a73*kx3+a74*kx4+a75*kx5+a76*kx6), y0+h*(a71*ky1+a72*ky2+a73*ky3+a74*ky4+a75*ky5+a76*ky6));
        ky7=fy(count+c7*h, x0+h*(a71*kx1+a72*kx2+a73*kx3+a74*kx4+a75*kx5+a76*kx6), y0+h*(a71*ky1+a72*ky2+a73*ky3+a74*ky4+a75*ky5+a76*ky6));

        x1=x0+h*(b1*kx1+b2*kx2+b3*kx3+b4*kx4+b5*kx5+b6*kx6+b7*kx7);
        x2=x0+h*(p1*kx1+p2*kx2+p3*kx3+p4*kx4+p5*kx5+p6*kx6+p7*kx7);
        errx=fabs(x2-x1);

        y1=y0+h*(b1*ky1+b2*ky2+b3*ky3+b4*ky4+b5*ky5+b6*ky6+b7*ky7);
        y2=y0+h*(p1*ky1+p2*ky2+p3*ky3+p4*ky4+p5*ky5+p6*ky6+p7*ky7);
        erry=fabs(y2-y1);

        err=sqrt((pow(errx,2)+pow(erry,2))/2);
        //err=fabs(errx-erry)/31;

        if (err<=tol)
        {
            if (step==13)
            {
                fout <<x0<<" "<< y0<< "\n";
            }

            x0=x1;
            y0=y1;

            count+=h;
            //n++;
            if(T-count<epsilon)
            {
                x=x0;
                y=y0;
                break;
            }
        }
        h=h*fmin(facmax, fmax(facmin, fac*pow(tol/err,1.0/6)));

    }
    //x=x0;
    //y=y0;
    //fout << x0 << " " << y0 << "\n";
    //return n;
}

int main()
{
    int N=0;
    double a=0.0, b=0.0, x0=0, y0=0, x11=0, x12=0, y11=0, y12=0;
    double alpha1=1.5, alpha2=50.0;
    double alpha10=0, alpha20=0;
    double X1=0, X2=0;
    double X11=0, X12=0, X21=0, X22=0;
    double h=1.0;
    ofstream data("data.txt");

    graph(alpha2,tol,alpha1,b,data);
    X1=x-alpha1;
    X2=y-b;
    cout << X1 << " " << X2 << endl;
    //while (fmax(fabs(X1),fabs(X2))>epsilon)
    while (sqrt(pow(X1,2)+pow(X2,2))>epsilon)
    {
        step++;
        cout<<"step="<<step<<"\n";

        x0=X1;
        y0=X2;
        graph(alpha2,tol,alpha1+delta,b,data);
        x11=x-alpha1-delta;
        y11=y-b;
        X11=(x11-x0)/delta;
        X21=(y11-y0)/delta;
        graph(alpha2+delta,tol,alpha1,b,data);
        cout <<x<<endl;
        x12=x-alpha1;
        cout <<x12<<endl;
        y12=y-b;
        cout <<y12<<endl;
        X12=(x12-x0)/delta;
        X22=(y12-y0)/delta;
        cout <<X11 << " " << X12 << " " <<X21 << " " <<X22 <<endl;
        alpha1=alpha1+h*(-x0*X22+y0*X12)/(X11*X22-X12*X21);//По методу Крамера
        alpha10=alpha1;
        alpha2=alpha2+h*(-y0*X11+x0*X21)/(X11*X22-X12*X21);
        alpha20=alpha2;
        cout <<alpha1 << " " << alpha2 <<endl;
        graph(alpha2,tol,alpha1,b,data);
        X1=x-alpha1;
        X2=y-b;
        N++;
        while (sqrt(pow(X1,2)+pow(X2,2)) > (sqrt(pow(x0,2)+pow(y0,2))) and h>pow(2,-30))
        {

            h=h/2.0;
            alpha1=alpha10+h*(-x0*X22+y0*X12)/(X11*X22-X12*X21);//По методу Крамера
            alpha2=alpha20+h*(-y0*X11+x0*X21)/(X11*X22-X12*X21);
            //N++;
            graph(alpha2,tol,alpha1,b,data);
            X1=x-alpha1;
            X2=y-b;

        }
        h=1.0;

    }
    cout <<"Number of iterations="<<N<<"\n";
    return 0;
}
