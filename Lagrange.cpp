#include <iostream>
#include <math.h>
#include <fstream>

double tol=pow(10,-8);
double x=0, y=0, px=0, py=0, B=0;
double per=0;
double epsilon=pow(10,-7);
double delta=pow(10,-8);
int step=0;
double x_next = 0;
double tmp=0;
double x_prev =0;
double x_curr=0;
int num=0;
int coco=0;
int m=0;
//int papa=0;

using namespace std;

//double method_chord(double x_prev, double x_curr, p1, p2, double e);

double fx(double t, double x, double y, double px, double py)
{
    double res;
    res=y+0*t+0*x+0*px+0*py;
    return res;
}

double fy(double t, double x, double y, double px, double py)
{
    double res;
    if (py>1.0)
    {
        res=1.0+0*t+0*x+0*y+0*px+0*py;
    }
    else if (py<-1.0)
    {
        res=-1.0+0*t+0*x+0*y+0*px+0*py;
    }
    else
    {
        res=py+0*t+0*x+0*y+0*px+0*py;
    }
    return res;
}

double fpx(double t, double x, double y, double px, double py)
{
    double res;
    res=-x+0*t+0*y+0*px+0*py;
    return res;
}

double fpy(double t, double x, double y, double px, double py)
{
    double res;
    res=-px+0*t+0*x+0*y+0*py;
    return res;
}

double fB(double t, double x, double y, double px, double py)
{
    double res;
    res=pow(fy(t,x,y,px,py),2)-x*x;
    return res;
}

void graph(double count, double T,double tol, double x0, double y0, double px0, double py0, double B0, ofstream &fout)
{
    int n=0;
    int papa=0;
    double /*x0=0, y0=1.0,*/ x1=0, x2=0, y1=0, y2=0, px1=0, px2=0, py1=0, py2=0, /*count=0.0,*/ h=0.0001, facmax=1.3, facmin=0.7, fac=0.9, t_next=0, t_curr=-1;
    double /*B0=0,*/ B1=0, B2=0;
    double x10=0, y10=0, px10=0, py10=0;
    double errx=0, erry=0, errpx=0, errpy=0, errB=0, err=0;
    double var1=0, var2=0, var3=0, var4=0;
    double kx1=0, kx2=0, kx3=0, kx4=0, kx5=0, kx6=0, kx7=0;
    double ky1=0, ky2=0, ky3=0, ky4=0, ky5=0, ky6=0, ky7=0;
    double kpx1=0, kpx2=0, kpx3=0, kpx4=0, kpx5=0, kpx6=0, kpx7=0;
    double kpy1=0, kpy2=0, kpy3=0, kpy4=0, kpy5=0, kpy6=0, kpy7=0;
    double kB1=0, kB2=0, kB3=0, kB4=0, kB5=0, kB6=0, kB7=0;
    double c2=1.0/5, c3=3.0/10, c4=4.0/5, c5=8.0/9, c6=1.0, c7=1.0;
    double a21=1.0/5;
    double a31=3.0/40, a32=9.0/40;
    double a41=44.0/45, a42=-56.0/15, a43=32.0/9;
    double a51=19372.0/6561, a52=-25360.0/2187, a53=64448.0/6561, a54=-212.0/729;
    double a61=9017.0/3168, a62=-355.0/33, a63=46732.0/5247, a64=49.0/176, a65=-5103.0/18656;
    double a71=35.0/384, a72=0.0, a73=500.0/1113, a74=125.0/192, a75=-2187.0/6784, a76=11.0/84;
    double b1=35.0/384, b2=0.0, b3=500.0/1113, b4=125.0/192, b5=-2187.0/6784, b6=11.0/84, b7=0.0;
    double p1=5179.0/57600, p2=0.0, p3=7571.0/16695, p4=393.0/640, p5=-92097.0/339200, p6=187.0/2100, p7=1.0/40;
    ofstream data1("data1.txt");

    while(count<T)
    {
        if(count+h>T)
        {
            h=T-count;
        }
        kx1=fx(count,x0,y0,px0,py0);
        ky1=fy(count,x0,y0,px0,py0);
        kpx1=fpx(count,x0,y0,px0,py0);
        kpy1=fpy(count,x0,y0,px0,py0);
        kB1=fB(count,x0,y0,px0,py0);

        kx2=fx(count+c2*h, x0+h*a21*kx1, y0+h*a21*ky1, px0+h*a21*kpx1, py0+h*a21*kpy1);
        ky2=fy(count+c2*h, x0+h*a21*kx1, y0+h*a21*ky1, px0+h*a21*kpx1, py0+h*a21*kpy1);
        kpx2=fpx(count+c2*h, x0+h*a21*kx1, y0+h*a21*ky1, px0+h*a21*kpx1, py0+h*a21*kpy1);
        kpy2=fpy(count+c2*h, x0+h*a21*kx1, y0+h*a21*ky1, px0+h*a21*kpx1, py0+h*a21*kpy1);
        kB2=fB(count+c2*h, x0+h*a21*kx1, y0+h*a21*ky1, px0+h*a21*kpx1, py0+h*a21*kpy1);

        kx3=fx(count+c3*h, x0+h*(a31*kx1+a32*kx2), y0+h*(a31*ky1+a32*ky2), px0+h*(a31*kpx1+a32*kpx2), py0+h*(a31*kpy1+a32*kpy2));
        ky3=fy(count+c3*h, x0+h*(a31*kx1+a32*kx2), y0+h*(a31*ky1+a32*ky2), px0+h*(a31*kpx1+a32*kpx2), py0+h*(a31*kpy1+a32*kpy2));
        kpx3=fpx(count+c3*h, x0+h*(a31*kx1+a32*kx2), y0+h*(a31*ky1+a32*ky2), px0+h*(a31*kpx1+a32*kpx2), py0+h*(a31*kpy1+a32*kpy2));
        kpy3=fpy(count+c3*h, x0+h*(a31*kx1+a32*kx2), y0+h*(a31*ky1+a32*ky2), px0+h*(a31*kpx1+a32*kpx2), py0+h*(a31*kpy1+a32*kpy2));
        kB3=fB(count+c3*h, x0+h*(a31*kx1+a32*kx2), y0+h*(a31*ky1+a32*ky2), px0+h*(a31*kpx1+a32*kpx2), py0+h*(a31*kpy1+a32*kpy2));

        kx4=fx(count+c4*h, x0+h*(a41*kx1+a42*kx2+a43*kx3), y0+h*(a41*ky1+a42*ky2+a43*ky3), px0+h*(a41*kpx1+a42*kpx2+a43*kpx3), py0+h*(a41*kpy1+a42*kpy2+a43*kpy3));
        ky4=fy(count+c4*h, x0+h*(a41*kx1+a42*kx2+a43*kx3), y0+h*(a41*ky1+a42*ky2+a43*ky3), px0+h*(a41*kpx1+a42*kpx2+a43*kpx3), py0+h*(a41*kpy1+a42*kpy2+a43*kpy3));
        kpx4=fpx(count+c4*h, x0+h*(a41*kx1+a42*kx2+a43*kx3), y0+h*(a41*ky1+a42*ky2+a43*ky3), px0+h*(a41*kpx1+a42*kpx2+a43*kpx3), py0+h*(a41*kpy1+a42*kpy2+a43*kpy3));
        kpy4=fpy(count+c4*h, x0+h*(a41*kx1+a42*kx2+a43*kx3), y0+h*(a41*ky1+a42*ky2+a43*ky3), px0+h*(a41*kpx1+a42*kpx2+a43*kpx3), py0+h*(a41*kpy1+a42*kpy2+a43*kpy3));
        kB4=fB(count+c4*h, x0+h*(a41*kx1+a42*kx2+a43*kx3), y0+h*(a41*ky1+a42*ky2+a43*ky3), px0+h*(a41*kpx1+a42*kpx2+a43*kpx3), py0+h*(a41*kpy1+a42*kpy2+a43*kpy3));

        kx5=fx(count+c5*h, x0+h*(a51*kx1+a52*kx2+a53*kx3+a54*kx4), y0+h*(a51*ky1+a52*ky2+a53*ky3+a54*ky4), px0+h*(a51*kpx1+a52*kpx2+a53*kpx3+a54*kpx4), py0+h*(a51*kpy1+a52*kpy2+a53*kpy3+a54*kpy4));
        ky5=fy(count+c5*h, x0+h*(a51*kx1+a52*kx2+a53*kx3+a54*kx4), y0+h*(a51*ky1+a52*ky2+a53*ky3+a54*ky4), px0+h*(a51*kpx1+a52*kpx2+a53*kpx3+a54*kpx4), py0+h*(a51*kpy1+a52*kpy2+a53*kpy3+a54*kpy4));
        kpx5=fpx(count+c5*h, x0+h*(a51*kx1+a52*kx2+a53*kx3+a54*kx4), y0+h*(a51*ky1+a52*ky2+a53*ky3+a54*ky4), px0+h*(a51*kpx1+a52*kpx2+a53*kpx3+a54*kpx4), py0+h*(a51*kpy1+a52*kpy2+a53*kpy3+a54*kpy4));
        kpy5=fpy(count+c5*h, x0+h*(a51*kx1+a52*kx2+a53*kx3+a54*kx4), y0+h*(a51*ky1+a52*ky2+a53*ky3+a54*ky4), px0+h*(a51*kpx1+a52*kpx2+a53*kpx3+a54*kpx4), py0+h*(a51*kpy1+a52*kpy2+a53*kpy3+a54*kpy4));
        kB5=fB(count+c5*h, x0+h*(a51*kx1+a52*kx2+a53*kx3+a54*kx4), y0+h*(a51*ky1+a52*ky2+a53*ky3+a54*ky4), px0+h*(a51*kpx1+a52*kpx2+a53*kpx3+a54*kpx4), py0+h*(a51*kpy1+a52*kpy2+a53*kpy3+a54*kpy4));

        var1=px0+h*(a61*kpx1+a62*kpx2+a63*kpx3+a64*kpx4+a65*kpx5);
        var2=py0+h*(a61*kpy1+a62*kpy2+a63*kpy3+a64*kpy4+a65*kpy5);
        kx6=fx(count+c6*h, x0+h*(a61*kx1+a62*kx2+a63*kx3+a64*kx4+a65*kx5), y0+h*(a61*ky1+a62*ky2+a63*ky3+a64*ky4+a65*ky5), var1, var2);
        ky6=fy(count+c6*h, x0+h*(a61*kx1+a62*kx2+a63*kx3+a64*kx4+a65*kx5), y0+h*(a61*ky1+a62*ky2+a63*ky3+a64*ky4+a65*ky5), var1, var2);
        kpx6=fpx(count+c6*h, x0+h*(a61*kx1+a62*kx2+a63*kx3+a64*kx4+a65*kx5), y0+h*(a61*ky1+a62*ky2+a63*ky3+a64*ky4+a65*ky5), var1, var2);
        kpy6=fpy(count+c6*h, x0+h*(a61*kx1+a62*kx2+a63*kx3+a64*kx4+a65*kx5), y0+h*(a61*ky1+a62*ky2+a63*ky3+a64*ky4+a65*ky5), var1, var2);
        kB6=fB(count+c6*h, x0+h*(a61*kx1+a62*kx2+a63*kx3+a64*kx4+a65*kx5), y0+h*(a61*ky1+a62*ky2+a63*ky3+a64*ky4+a65*ky5), var1, var2);

        var3=px0+h*(a71*kpx1+a72*kpx2+a73*kpx3+a74*kpx4+a75*kpx5+a76*kpx6);
        var4=py0+h*(a71*kpy1+a72*kpy2+a73*kpy3+a74*kpy4+a75*kpy5+a76*kpy6);
        kx7=fx(count+c7*h, x0+h*(a71*kx1+a72*kx2+a73*kx3+a74*kx4+a75*kx5+a76*kx6), y0+h*(a71*ky1+a72*ky2+a73*ky3+a74*ky4+a75*ky5+a76*ky6), var3, var4);
        ky7=fy(count+c7*h, x0+h*(a71*kx1+a72*kx2+a73*kx3+a74*kx4+a75*kx5+a76*kx6), y0+h*(a71*ky1+a72*ky2+a73*ky3+a74*ky4+a75*ky5+a76*ky6), var3, var4);
        kpx7=fpx(count+c7*h, x0+h*(a71*kx1+a72*kx2+a73*kx3+a74*kx4+a75*kx5+a76*kx6), y0+h*(a71*ky1+a72*ky2+a73*ky3+a74*ky4+a75*ky5+a76*ky6), var3, var4);
        kpy7=fpy(count+c7*h, x0+h*(a71*kx1+a72*kx2+a73*kx3+a74*kx4+a75*kx5+a76*kx6), y0+h*(a71*ky1+a72*ky2+a73*ky3+a74*ky4+a75*ky5+a76*ky6), var3, var4);
        kB7=fB(count+c7*h, x0+h*(a71*kx1+a72*kx2+a73*kx3+a74*kx4+a75*kx5+a76*kx6), y0+h*(a71*ky1+a72*ky2+a73*ky3+a74*ky4+a75*ky5+a76*ky6), var3, var4);

        x1=x0+h*(b1*kx1+b2*kx2+b3*kx3+b4*kx4+b5*kx5+b6*kx6+b7*kx7);
        x2=x0+h*(p1*kx1+p2*kx2+p3*kx3+p4*kx4+p5*kx5+p6*kx6+p7*kx7);
        errx=fabs(x2-x1);

        y1=y0+h*(b1*ky1+b2*ky2+b3*ky3+b4*ky4+b5*ky5+b6*ky6+b7*ky7);
        y2=y0+h*(p1*ky1+p2*ky2+p3*ky3+p4*ky4+p5*ky5+p6*ky6+p7*ky7);
        erry=fabs(y2-y1);

        px1=px0+h*(b1*kpx1+b2*kpx2+b3*kpx3+b4*kpx4+b5*kpx5+b6*kpx6+b7*kpx7);
        px2=px0+h*(p1*kpx1+p2*kpx2+p3*kpx3+p4*kpx4+p5*kpx5+p6*kpx6+p7*kpx7);
        errpx=fabs(px2-px1);

        py1=py0+h*(b1*kpy1+b2*kpy2+b3*kpy3+b4*kpy4+b5*kpy5+b6*kpy6+b7*kpy7);
        py2=py0+h*(p1*kpy1+p2*kpy2+p3*kpy3+p4*kpy4+p5*kpy5+p6*kpy6+p7*kpy7);
        errpy=fabs(py2-py1);

        B1=B0+h*(b1*kB1+b2*kB2+b3*kB3+b4*kB4+b5*kB5+b6*kB6+b7*kB7);
        B2=B0+h*(p1*kB1+p2*kB2+p3*kB3+p4*kB4+p5*kB5+p6*kB6+p7*kB7);
        errB=fabs(B2-B1);

        err=sqrt((pow(errx,2)+pow(erry,2)+(pow(errpx,2)+pow(errpy,2)+pow(errB,2))/5));
        //err=fabs(errx-erry)/31;

        if (err<=tol)
        {
            if (step==1)
            {
                fout<<count<<" "<<x0<<" "<<y0<<" "<<px0<<" "<<py0<<" "<<B0<<"\n";
            }

            if (((py0<-1.0 and py1<1.0 and py1>-1.0) or (py1<-1.0 and py0<1.0 and py0>-1.0) or (py0>1.0 and py1<1.0 and py1>-1.0) or (py1>1.0 and py0<1.0 and py0>-1.0) or (py0<-1.0 and py1>1.0) or (py0>1.0 and py1<-1.0)) and (num==0) and (m!=2) and (fabs(py0+1)>pow(10,-12)) and (fabs(py0-1)>pow(10,-12)))
            {
                //hor=method_chord(count, count+h, py0,py1, pow(10,-13));
                //x10=x0;
                //y10=y0;
                //px10=px0;
                //py10=py0;
                printf("%.25f \t", py0);
                cout<<py1<<"\t"<<count<<"\t"<<count+h<<endl;
                m=0;
                num=1;
                t_next=count+h;
                do
                {
                    t_curr=t_next;
                    if (py0+py1<-1)
                    {
                        t_next = t_next + (-1.0-py1)*(t_next-count)/(py1-py0);
                    }
                    if (py0+py1>1)
                    {
                        t_next = t_next + (1.0-py1)*(t_next-count)/(py1-py0);
                    }
                    cout<<t_next<<endl;
                    graph(count, t_next, tol, x0, y0, px0, py0, B0, data1);
                    py1=py;
                    cout<<py1<<endl;
                    printf("fabs=%.15f \n", fabs(t_next - t_curr));
                }while (fabs(t_next - t_curr) > pow(10,-13));
                num=0;
                per=t_next;
                //cout<<"per="<<per<<endl;
                printf("per= %.15f \n",per);
                coco=1;
                papa++;

            }
            cout<<"papa="<<papa<<endl;
            if (coco!=1)
            {
                x0=x1;
                y0=y1;
                px0=px1;
                py0=py1;
                B0=B1;

                count+=h;
                n++;
            }
            if(T-count<epsilon)
            {
                x=x0;
                y=y0;
                px=px0;
                py=py0;
                B=B0;
                break;
            }
        }
        if (coco!=1)
        {
            h=h*fmin(facmax, fmax(facmin, fac*pow(tol/err,1.0/6)));
            m=0;
        }
        else
        {
            cout<<count<<endl;
            h=per-count;
            m=2;
        }
        coco=0;

    }
    if (step==1)
    {
        fout <<"\n";
    }
    //cout <<n<<endl;
    //x=x0;
    //y=y0;
    //fout << x0 << " " << y0 << "\n";
    //return n;
}

/*double method_chord(double x_prev, double x_curr, p1, p2, double e)
{
    double x_next = 0;
    double tmp;
    double t1=0, t2=0, p01=0;
    t1=x_prev;
    p01=p1;

    do
    {
        tmp = x_next;
        x_next = x_curr - p2 * (x_prev - x_curr) / (p1 - p2);
        x_prev = x_curr;
        x_curr = tmp;
        p1=p2;
        graph(t1,x_next,tol,)

    } while (abs(x_next - x_curr) > e);

    return x_next;
}*/





int main()
{
    int N=0;
    double /*a=0.0, b=0.0,*/ x0=0, y0=0, x11=0, x12=0, y11=0, y12=0, T=10.0;
    double alpha1=36.1916, alpha2=78.5923;
    double alpha10=0, alpha20=0;
    double X1=0, X2=0;
    double X11=0, X12=0, X21=0, X22=0;
    double h=1.0;
    ofstream data("data.txt");

    graph(0,T,tol,0,1,alpha1,alpha2,0,data);
    X1=x;
    X2=py;
    cout << X1 << " " << X2 << endl;
    //while (fmax(fabs(X1),fabs(X2))>epsilon)
    while (sqrt(pow(X1,2)+pow(X2,2))>epsilon)
    {
        step++;
        cout<<"step="<<step<<"\n";

        x0=X1;
        y0=X2;
        graph(0,T,tol,0,1,alpha1+delta,alpha2,0,data);
        x11=x;
        //printf("%.15f \n",x11);
        //cout <<x11<<endl;
        y11=py;
        //cout <<y11<<endl;
        X11=(x11-x0)/delta;
        X21=(y11-y0)/delta;
        graph(0,T,tol,0,1,alpha1,alpha2+delta,0,data);
        x12=x;
        //cout <<x12<<endl;
        //printf("%.15f \n",x12);
        y12=py;
        //cout <<y12<<endl;
        X12=(x12-x0)/delta;
        X22=(y12-y0)/delta;
        cout <<X11 << " " << X12 << " " <<X21 << " " <<X22 <<endl;
        alpha1=alpha1+h*(-x0*X22+y0*X12)/(X11*X22-X12*X21);//По методу Крамера
        alpha10=alpha1;
        alpha2=alpha2+h*(-y0*X11+x0*X21)/(X11*X22-X12*X21);
        alpha20=alpha2;
        cout <<alpha1 << " " << alpha2 <<endl;
        graph(0,T,tol,0,1,alpha1,alpha2,0,data);
        X1=x;
        X2=py;
        N++;
        while (sqrt(pow(X1,2)+pow(X2,2)) > (sqrt(pow(x0,2)+pow(y0,2))) and h>pow(2,-30))
        {

            h=h/2.0;
            alpha1=alpha10+h*(-x0*X22+y0*X12)/(X11*X22-X12*X21);//По методу Крамера
            alpha2=alpha20+h*(-y0*X11+x0*X21)/(X11*X22-X12*X21);
            //N++;
            graph(0,T,tol,0,1,alpha1,alpha2,0,data);
            X1=x;
            X2=py;

        }
        h=1.0;

    }
    //cout<<alpha1<<"\t"<<alpha2<<endl;
    printf("alpha1= %.15f \t alpha2= %.15f \n", alpha1, alpha2);
    cout <<"Number of iterations="<<N<<"\n";
    cout<<step<<"\n";
    //cout<<B<<"\n";
    printf("B0= %.15f \n", B);
    return 0;
}
