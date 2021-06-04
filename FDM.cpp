#include <iostream>
#include <math.h>
#include <fstream>

double zero=0.000001;
double tau=0.01;
double h=0.001;

using namespace std;


double f0(double x)
{
    if (x<=zero) return 0.0;
    else return 1.0;
}

double sol (double t, double x)
{
    if (x<=zero) return 0;
    if (x>=t) return 1;
    else return x/t;
}

int main()
{
    double tmin=0.0, tmax=1.0, xmin=-0.5, xmax=1.5, u, t, x, max1=-99, max2=-99, sum1=0, sum2=0;
    double tau1, tau2, tau3, tau4, h1, h2, h3, h4;
    int m,n,M,N;
    ofstream fout("data4.txt");
    n=(tmax-tmin)/tau+1;
    m=(xmax-xmin)/h+1;
    M=m;
    double *low = new double [m];
    double *up = new double [m];

    for (int i=0; i<m; i++)
    {
        low[i]=f0(xmin+i*h);
    }

    for (int i=0; i<n; i++)
    {
        up[0]=low[0];
        up[m-1]=low[m-1];

        for (int j=0; j<m-2; j++)
        {
            up[j+1]=sqrt(pow(up[j],2)-2*h/tau*(up[j]-low[j]));
        }
        for (int j=0; j<m; j++)
        {
            low[j]=up[j];
        }
    }
    for (int j=0; j<m; j++)
    {
        fout << xmin+j*h << " " << up[j] << "\n";
    }


    t=tmin;
    x=xmin;
    ofstream file("datau4.txt");
    ofstream delta("delta4.txt");
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            u=sol(t,x);
            if (i==n-1)
            file << x << " " << u << "\n";
            x=x+h;
        }
        t=t+tau;
        x=xmin;
    }

    for (int j=0; j<m; j++)
    {
        sum1=sum1+abs(up[j]-sol(1.0,xmin+h*j));
        delta<<xmin+h*j<<" "<<up[j]-sol(1.0,xmin+h*j)<<"\n";
        sum2=sum2+abs(up[j]);
        if (abs(up[j]-sol(1.0,xmin+h*j))>max1) max1=abs(up[j]-sol(1.0,xmin+h*j));
        if (abs(up[j]>max2)) max2=abs(up[j]);
    }
    sum1=sum1*h;
    sum2=sum2*h;
    cout<<"tau="<<tau<<" h="<<h<<"\n";
    cout<<max1<<" "<<sum1<<" "<<max1/max2<<" "<<sum1/sum2<<endl;


    sum1=0;
    max1=-99;
    tau1=tau/2;
    h1=h/2;
    n=(tmax-tmin)/tau1+1;
    m=(xmax-xmin)/h1+1;
    double *low1 = new double [m];
    double *up1 = new double [m];

    for (int i=0; i<m; i++)
    {
        low1[i]=f0(xmin+i*h1);
    }

    for (int i=0; i<n; i++)
    {
        up1[0]=low1[0];
        up1[m-1]=low1[m-1];

        for (int j=0; j<m-2; j++)
        {
            up1[j+1]=sqrt(pow(up1[j],2)-2*h/tau*(up1[j]-low1[j]));
        }
        for (int j=0; j<m; j++)
        {
            low1[j]=up1[j];
        }
    }

    for (int j=0; j<M; j++)
    {
        sum1=sum1+abs(up[j]-up1[2*j]);
        if (abs(up[j]-up1[2*j])>max1) max1=abs(up[j]-up1[2*j]);
    }
    sum1=sum1*h;
    cout<<max1<<" "<<sum1<<" "<<max1/max2<<" "<<sum1/sum2<<endl;


    sum1=0;
    max1=-99;
    tau2=tau/4;
    h2=h/4;
    n=(tmax-tmin)/tau2+1;
    m=(xmax-xmin)/h2+1;
    double *low2 = new double [m];
    double *up2 = new double [m];

    for (int i=0; i<m; i++)
    {
        low2[i]=f0(xmin+i*h2);
    }

    for (int i=0; i<n; i++)
    {
        up2[0]=low2[0];
        up2[m-1]=low2[m-1];

        for (int j=0; j<m-2; j++)
        {
            up2[j+1]=sqrt(pow(up2[j],2)-2*h/tau*(up2[j]-low2[j]));
        }
        for (int j=0; j<m; j++)
        {
            low2[j]=up2[j];
        }
    }

    for (int j=0; j<M; j++)
    {
        sum1=sum1+abs(up[j]-up2[4*j]);
        if (abs(up[j]-up2[4*j])>max1) max1=abs(up[j]-up2[4*j]);
    }
    sum1=sum1*h;
    cout<<max1<<" "<<sum1<<" "<<max1/max2<<" "<<sum1/sum2<<endl;


    sum1=0;
    max1=-99;
    tau3=tau/8;
    h3=h/8;
    n=(tmax-tmin)/tau3+1;
    m=(xmax-xmin)/h3+1;
    double *low3 = new double [m];
    double *up3 = new double [m];

    for (int i=0; i<m; i++)
    {
        low3[i]=f0(xmin+i*h3);
    }

    for (int i=0; i<n; i++)
    {
        up3[0]=low3[0];
        up3[m-1]=low3[m-1];

        for (int j=0; j<m-2; j++)
        {
            up3[j+1]=sqrt(pow(up3[j],2)-2*h/tau*(up3[j]-low3[j]));
        }
        for (int j=0; j<m; j++)
        {
            low3[j]=up3[j];
        }
    }

    for (int j=0; j<M; j++)
    {
        sum1=sum1+abs(up[j]-up3[8*j]);
        if (abs(up[j]-up3[8*j])>max1) max1=abs(up[j]-up3[8*j]);
    }
    sum1=sum1*h;
    cout<<max1<<" "<<sum1<<" "<<max1/max2<<" "<<sum1/sum2<<endl;


    sum1=0;
    max1=-99;
    tau4=tau/16;
    h4=h/16;
    n=(tmax-tmin)/tau4+1;
    m=(xmax-xmin)/h4+1;
    double *low4 = new double [m];
    double *up4 = new double [m];

    for (int i=0; i<m; i++)
    {
        low4[i]=f0(xmin+i*h4);
    }

    for (int i=0; i<n; i++)
    {
        up4[0]=low1[0];
        up4[m-1]=low4[m-1];

        for (int j=0; j<m-2; j++)
        {
            up4[j+1]=sqrt(pow(up4[j],2)-2*h/tau*(up4[j]-low4[j]));
        }
        for (int j=0; j<m; j++)
        {
            low4[j]=up4[j];
        }
    }

    for (int j=0; j<M; j++)
    {
        sum1=sum1+abs(up[j]-up4[16*j]);
        if (abs(up[j]-up4[16*j])>max1) max1=abs(up[j]-up4[16*j]);
    }
    sum1=sum1*h;
    cout<<max1<<" "<<sum1<<" "<<max1/max2<<" "<<sum1/sum2<<endl;


    return 0;

}



