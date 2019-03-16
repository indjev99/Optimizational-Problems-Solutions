#include<iostream>
#include<fstream>
#include<vector>
#include<queue>
#include<stack>
#include<algorithm>
#include<math.h>
#include<unordered_map>
#include<conio.h>
#include<chrono>
using namespace std;
using namespace chrono;

#define official

#ifdef official
#define cin inF
#define cout outF
#endif

ifstream inF("equations.in");
ofstream outF("equations.out");

const int MAX_M=2000;
const int MAX_N=1000;
const int MAX_S=20;
const long long MAX_X=1e6;
const long long MAX_X1=1e6;
const long long MAX_X_GEN=1e4;
const int TRIES=2000;

const double EPS=1e-6;

const int EXACT_SCALE=0;

int n,m,m2,s;

struct Equation
{
    long long a[MAX_N+100];
    long long b[MAX_S+10];
    long long aprox;
    long long exact;
    int num;
    void add (Equation& other, long long scalar)
    {
        for (int i=0; i<n; ++i)
        {
            a[i]+=other.a[i]*scalar;
        }
        for (int i=0; i<s; ++i)
        {
            b[i]+=other.b[i]*scalar;
        }
    }
    void scale (long long scalar)
    {
        for (int i=0; i<n; ++i)
        {
            a[i]*=scalar;
        }
        for (int i=0; i<s; ++i)
        {
            b[i]*=scalar;
        }
    }
};

bool cmp(Equation a, Equation b)
{
    return a.aprox+a.exact>b.aprox+b.exact;
}

int big_rand()
{
    int r=rand()*RAND_MAX;
    return r+rand();
}

Equation eqs[MAX_M];

long long xGCD(long long a, long long b, long long &x, long long &y)
{
    long long lastx,lasty,temp,q;
    bool an=0,bn=0,sw=0;
    if (a<0)
    {
        a=-a;
        an=1;
    }
    if (b<0)
    {
        b=-b;
        bn=1;
    }
    if (a<b)
    {
        swap(a,b);
        sw=1;
    }
    x=0;
    y=1;
    lastx=1;
    lasty=0;
    while (b!=0)
    {
        q=a/b;
        temp=a%b;
        a=b;
        b=temp;

        temp=x;
        x=lastx-q*x;
        lastx =temp;

        temp=y;
        y=lasty-q*y;
        lasty=temp;
    }
    x=lastx;
    y=lasty;
    if (sw) swap(x,y);
    if (an) x=-x;
    if (bn) y=-y;
    return a;
}
long long GCD(long long a, long long b)
{
    long long temp;
    if (a<0)
    {
        a=-a;
    }
    if (b<0)
    {
        b=-b;
    }
    if (a<b)
    {
        swap(a,b);
    }
    while (b!=0)
    {
        temp=a%b;
        a=b;
        b=temp;
    }
    return a;
}
long long mod(long long a, long long m)
{
    bool n=0;
    if (a<0)
    {
        a=-a;
        n=1;
    }
    a%=m;
    if (n) a=m-a;
}

struct Sol
{
    int chm;
    long long x[MAX_N+100];
    long long sums[MAX_M+100];
    double score=0;
    void change(int ind, long long nx)
    {
        long long diff=nx-x[ind];
        x[ind]=nx;
        score=0;
        for (int i=0; i<m; ++i)
        {
            sums[i]+=diff*eqs[i].a[ind];
            if (sums[i]==eqs[i].b[chm]) score+=eqs[i].exact;
            score+=(1.0/(fabs(sums[i]-eqs[i].b[chm])+1.0))*eqs[i].aprox;
        }
    }
    void eval()
    {
        //cerr<<"eval"<<endl;
        long long sum;
        score=0;
        for (int i=0; i<m; ++i)
        {
            sum=0;
            for (int j=0; j<n; ++j)
            {
                sum+=round(x[j])*eqs[i].a[j];
            }
            //cerr<<sum<<" "<<eqs[i].b[chm]<<" "<<1.0/(fabs(sum-eqs[i].b[chm])+1.0)*eqs[i].aprox<<'\n';
            sums[i]=sum;
            if (sum==eqs[i].b[chm]) score+=eqs[i].exact;
            score+=(1.0/(fabs(sum-eqs[i].b[chm])+1.0))*eqs[i].aprox;
        }
        //cerr<<endl;
    }
};

int chm;
int x[MAX_N+100];

void input()
{
    cin>>m>>n;
    for (int i=0; i<m; ++i)
    {
        eqs[i].num=i;
        for (int j=0; j<n; ++j)
        {
            cin>>eqs[i].a[j];
        }
        cin>>eqs[i].exact;
        cin>>eqs[i].aprox;
    }
    cin>>s;
    for (int i=0; i<s; ++i)
    {
        for (int j=0; j<m; ++j)
        {
            cin>>eqs[j].b[i];
        }
    }
}

void output()
{
    cout<<chm<<'\n';
    for (int i=0; i<n; ++i)
    {
        if (i) cout<<' ';
        cout<<x[i];
    }
    cout<<'\n';
}

double bps=0;

void find_bps()
{
    bps=0;
    for (int i=0; i<m; ++i)
    {
        bps+=eqs[i].exact+eqs[i].aprox;
    }
}

high_resolution_clock::time_point startT,currT,start2T;
double reserve=0;

Sol maxSol;
Sol curr;
bool succ=0;

void try_solve1(Equation& eq)
{
    long long x,y,z,q,xh,yh;
    long long currcoeff;
    long long b,a1,a2,b0;
    long long gcd;
    int p3,p4;
    currT=high_resolution_clock::now();
    for (int p1=0; p1<n-1 && (!succ || duration_cast<duration<double>>(currT-startT).count()<4); ++p1)
    {
        a1=eq.a[p1];
        if (a1==0) continue;
        gcd=a1;
        xh=1;
        for (int i=0; i<s; ++i)
        {
            b0=eq.b[i];
            b=b0;
            while (1)
            {
                if (b%gcd) break;
                x=xh*(b/gcd);
                if (x>MAX_X || x<-MAX_X) break;
                for (int j=0; j<n; ++j)
                {
                    curr.x[j]=0;
                }
                curr.x[p1]=x;
                curr.chm=i;
                curr.eval();
                if (curr.score>maxSol.score) maxSol=curr;
                //cerr<<endl;
                //cerr<<maxSol.chm<<endl;
                //cerr<<"1. x["<<p1<<"] = "<<double(maxSol.x[p1])<<endl;
                succ=1;
                break;
            }
            for (p3=0; p3<n; ++p3)
            {
                p4=rand()%n;
                if (p3==p1 || p4==p1 || p4==p3) continue;
                z=big_rand()%(2*MAX_X_GEN+1)-MAX_X_GEN;
                q=big_rand()%(2*MAX_X_GEN+1)-MAX_X_GEN;
                b=b0-eq.a[p3]*z-eq.a[p4]*q;
                while (1)
                {
                    if (b%gcd) break;
                    x=xh*(b/gcd);
                    if (x>MAX_X || x<-MAX_X) break;
                    for (int j=0; j<n; ++j)
                    {
                        curr.x[j]=0;
                    }
                    curr.x[p1]=x;
                    curr.x[p3]=z;
                    curr.x[p4]=q;
                    curr.chm=i;
                    curr.eval();
                    if (curr.score>maxSol.score) maxSol=curr;
                    //cerr<<endl;
                    //cerr<<maxSol.chm<<endl;
                    //cerr<<"1. x["<<p1<<"] = "<<double(maxSol.x[p1])<<endl;
                    succ=1;
                    break;
                }
            }
        }
        currT=high_resolution_clock::now();
        for (int p2=p1+1; p2<n && (!succ || duration_cast<duration<double>>(currT-startT).count()<4.5-2*reserve); ++p2)
        {
            start2T=high_resolution_clock::now();
            a2=eq.a[p2];
            if (a2==0) continue;
            gcd=xGCD(a1,a2,xh,yh);
            if (xh>MAX_X1 || xh<-MAX_X1 || yh>MAX_X1 || yh<-MAX_X1) break;
            //cerr<<double(xh)<<" "<<double(yh)<<" "<<gcd<<'\n';
            for (int i=0; i<s; ++i)
            {
                b0=eq.b[i];
                b=b0;
                while (1)
                {
                    if (b%gcd) break;
                    x=xh*(b/gcd);
                    y=yh*(b/gcd);
                    //cerr<<double(x)<<" "<<double(y)<<" "<<gcd<<'\n';
                    if (x>MAX_X || x<-MAX_X || y>MAX_X || y<-MAX_X) break;
                    for (int j=0; j<n; ++j)
                    {
                        curr.x[j]=0;
                    }
                    curr.x[p1]=x;
                    curr.x[p2]=y;
                    curr.chm=i;
                    curr.eval();
                    if (curr.score>maxSol.score) maxSol=curr;
                    //cerr<<endl;
                    //cerr<<maxSol.chm<<endl;
                    //cerr<<"p1 x["<<p1<<"] = "<<double(x)<<" "<<double(x*eq.a[p1])<<endl;
                    //cerr<<"p2 x["<<p2<<"] = "<<double(y)<<" "<<double(y*eq.a[p2])<<endl;
                    //cerr<<"sum: "<<double(x*eq.a[p1]+y*eq.a[p2])<<" b = "<<double(b0)<<endl;
                    succ=1;
                    break;
                }
                for (p3=0; p3<n; ++p3)
                {
                    p4=rand()%n;
                    if (p3==p1 || p3==p2 || p4==p1 || p4==p2 || p4==p3) continue;
                    z=big_rand()%(2*MAX_X_GEN+1)-MAX_X_GEN;
                    q=big_rand()%(2*MAX_X_GEN+1)-MAX_X_GEN;
                    b=b0-eq.a[p3]*z-eq.a[p4]*q;
                    while (1)
                    {
                        if (b%gcd) break;
                        x=xh*(b/gcd);
                        y=yh*(b/gcd);
                        //cerr<<double(x)<<" "<<double(y)<<" "<<gcd<<'\n';
                        if (x>MAX_X || x<-MAX_X || y>MAX_X || y<-MAX_X) break;
                        for (int j=0; j<n; ++j)
                        {
                            curr.x[j]=0;
                        }
                        curr.x[p1]=x;
                        curr.x[p2]=y;
                        curr.x[p3]=z;
                        curr.x[p4]=q;
                        curr.chm=i;
                        curr.eval();
                        if (curr.score>maxSol.score) maxSol=curr;
                        //cerr<<endl;
                        //cerr<<maxSol.chm<<endl;
                        //cerr<<"p1 x["<<p1<<"] = "<<double(x)<<" "<<double(x*eq.a[p1])<<endl;
                        //cerr<<"p2 x["<<p2<<"] = "<<double(y)<<" "<<double(y*eq.a[p2])<<endl;
                        //cerr<<"p3 x["<<p3<<"] = "<<double(z)<<" "<<double(z*eq.a[p3])<<endl;
                        //cerr<<"sum: "<<double(x*eq.a[p1]+y*eq.a[p2]+z*eq.a[p3])<<" b = "<<double(b0)<<endl;
                        //cerr<<"sumxy: "<<double(x*eq.a[p1]+y*eq.a[p2])<<" b = "<<double(b)<<endl;
                        //cerr<<double(b0-eq.a[p3]*z)<<" "<<double(b)<<endl;
                        //cerr<<double(b+eq.a[p3]*z)<<" "<<double(b0)<<endl;
                        succ=1;
                        break;
                    }
                }
            }
            currT=high_resolution_clock::now();
            reserve=duration_cast<duration<double>>(currT-start2T).count();
        }
        currT=high_resolution_clock::now();
    }
}
void smart_solve()
{
    sort(eqs,eqs+m,cmp);
    maxSol.score=0;
    succ=0;

    for (int i=0; i<m && !succ; ++i)
    {
        try_solve1(eqs[i]);
    }
    cerr<<maxSol.score<<" out of "<<bps<<endl;
}

void solve()
{
    find_bps();
    smart_solve();

    chm=maxSol.chm+1;
    for (int i=0; i<n; ++i)
    {
        x[i]=maxSol.x[i];
    }
}

int main()
{
    startT=high_resolution_clock::now();
    currT=high_resolution_clock::now();
    srand(0);

    input();
    solve();
    output();

    currT=high_resolution_clock::now();

    cerr<<"\n TOTAL TIME: "<<duration_cast<duration<double>>(currT-startT).count()<<endl;

    return 0;
}
