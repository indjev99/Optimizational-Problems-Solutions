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

ifstream inF("drawing.in");
ofstream outF("drawing.out");

const int MAX_N=512;
const int MAX_K=16;
int RES=900;
int RES0=0;
const int MAX_RES=1024;
int STEPS=5;
const int MAX_STEPS=16;
const double EPS=1e-12;
double POW1=1.5;//2;
double POW2=3;//4;
int PROB1=3;

int startPoint=-1;

struct Color
{
    double r,g,b;
};
bool operator<(const Color& a, const Color& b)
{
    return a.r+a.g+b.g>b.r+b.g+b.b;
}
Color operator+(const Color& a, const Color& b)
{
    Color r;
    r.r=a.r+b.r;
    r.g=a.g+b.g;
    r.b=a.b+b.b;
    return r;
}
Color operator*(const Color& a, double b)
{
    Color r;
    r.r=a.r*b;
    r.g=a.g*b;
    r.b=a.b*b;
    return r;
}
Color operator/(const Color& a, double b)
{
    Color r;
    r.r=a.r/b;
    r.g=a.g/b;
    r.b=a.b/b;
    return r;
}
Color floor(const Color& a)
{
    Color r;
    r.r=floor(a.r);
    r.g=floor(a.g);
    r.b=floor(a.b);
    return r;
}
double diff2(const Color& a, const Color& b)
{
    double d=(a.r-b.r)*(a.r-b.r)+(a.g-b.g)*(a.g-b.g)+(a.b-b.b)*(a.b-b.b);
    return d;
}
double diff(const Color& a, const Color& b)
{
    return sqrt(diff2(a,b));
}

struct ResultingColor
{
    Color c;
    int type;
    int steps[MAX_STEPS*2];
};
bool operator<(const ResultingColor& a, const ResultingColor& b)
{
    return a.c<b.c;
}
struct Unit
{
    int x,y;
    int type;
    double diff;
};
bool operator<(const Unit& a, const Unit& b)
{
    return a.diff<b.diff;
}

int n,k;
Color p[MAX_K];
Color d[MAX_N][MAX_N];
ResultingColor r[MAX_RES];
vector<Unit> units;
vector<Unit> units2;
vector<Unit> inGroup[MAX_RES];
int wx,wy,ox,oy;
int wi=-1;
bool used[MAX_N][MAX_N];

int OPERATIONS=0;

void input()
{
    cin>>n>>k;
    for (int i=0; i<k; ++i)
    {
        cin>>p[i].r>>p[i].g>>p[i].b;
    }
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<n; ++j)
        {
            cin>>d[i][j].r;
        }
    }
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<n; ++j)
        {
            cin>>d[i][j].g;
        }
    }
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<n; ++j)
        {
            cin>>d[i][j].b;
        }
    }
}
void getResulting1(int num, const Color& curr, int x, int y)
{
    int mj;
    double cd;
    double pr[MAX_K];
    double cpr;
    double tp;
    Color cc=curr/(1<<STEPS);
    Color ccc;
    double coeff=1;
    double total=1.0/(1<<STEPS);
    r[num].type=1;
    for (int i=STEPS-1; i>=0; --i)
    {
        coeff*=2;
        tp=0;
        total+=1.0/coeff;
        for (int j=startPoint; j<k; ++j)
        {
            if (j==-1) ccc= {255,255,255};
            else ccc=p[j];
            ccc=ccc/coeff;
            cd=diff2(d[x][y],(cc+ccc)/total);
            pr[j+1]=100000.0/pow(cd+1,POW1/2);
            tp+=pr[j+1];
        }
        for (int j=startPoint; j<k; ++j)
        {
            cpr=pr[j+1]/tp;
            if (rand()<=cpr*(RAND_MAX+1))
            {
                mj=j;
                break;
            }
            tp-=pr[j+1];
        }
        r[num].steps[i]=mj;
        if (mj==-1) ccc= {255,255,255};
        else ccc=p[mj];
        ccc=ccc/coeff;
        cc=cc+ccc;
    }
    //cerr<<cc.r<<" "<<cc.g<<" "<<cc.b<<endl;
    r[num].c=curr;
    for (int i=0; i<STEPS; ++i)
    {
        ++OPERATIONS;
        //cerr<<r[num].steps[i]<<" ";
        if (r[num].steps[i]==-1) ccc= {255,255,255};
        else ccc=p[r[num].steps[i]];
        r[num].c=floor((r[num].c+ccc)/2);
    }
}
void getResulting2(int num, const Color& curr, int x, int y)
{
    int mj;
    double cd;
    double pr[MAX_K];
    double cpr;
    double tp;
    Color cc=curr/(1<<STEPS);
    Color ccc;
    double coeff=2;
    double total=1.0/(1<<STEPS);
    r[num].type=2;
    for (int i=STEPS-1; i>=0; --i)
    {
        coeff*=2;
        tp=0;
        total+=1.0/coeff;
        for (int j=startPoint; j<k; ++j)
        {
            if (j==-1) ccc= {255,255,255};
            else ccc=p[j];
            ccc=ccc/coeff;
            cd=diff2(d[x][y],(cc+ccc)/total);
            pr[j+1]=100000.0/pow(cd+1,POW2/2);
            tp+=pr[j+1];
        }
        for (int j=startPoint; j<k; ++j)
        {
            cpr=pr[j+1]/tp;
            if (rand()<=cpr*(RAND_MAX+1))
            {
                mj=j;
                break;
            }
            tp-=pr[j+1];
        }
        r[num].steps[i]=mj;
        if (mj==-1) ccc= {255,255,255};
        else ccc=p[mj];
        ccc=ccc/coeff;
        cc=cc+ccc;

        tp=0;
        total+=1.0/coeff;
        for (int j=startPoint; j<k; ++j)
        {
            if (j==-1) ccc= {255,255,255};
            else ccc=p[j];
            ccc=ccc/coeff;
            cd=diff2(d[x][y],(cc+ccc)/total);
            pr[j+1]=100000.0/pow(cd+1,POW2/2);
            tp+=pr[j+1];
        }
        for (int j=startPoint; j<k; ++j)
        {
            cpr=pr[j+1]/tp;
            if (rand()<=cpr*(RAND_MAX+1))
            {
                mj=j;
                break;
            }
            tp-=pr[j+1];
        }
        r[num].steps[STEPS+i]=mj;
        if (mj==-1) ccc= {255,255,255};
        else ccc=p[mj];
        ccc=ccc/coeff;
        cc=cc+ccc;
    }
    //cerr<<cc.r<<" "<<cc.g<<" "<<cc.b<<endl;
    r[num].c=curr;
    for (int i=0; i<STEPS; ++i)
    {
        ++OPERATIONS;
        //cerr<<r[num].steps[i]<<" ";
        if (r[num].steps[i]==-1) ccc= {255,255,255};
        else ccc=p[r[num].steps[i]];
        r[num].c=floor((r[num].c+ccc)/2);
    }
    ++OPERATIONS;
    cc=r[num].c;
    for (int i=0; i<STEPS; ++i)
    {
        ++OPERATIONS;
        //cerr<<r[num].steps[i]<<" ";
        if (r[num].steps[STEPS+i]==-1) ccc= {255,255,255};
        else ccc=p[r[num].steps[STEPS+i]];
        cc=floor((cc+ccc)/2);
    }
    ++OPERATIONS;
    r[num].c=floor((cc+r[num].c)/2);
}
Color getResultingForAll(int num, const Color& curr)
{
    int cnum=num;
    for (int i=0; i<STEPS; ++i)
    {
        if (cnum%2==0) r[num].steps[i]=wi;
        else r[num].steps[i]=1-wi;
        cnum/=2;
    }

    Color ccc;
    r[num].type=1;
    r[num].c=curr;
    for (int i=0; i<STEPS; ++i)
    {
        ++OPERATIONS;
        //cerr<<r[num].steps[i]<<" ";
        if (r[num].steps[i]==-1) ccc= {255,255,255};
        else ccc=p[r[num].steps[i]];
        r[num].c=floor((r[num].c+ccc)/2);
    }
    return r[num].c;
}
Color getResultingForAll2(int num, const Color& curr)
{
    int cnum=num;
    for (int i=0; i<STEPS; ++i)
    {
        if (cnum%3==0) r[num].steps[i]=startPoint;
        else if (cnum%3==1) r[num].steps[i]=startPoint+1;
        else r[num].steps[i]=startPoint+2;
        cnum/=3;
    }

    Color ccc;
    r[num].type=1;
    r[num].c=curr;
    for (int i=0; i<STEPS; ++i)
    {
        ++OPERATIONS;
        //cerr<<r[num].steps[i]<<" ";
        if (r[num].steps[i]==-1) ccc= {255,255,255};
        else ccc=p[r[num].steps[i]];
        r[num].c=floor((r[num].c+ccc)/2);
    }
    return r[num].c;
}
int doAll=0;
Color getResulting(int num, const Color& curr)
{
    if (doAll==1)
    {
        return getResultingForAll(num,curr);
    }
    if (doAll==2)
    {
        return getResultingForAll2(num,curr);
    }

    int x,y;
    x=rand()%n;
    y=rand()%n;

    if (PROB1 && rand()%PROB1==0) getResulting1(num,curr,x,y);
    else getResulting2(num,curr,x,y);
    //getResulting2(num,curr,x,y);

    /*cerr<<endl;
    cerr<<d[x][y].r<<" "<<d[x][y].g<<" "<<d[x][y].b<<endl;
    cerr<<d[x+1][y].r<<" "<<d[x+1][y].g<<" "<<d[x+1][y].b<<endl;
    cerr<<r[num].c.r<<" "<<r[num].c.g<<" "<<r[num].c.b<<endl;
    cerr<<endl;*/
    //getch();
    return r[num].c;
}

double avgdiff=0;
int zerodiffs=0;

void findUnits1()
{
    double currdiff;
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<n; ++j)
        {
            if (i<n-1)
            {
                currdiff=diff2(d[i][j],d[i+1][j]);
                if (currdiff<EPS) ++zerodiffs;
                avgdiff+=currdiff;
                units2.push_back({i,j,1,currdiff});
            }

            if (j<n-1)
            {
                currdiff=diff2(d[i][j],d[i][j+1]);
                if (currdiff<EPS) ++zerodiffs;
                avgdiff+=currdiff;
                units2.push_back({i,j,2,currdiff});
            }
        }
    }
    avgdiff/=2*n*(n-1);
    avgdiff=sqrt(avgdiff);
}
void findUnits2()
{
    sort(units2.begin(),units2.end());
    int cnt=0;
    int x,y,x2,y2;
    for (int i=0; i<units2.size() && cnt<OPERATIONS; ++i)
    {
        auto& a=units2[i];
        x=a.x;
        y=a.y;
        if (a.type==1)
        {
            x2=x+1;
            y2=y;
        }
        else
        {
            x2=x;
            y2=y+1;
        }
        if (used[x][y] || used[x2][y2]) continue;
        ++cnt;
        used[x][y]=1;
        used[x2][y2]=1;
        units.push_back(a);
    }
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<n; ++j)
        {
            if (used[i][j]==0) units.push_back({i,j,0,0});
        }
    }
}

void findGroup(int num)
{
    int x,y,x2,y2,type;
    x=units[num].x;
    y=units[num].y;
    type=units[num].type;
    double md=-1,cd;
    int mg;
    if (type!=0)
    {
        if (type==1)
        {
            x2=x+1;
            y2=y;
        }
        else if (type==2)
        {
            x2=x;
            y2=y+1;
        }
        for (int i=0; i<RES; ++i)
        {
            cd=diff(d[x][y],r[i].c)+diff(d[x2][y2],r[i].c);
            if (md<0 || cd<md)
            {
                md=cd;
                mg=i;
            }
        }
    }
    else
    {
        for (int i=0; i<RES; ++i)
        {
            cd=diff(d[x][y],r[i].c);
            if (md<0 || cd<md)
            {
                md=cd;
                mg=i;
            }
        }
    }
    inGroup[mg].push_back(units[num]);
}

Color last;
int lastGroup;
void findWhite()
{
    for (int i=RES-1; i>=0; --i)
    {
        if (!inGroup[i].empty())
        {
            lastGroup=i;
            last=r[i].c;
            wx=inGroup[i][0].x;
            wy=inGroup[i][0].y;
            if (inGroup[i].size()>1)
            {
                ox=inGroup[i][1].x;
                oy=inGroup[i][1].y;
            }
            else if (inGroup[i][0].type!=0)
            {
                if (inGroup[i][0].type==1)
                {
                    ox=wx+1;
                    oy=wy;
                }
                else
                {
                    ox=wx;
                    oy=wy+1;
                }
            }
            else
            {
                for (int j=i-1; j>=0; --j)
                {
                    if (!inGroup[j].empty())
                    {
                        ox=inGroup[i][0].x;
                        oy=inGroup[i][0].y;
                        break;
                    }
                }
            }
            break;
        }
    }
}
void findOperations()
{
    OPERATIONS=units.size();
    for (int i=0; i<=lastGroup; ++i)
    {
        //if (!inGroup[i].empty())
        {
            if (r[i].type==1) OPERATIONS+=STEPS;
            else OPERATIONS+=2*STEPS+2;
        }
    }
}
void outputGroup(int num)
{
    //if (inGroup[num].empty()) return;
    if (r[num].type==1)
    {
        for (int i=0; i<STEPS; ++i)
        {
            if (r[num].steps[i]==-1) cout<<"2 "<<wx<<' '<<wy<<'\n';
            else cout<<"1 "<<r[num].steps[i]<<'\n';
        }
    }
    else
    {
        for (int i=0; i<STEPS; ++i)
        {
            if (r[num].steps[i]==-1) cout<<"2 "<<wx<<' '<<wy<<'\n';
            else cout<<"1 "<<r[num].steps[i]<<'\n';
        }
        cout<<"3 "<<ox<<' '<<oy<<' '<<ox<<' '<<oy<<'\n';
        for (int i=0; i<STEPS; ++i)
        {
            if (r[num].steps[STEPS+i]==-1) cout<<"2 "<<wx<<' '<<wy<<'\n';
            else cout<<"1 "<<r[num].steps[STEPS+i]<<'\n';
        }
        cout<<"2 "<<ox<<' '<<oy<<'\n';
    }
    int x,y,x2,y2,type;
    for (int i=0; i<inGroup[num].size(); ++i)
    {
        x=inGroup[num][i].x;
        y=inGroup[num][i].y;
        type=inGroup[num][i].type;
        if (type==0)
        {
            x2=x;
            y2=y;
        }
        else if (type==1)
        {
            x2=x+1;
            y2=y;
        }
        else if (type==2)
        {
            x2=x;
            y2=y+1;
        }
        cout<<"3 "<<x<<' '<<y<<' '<<x2<<' '<<y2<<'\n';
    }
}
void fudgeValues()
{
    RES0=0;
    if (avgdiff>150) //rand
    {
        RES=900;
        STEPS=4;
        POW1=7;
        PROB1=1;

        if (k==2 && startPoint==0) //bw
        {
            doAll=1;
            STEPS=7;
            RES=1<<STEPS;
            POW1=0;
        }
        else if (k==2 || (k==3 && startPoint==0)) //3 colours
        {
            doAll=2;
            RES=729;
            STEPS=6;
        }
    }
    else if (k==2 && startPoint==0)
    {
        //bw
        doAll=1;
        STEPS=9;
        RES=1<<STEPS;
        PROB1=1;
    }
    /*else if (k==2 || (k==3 && startPoint==0)) //3 colours
    {
        doAll=2;
        RES=729;
        STEPS=6;
    }*/
    else if (avgdiff>32.5 && avgdiff<42.5 && zerodiffs>70000 && zerodiffs<80000)
    {
        //2rand
        RES=900;
        STEPS=8;
        POW1=10;
        POW2=10;
        PROB1=100;
        if (k==2)
        {
            startPoint=0;
            doAll=1;
            RES0=1<<STEPS;
        }
    }
    else if (avgdiff>65 && avgdiff<75 && zerodiffs>70000 && zerodiffs<80000)
    {
        //2rand rgb
        RES=900;
        STEPS=8;
        POW1=3.5;
        POW2=7;
        PROB1=20;
    }
    else //if (k==5 && startPoint==0)
    {
        //rgb
        RES=900;
        STEPS=5;
        POW1=8;
        POW2=8;
        PROB1=100;

        if (k==2)
        {
            doAll=2;
            RES0=243;
        }
    }
}
void isWhite()
{
    startPoint=-1;
    for (int i=0; i<k; ++i)
    {
        if (p[i].r==255 && p[i].g==255 && p[i].b==255)
        {
            wi=i;
            startPoint=0;
            break;
        }
    }
    if (startPoint==-1)
    {
        wi=0;
        for (int i=0; i<k; ++i)
        {
            if (p[i].r+p[i].g+p[i].b>p[wi].r+p[wi].g+p[wi].b)
            {
                wi=i;
            }
        }
    }
}
void solve()
{
    isWhite();
    findUnits1();
    cerr<<"average difference: "<<avgdiff<<endl;
    cerr<<"zero differences: "<<zerodiffs<<endl;
    fudgeValues();
    Color curr= {255,255,255};
    if (RES0)
    {
        for (int i=0; i<RES0; ++i)
        {
            curr=getResulting(i,curr);
        }
        doAll=0;
        startPoint=-1;
    }
    for (int i=RES0; i<RES; ++i)
    {
        curr=getResulting(i,curr);
    }
    cerr<<"operations: "<<OPERATIONS<<endl;
    findUnits2();
    //sort(r,r+RES);
    for (int i=0; i<units.size(); ++i)
    {
        findGroup(i);
    }
    findWhite();
    if (startPoint==0)
    {
        swap(ox,wx);
        swap(oy,wy);
    }
    findOperations();
    cerr<<last.r<<" "<<last.g<<" "<<last.b<<endl;
    //cerr<<wx<<" "<<wy<<endl;
    //cerr<<ox<<" "<<oy<<endl;
    cout<<OPERATIONS<<endl;
    cerr<<"operations: "<<OPERATIONS<<endl;
    cerr<<"units: "<<units.size()<<endl;
    for (int i=0; i<=lastGroup; ++i)
    {
        outputGroup(i);
    }
    int cnt=0;
    //sort(r,r+RES);
    for (int i=0; i<RES; ++i)
    {
        if (inGroup[i].empty()) continue;
        ++cnt;
        //cerr<<r[i].c.r<<" "<<r[i].c.g<<" "<<r[i].c.b<<endl;
    }
    cerr<<cnt<<endl;
}

high_resolution_clock::time_point startT,currT;
int main()
{
    startT=high_resolution_clock::now();
    currT=high_resolution_clock::now();
    srand(0);

    input();
    solve();

    currT=high_resolution_clock::now();

    cerr<<"\n TOTAL TIME: "<<duration_cast<duration<double>>(currT-startT).count()<<endl;

    return 0;
}
