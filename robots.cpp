#include<iostream>
#include<queue>
#include<fstream>
#include<stdlib.h>
#include<time.h>
using namespace std;

#define official
#ifdef official
#define cin inF
#define cout outF
#endif // official

ifstream inF("robots.in");
ofstream outF("robots.out");

const int MAX_N=160;
const int MAX_T=12;
const int MAX_K=22;
int MAX_GOALS=500;
const int MAX_GOALS2=10000;

struct Point
{
    int x,y;
};
struct Level
{
    Point r[MAX_K];
    Point g[MAX_GOALS2];
    Point r2[MAX_K];
    string ansS[MAX_GOALS2];
    int goals;
};

int n,m;
int t,k;
char brd[MAX_N][MAX_N];
Level l[MAX_T];
int ansT;
string ansS;

int dist[MAX_N][MAX_N];
char dirr[MAX_N][MAX_N];

void input()
{
    cin>>n>>m;
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<m; ++j)
        {
            cin>>brd[i][j];
        }
    }
    cin>>t>>k;
    for (int i=0; i<t; ++i)
    {
        for (int j=0; j<k; ++j)
        {
            cin>>l[i].r[j].x>>l[i].r[j].y;
            --l[i].r[j].x;
            --l[i].r[j].y;
        }
    }
}
void output()
{
    cout<<ansT+1<<endl;
    cout<<ansS<<endl;
    cerr<<ansS.size()<<endl;
}
bool v[MAX_N][MAX_N];
void DFS(const Level& l, const Point& curr)
{
    int x,y;
    x=curr.x;
    y=curr.y;
    v[x][y]=1;
    if (x>0 && brd[x-1][y]=='.' && v[x-1][y]==0) DFS(l, {x-1,y});
    if (x<n-1 && brd[x+1][y]=='.' && v[x+1][y]==0) DFS(l, {x+1,y});
    if (y>0 && brd[x][y-1]=='.' && v[x][y-1]==0) DFS(l, {x,y-1});
    if (y<m-1 && brd[x][y+1]=='.' && v[x][y+1]==0) DFS(l, {x,y+1});
}
vector<Point> potGoals;
void selectGoals(Level& l, int fp)
{
    int pgs=potGoals.size();
    int num;
    while (l.goals<MAX_GOALS && l.goals<pgs)
    {
        if (l.goals<fp) num=l.goals;
        else num=rand()%(pgs-l.goals);
        l.g[l.goals]=potGoals[num];
        swap(potGoals[num],potGoals[pgs-l.goals]);
        ++l.goals;
    }
}
void findGoals(Level& l)
{
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<m; ++j)
        {
            v[i][j]=0;
        }
    }
    DFS(l,l.r[0]);
    bool first=1;
    Point cp;
    l.goals=0;
    int x,y;
    bool op1,op2,op3,op4;
    potGoals.resize(0);
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<m; ++j)
        {
            if (v[i][j]==0) continue;
            x=i;
            y=j;
            op1=0;
            op2=0;
            op3=0;
            op4=0;
            if (x>0 && brd[x-1][y]=='.') op1=1;
            if (y>0 && brd[x][y-1]=='.') op2=1;
            if (x<n-1 && brd[x+1][y]=='.') op3=1;
            if (y<m-1 && brd[x][y+1]=='.') op4=1;
            if ((op1==0 && op2==0 && op3==0) || (op2==0 && op3==0 && op4==0) || (op3==0 && op4==0 && op1==0) || (op4==0 && op1==0 && op2==0))
            {
                potGoals.push_back({x,y});
            }
        }
    }
    int pgs;
    pgs=potGoals.size();
    if (pgs<MAX_GOALS)
    {
        for (int i=0; i<n; ++i)
        {
            for (int j=0; j<m; ++j)
            {
                if (v[i][j]==0) continue;
                x=i;
                y=j;
                op1=0;
                op2=0;
                op3=0;
                op4=0;
                if (x>0 && brd[x-1][y]=='.') op1=1;
                if (y>0 && brd[x][y-1]=='.') op2=1;
                if (x<n-1 && brd[x+1][y]=='.') op3=1;
                if (y<m-1 && brd[x][y+1]=='.') op4=1;
                if ((op1==0 && op2==0) || (op2==0 && op3==0) || (op3==0 && op4==0) || (op4==0 && op1==0))
                {
                    potGoals.push_back({x,y});
                }
            }
        }
        selectGoals(l,pgs);
    }
    selectGoals(l,0);
}
void findDistAndDirr(Level& l, int gn)
{
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<m; ++j)
        {
            dist[i][j]=-1;
            dirr[i][j]=' ';
        }
    }
    queue<Point> q;
    int x,y;
    q.push(l.g[gn]);
    dist[l.g[gn].x][l.g[gn].y]=0;
    dirr[l.g[gn].x][l.g[gn].y]='-';
    while (!q.empty())
    {
        x=q.front().x;
        y=q.front().y;
        q.pop();
        if (x>0 && brd[x-1][y]=='.' && dist[x-1][y]==-1)
        {
            dist[x-1][y]=dist[x][y]+1;
            dirr[x-1][y]='S';
            q.push({x-1,y});
        }
        if (x<n-1 && brd[x+1][y]=='.' && dist[x+1][y]==-1)
        {
            dist[x+1][y]=dist[x][y]+1;
            dirr[x+1][y]='N';
            q.push({x+1,y});
        }
        if (y>0 && brd[x][y-1]=='.' && dist[x][y-1]==-1)
        {
            dist[x][y-1]=dist[x][y]+1;
            dirr[x][y-1]='E';
            q.push({x,y-1});
        }
        if (y<m-1 && brd[x][y+1]=='.' && dist[x][y+1]==-1)
        {
            dist[x][y+1]=dist[x][y]+1;
            dirr[x][y+1]='W';
            q.push({x,y+1});
        }
    }

    /*cerr<<endl;
    for (int i=0;i<n;++i)
    {
        for (int j=0;j<m;++j)
        {
            cerr<<dirr[i][j];
        }
        cerr<<endl;
    }
    cerr<<"  "<<l.sumdist[gn]<<endl;*/
}
void performMove(Level& l, char mov)
{
    ansS+=mov;
    int dx,dy;
    if (mov=='N')
    {
        dx=-1;
        dy=0;
    }
    else if (mov=='S')
    {
        dx=1;
        dy=0;
    }
    else if (mov=='W')
    {
        dx=0;
        dy=-1;
    }
    else if (mov=='E')
    {
        dx=0;
        dy=1;
    }
    else
    {
        dx=0;
        dy=0;
    }
    int x,y;
    for (int i=0; i<k; ++i)
    {
        x=l.r2[i].x;
        y=l.r2[i].y;
        if (x+dx>=0 && x+dx<n && y+dy>=0 && y+dy<m && brd[x+dx][y+dy]=='.') l.r2[i]= {x+dx,y+dy};
    }
}
int minMoves=-1;
int moveRobot(Level& l, int gn)
{
    int maxD=-1;
    int maxR;
    int x,y;
    int d;
    for (int i=0; i<k; ++i)
    {
        x=l.r2[i].x;
        y=l.r2[i].y;
        d=dist[x][y];
        if (maxD==-1 || d>maxD)
        {
            maxD=d;
            maxR=i;
        }
    }

    x=l.r2[maxR].x;
    y=l.r2[maxR].y;
    //cerr<<maxR<<endl;
    while (dirr[x][y]!='-' && (minMoves==-1 || ansS.size()<=minMoves) && ansS.size()<=10000)
    {
        /*cerr<<" "<<x<<" "<<y<<endl;
        cerr<<dirr[x][y]<<endl;
        cerr<<dist[x][y]<<endl;*/
        performMove(l,dirr[x][y]);
        x=l.r2[maxR].x;
        y=l.r2[maxR].y;
    }
    return maxD;
}
void solveForLevelGoal(Level& l, int gn)
{
    ansS="";
    for (int i=0; i<k; ++i)
    {
        l.r2[i]=l.r[i];
    }
    int maxD;
    do
    {
        maxD=moveRobot(l,gn);
    }
    while (maxD>0 && (minMoves==-1 || ansS.size()<=minMoves) && ansS.size()<=10000);
    l.ansS[gn]=ansS;
}
void solve()
{
    if (n<=50 && m<=50) MAX_GOALS=6500;
    else if (n<=100 && m<=100) MAX_GOALS=1600;
    else if (k==2) MAX_GOALS=1320;
    else MAX_GOALS=500;

    for (int i=0; i<t; ++i)
    {
        findGoals(l[i]);
        //cerr<<endl;
    }
    int minL, minG;
    for (int i=0; i<t && clock()/CLOCKS_PER_SEC<4.9; ++i)
    {
        for (int j=0; j<l[i].goals && clock()/CLOCKS_PER_SEC<4.9; ++j)
        {
            //cerr<<i<<" "<<j<<endl;
            findDistAndDirr(l[i],j);
            solveForLevelGoal(l[i],j);
            if (minMoves==-1 || minMoves>l[i].ansS[j].size())
            {
                minL=i;
                minG=j;
                minMoves=l[i].ansS[j].size();
            }
        }
    }
    ansT=minL;
    ansS=l[minL].ansS[minG];
}

int main()
{
    srand(0);
    input();
    solve();
    output();
}
