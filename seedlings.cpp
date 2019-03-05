#include<iostream>
#include<conio.h>
#include<queue>
#include<vector>
using namespace std;

const int MAX_N=64;
int n,m;
bool taken[MAX_N][MAX_N];
int state[MAX_N][MAX_N];
int state2[MAX_N][MAX_N];
int pathTo[MAX_N][MAX_N];
int d;

struct Coord
{
    int w;
    int k;
};
struct Block
{
    int r,o;
    vector<Coord> b;
};
vector<Block> bs= {  {1,0,{{0,0},{1,0},{2,0},{3,0}}},
    {1,1,{{0,0},{0,-1},{0,-2},{0,-3}}},
    {2,0,{{0,0},{1,0},{0,1},{1,1}}},
    {3,0,{{0,0},{1,0},{2,0},{1,1}}},
    {3,1,{{0,0},{0,-1},{0,-2},{1,-1}}},
    {3,2,{{0,0},{-1,0},{-2,0},{-1,-1}}},
    {3,3,{{0,0},{0,1},{0,2},{-1,1}}},
    {4,0,{{0,0},{1,0},{2,0},{0,1}}},
    {4,1,{{0,0},{0,-1},{0,-2},{1,0}}},
    {4,2,{{0,0},{-1,0},{-2,0},{0,-1}}},
    {4,3,{{0,0},{0,1},{0,2},{-1,0}}},
    {5,0,{{0,0},{0,1},{1,1},{2,1}}},
    {5,1,{{0,0},{1,0},{1,-1},{1,-2}}},
    {5,2,{{0,0},{0,-1},{-1,-1},{-2,-1}}},
    {5,3,{{0,0},{-1,0},{-1,1},{-1,2}}},
    {6,0,{{0,0},{1,0},{1,1},{2,1}}},
    {6,1,{{0,0},{0,-1},{1,-1},{1,-2}}},
    {7,0,{{0,0},{0,1},{1,1},{1,2}}},
    {7,1,{{0,0},{1,0},{1,-1},{2,-1}}}
};
void initTest()
{
    bool t[7][7];
    for (int k=0; k<bs.size(); ++k)
    {
        for (int i=0; i<7; ++i)
        {
            for (int j=0; j<7; ++j)
            {
                t[i][j]=0;
            }
        }
        for (int i=0; i<4; ++i)
        {
            t[3+bs[k].b[i].w][3+bs[k].b[i].k]=1;
        }
        cerr<<bs[k].r<<' '<<bs[k].o<<endl;
        for (int i=0; i<7; ++i)
        {
            for (int j=0; j<7; ++j)
            {
                if (t[i][j]==0) cerr<<' ';
                else cerr<<'#';
            }
            cerr<<endl;
        }
    }
}
struct Shelf
{
    int w,k,r,o;
};
vector<Shelf> s;
struct Step
{
    int w,k;
    int seen;
};
vector<Step> steps;
int added[MAX_N][MAX_N];

void input()
{
    cin>>n>>m;
    char c;
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<m; ++j)
        {
            cin>>c;
            if (c=='X') taken[i][j]=1;
            else taken[i][j]=0;
        }
    }
}
void output()
{
    cout<<s.size()<<' '<<d<<'\n';
    for (int i=0; i<s.size(); ++i)
    {
        cout<<s[i].w+1<<' '<<s[i].k+1<<' '<<s[i].r<<' '<<s[i].o<<'\n';
    }
}

void update(int w, int k, int& cseen, int dir)
{
    if (taken[w][k]==1) return;
    if (state[w][k]==1 || state[w][k]==4) return;
    if (state[w][k]==0) ++cseen;
    state[w][k]=2;
    if (w<n-1)
    {
        if (taken[w+1][k]==0 && state[w+1][k]==0)
        {
            ++cseen;
            state[w+1][k]=3;
        }
    }
    if (k<m-1)
    {
        if (taken[w][k+1]==0 && state[w][k+1]==0)
        {
            ++cseen;
            state[w][k+1]=3;
        }
    }
    if (w>0)
    {
        if (taken[w-1][k]==0 && state[w-1][k]==0)
        {
            ++cseen;
            state[w-1][k]=3;
        }
    }
    if (k>0)
    {
        if (taken[w][k-1]==0 && state[w][k-1]==0)
        {
            ++cseen;
            state[w][k-1]=3;
        }
    }
}
int t;
void add(int w, int k)
{
    added[w][k]=t++;
}
void rem(int w, int k)
{
    added[w][k]=0;
}
pair<int,int> getNext()
{
    int mp=-1,mt;
    int w,k;
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<m; ++j)
        {
            if (added[i][j]>0)
            {
                if (mp==-1 || pathTo[i][j]>mp || (pathTo[i][j]==mp && added[i][j]<mt))
                {
                    w=i;
                    k=j;
                    mp=pathTo[i][j];
                    mt=added[i][j];
                }
            }
        }
    }
    if (mp!=-1) return {w,k};
    else return {-1,-1};
}
bool vispt[MAX_N][MAX_N];
vector<Coord> cdrn[MAX_N][MAX_N];
vector<Coord> roots;
void DFS_pathTo(const Coord& a)
{
    pathTo[a.w][a.k]=1;
    for (int i=0;i<cdrn[a.w][a.k].size();++i)
    {
        const Coord& next=cdrn[a.w][a.k][i];
        DFS_pathTo(next);
        pathTo[a.w][a.k]+=pathTo[next.w][next.k];
    }
}
void getPathTo()
{
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<m; ++j)
        {
            pathTo[i][j]=0;
            vispt[i][j]=0;
            cdrn[i][j].resize(0);
        }
    }
    roots.resize(0);
    queue<Coord> q;
    Coord curr;
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<m; ++j)
        {
            if (state[i][j]==1)
            {
                if (i<n-1 && taken[i+1][j]==0 && (state[i+1][j]==0 || state[i+1][j]==2 || state[i+1][j]==3) && vispt[i+1][j]==0)
                {
                    q.push({i+1,j});
                    vispt[i+1][j]=1;
                    roots.push_back({i+1,j});
                }
                if (j<m-1 && taken[i][j+1]==0 && (state[i][j+1]==0 || state[i][j+1]==2 || state[i][j+1]==3) && vispt[i][j+1]==0)
                {
                    q.push({i,j+1});
                    vispt[i][j+1]=1;
                    roots.push_back({i,j+1});
                }
                if (i>0 && taken[i-1][j]==0 && (state[i-1][j]==0 || state[i-1][j]==2 || state[i-1][j]==3) && vispt[i-1][j]==0)
                {
                    q.push({i-1,j});
                    vispt[i-1][j]=1;
                    roots.push_back({i-1,j});
                }
                if (j>0 && taken[i][j-1]==0 && (state[i][j-1]==0 || state[i][j-1]==2 || state[i][j-1]==3) && vispt[i][j-1]==0)
                {
                    q.push({i,j-1});
                    vispt[i][j-1]=1;
                    roots.push_back({i,j-1});
                }
            }
        }
    }
    int i,j;
    while (!q.empty())
    {
        curr=q.front();
        q.pop();
        i=curr.w;
        j=curr.k;
        if (i<n-1 && taken[i+1][j]==0 && (state[i+1][j]==0 || state[i+1][j]==2 || state[i+1][j]==3) && vispt[i+1][j]==0)
        {
            q.push({i+1,j});
            vispt[i+1][j]=1;
            cdrn[i][j].push_back({i+1,j});
        }
        if (j<m-1 && taken[i][j+1]==0 && (state[i][j+1]==0 || state[i][j+1]==2 || state[i][j+1]==3) && vispt[i][j+1]==0)
        {
            q.push({i,j+1});
            vispt[i][j+1]=1;
            cdrn[i][j].push_back({i,j+1});
        }
        if (i>0 && taken[i-1][j]==0 && (state[i-1][j]==0 || state[i-1][j]==2 || state[i-1][j]==3) && vispt[i-1][j]==0)
        {
            q.push({i-1,j});
            vispt[i-1][j]=1;
            cdrn[i][j].push_back({i-1,j});
        }
        if (j>0 && taken[i][j-1]==0 && (state[i][j-1]==0 || state[i][j-1]==2 || state[i][j-1]==3) && vispt[i][j-1]==0)
        {
            q.push({i,j-1});
            vispt[i][j-1]=1;
            cdrn[i][j].push_back({i,j-1});
        }
    }
    for (int i=0;i<roots.size();++i)
    {
        DFS_pathTo(roots[i]);
    }
}
void solve()
{
    t=1;
    d=0;
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<m; ++j)
        {
            state[i][j]=0;
            added[i][j]=0;
            if (taken[i][j]) state2[i][j]=9;
            else state2[i][j]=0;
        }
    }
    s.resize(0);
    steps.resize(0);
    if (taken[0][0]==0) add(0,0);
    pair<int, int> curr;
    int cseen=0;
    int bseen=0;
    int w,k;
    int cnt=0;
    int type;
    double bscore=0;
    double cscore;
    int fll;
    int cw;
    int ck;
    bool fail;
    bool succ;
    while (1)
    {
        ++cnt;
        getPathTo();
        if (cnt%3<2)
        //if (cnt%5!=2 && cnt%5!=4)
        {
            curr=getNext();
            w=curr.first;
            k=curr.second;
            if (w==-1) break;
            rem(w,k);

            if (state[w][k]==2 || state[w][k]==3)
            {
                --cseen;
            }
            state[w][k]=1;
            if (w<n-1) update(w+1,k,cseen,0);
            if (k<m-1) update(w,k+1,cseen,1);
            if (w>0) update(w-1,k,cseen,2);
            if (k>0) update(w,k-1,cseen,3);

            steps.push_back({w,k,cseen});
            if (cseen>bseen) bseen=cseen;

            if (w<n-1 && state[w+1][k]==2) add(w+1,k);
            if (k<m-1 && state[w][k+1]==2) add(w,k+1);
            if (w>0 && state[w-1][k]==2) add(w-1,k);
            if (k>0 && state[w][k-1]==2) add(w,k-1);
        }
        else
        {
            bscore=-50000;
            for (int i=0; i<n; ++i)
            {
                for (int j=0; j<m; ++j)
                {
                    if (state[i][j]==1 || state[i][j]==4 || taken[i][j]==1) continue;
                    for (int t=0; t<bs.size(); ++t)
                    {
                        cscore=0;
                        fail=0;
                        succ=0;
                        fll=0;
                        for (int p=0; p<4; ++p)
                        {
                            cw=i+bs[t].b[p].w;
                            ck=j+bs[t].b[p].k;
                            if (cw<0 || cw>=n || ck<0 || ck>=m || state[cw][ck]==1 || state[cw][ck]==4 || taken[cw][ck]==1)
                            {
                                fail=1;
                                break;
                            }
                            else
                            {
                                cscore-=3*pathTo[cw][ck];
                                if (state[cw][ck]==2) succ=1;
                                else cscore+=15;
                                if (cw==0 || (cw>0 && (state[cw-1][ck]==1 || state[cw-1][ck]==4 || taken[cw-1][ck]==1))) ++fll;
                                if (cw==n-1 || (cw<n-1 && (state[cw+1][ck]==1 || state[cw+1][ck]==4 || taken[cw+1][ck]==1))) ++fll;
                                if (ck==0 || (ck>0 && (state[cw][ck-1]==1 || state[cw][ck-1]==4 || taken[cw][ck-1]==1))) ++fll;
                                if (ck==m-1 || (ck<m-1 && (state[cw][ck-1]==1 || state[cw][ck-1]==4 || taken[cw][ck-1]==1))) ++fll;
                            }
                        }
                        if (fail==0 && succ==1)
                        {
                            if (t==2) cscore+=10*fll;
                            else cscore+=8*fll;
                            if (cscore>bscore)
                            {
                                bscore=cscore;
                                type=t;
                                w=i;
                                k=j;
                            }
                        }
                    }
                }

            }
            if (bscore>-5) //10 - 1/3 3;
            {
                for (int p=0; p<4; ++p)
                {
                    cw=w+bs[type].b[p].w;
                    ck=k+bs[type].b[p].k;
                    if (state[cw][ck]==2 || state[cw][ck]==3) --cseen;
                    state[cw][ck]=4;
                    rem(cw,ck);
                }
                cseen+=6;
                if (cseen>bseen) bseen=cseen;
                steps[steps.size()-1].seen=cseen;
            }
            /*cerr<<bscore<<endl;
            for (int i=0; i<n; ++i)
            {
                for (int j=0; j<m; ++j)
                {
                    if (taken[i][j]==1) cout<<'X';
                    else if (state[i][j]==0 || state[i][j]==3) cout<<' ';
                    else if (state[i][j]==2) cout<<'.';
                    else if (state[i][j]==4) cout<<'-';
                    else if (state[i][j]==1) cout<<'o';
                }
                cerr<<endl;
            }
            cerr<<cseen<<endl<<endl;
            getch();*/
        }

        /*cerr<<endl;
        for (int i=0; i<n; ++i)
        {
            for (int j=0; j<m; ++j)
            {
                if (taken[i][j]==1) cout<<'X';
                else if (state[i][j]==0 || state[i][j]==3) cout<<' ';
                else if (state[i][j]==2) cout<<'.';
                else if (state[i][j]==4) cout<<'-';
                else if (state[i][j]==1) cout<<'o';
            }
            cerr<<endl;
        }
        cerr<<cseen<<endl<<endl;
        getch();*/
    }
    if (bseen==0) return;
    for (int i=0; i<steps.size(); ++i)
    {
        w=steps[i].w;
        k=steps[i].k;
        state2[w][k]=8;
        if (steps[i].seen==bseen) break;
    }
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<m; ++j)
        {
            if (state2[i][j]==0)
            {
                if (i>0 && state2[i-1][j]==8) state2[i][j]=1;
                if (i<n-1 && state2[i+1][j]==8) state2[i][j]=1;
                if (j>0 && state2[i][j-1]==8) state2[i][j]=1;
                if (j<m-1 && state2[i][j+1]==8) state2[i][j]=1;
            }
        }
    }

    bscore=0;
    while (bscore>=0)
    {
        bscore=-1;
        for (int i=0; i<n; ++i)
        {
            for (int j=0; j<m; ++j)
            {
                for (int t=0; t<bs.size(); ++t)
                {
                    cscore=0;
                    fail=0;
                    succ=0;
                    fll=0;
                    for (int p=0; p<4; ++p)
                    {
                        cw=i+bs[t].b[p].w;
                        ck=j+bs[t].b[p].k;
                        if (cw<0 || cw>=n || ck<0 || ck>=m || state2[cw][ck]>=7)
                        {
                            fail=1;
                            break;
                        }
                        else
                        {
                            if (state2[cw][ck]==1) succ=1;
                            else cscore+=15;
                            if (cw==0 || (cw>0 && state2[cw-1][ck]>=7)) ++fll;
                            if (cw==n-1 || (cw<n-1 && state2[cw+1][ck]>=7)) ++fll;
                            if (ck==0 || (ck>0 && state2[cw][ck-1]>=7)) ++fll;
                            if (ck==m-1 || (ck<m-1 && state2[cw][ck+1]>=7)) ++fll;
                        }
                    }
                    if (fail==0 && succ==1)
                    {
                        if (t==2) cscore+=10*fll;
                        else cscore+=8*fll;
                        if (cscore>bscore)
                        {
                            bscore=cscore;
                            type=t;
                            w=i;
                            k=j;
                        }
                    }
                }
            }

        }
        if (bscore>=0)
        {
            for (int p=0; p<4; ++p)
            {
                cw=w+bs[type].b[p].w;
                ck=k+bs[type].b[p].k;
                state2[cw][ck]=7;
            }
            s.push_back({w,k,bs[type].r,bs[type].o});
        }
    }

    d=s.size()*6;

    /*cerr<<endl;
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<m; ++j)
        {
            if (state2[i][j]==0) cout<<' ';
            else if (state2[i][j]==1) cout<<'.';
            else if (state2[i][j]==7) cout<<'-';
            else if (state2[i][j]==8) cout<<'o';
            else if (state2[i][j]==9) cout<<'X';
        }
        cerr<<endl;
    }
    cerr<<bseen<<endl;*/

    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<m; ++j)
        {
            if (state2[i][j]==1)
            {
                s.push_back({i,j,0,0});
                ++d;
            }
        }
    }
}

int main()
{
    cin.tie(0);
    ios::sync_with_stdio(0);

    //initTest();

    int tests;
    cin>>tests;

    for (int i=0; i<tests; ++i)
    {
        input();
        solve();
        output();
    }

    return 0;
}
