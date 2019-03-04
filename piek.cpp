#include<iostream>
#include<vector>
#include<stdlib.h>
#include<math.h>
#include<chrono>
using namespace std;
using namespace chrono;

const int MAX_N=512;
const int MAX_POP=2000;
int POP;
const double EVAP=0.4;
const double q0=0.35;
const double Q=1000;
const double ELIT=30;
const int res=RAND_MAX+1;
double avg;
const int REPS_ARR[8]= {9,28,9,9,9,9,9,2};

struct Sol
{
    int len;
    int seq[MAX_N];
};

int n;
int d[MAX_N][MAX_N];
double pher[MAX_N][MAX_N];
//double stPher[MAX_N];
Sol ans;
Sol s[MAX_POP];

void input()
{
    cin>>n;
    POP=50;
    avg=0;
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<n; ++j)
        {
            cin>>d[i][j];
            avg+=d[i][j];
        }
    }
    /*n=400;
    for (int i=0;i<n;++i)
    {
        for (int j=0;j<n;++j)
        {
            d[i][j]=rand()+5000;
            d[j][i]=d[i][j];
            d[j][j]=0;
            avg+=d[i][j];
        }
    }*/
    avg=avg/n/n;
}

void output()
{
    cout<<ans.len<<'\n';
    for (int i=0; i<n; ++i)
    {
        cout<<ans.seq[i]+1<<' ';
    }
    cout<<ans.seq[0]+1<<'\n';
}

vector<pair<int,int> > ff;
int nx[500];
void output2()
{
    int counter=0;
    int i;
    for(i=0; i<n-1; i++)
    {
        ff.push_back(make_pair(ans.seq[i],ans.seq[i+1]));
        nx[ans.seq[i]]=ans.seq[i+1];
    }
    ff.push_back(make_pair(ans.seq[n-1],ans.seq[0]));
    nx[ans.seq[n-1]]=ans.seq[0];
    int cnt=0;
    int sz=n;
    int j,x1,x2,y1,y2,pr,k;
    while(cnt<16)
    {
        for(i=0; i<sz; i++)
        {
            for(j=i-1; j>=0; j--)
            {
                if(ff[i].first==ff[j].first||ff[i].second==ff[j].second||ff[i].first==ff[j].second||ff[i].second==ff[j].first)
                    continue;
                if(nx[ff[i].first]!=ff[i].second)
                    swap(ff[i].first,ff[i].second);
                if(nx[ff[j].first]!=ff[j].second)
                    swap(ff[j].first,ff[j].second);
                x1=d[ff[i].first][ff[i].second];
                x2=d[ff[j].first][ff[j].second];
                y1=d[ff[i].first][ff[j].first];
                y2=d[ff[j].second][ff[i].second];
                if(x1+x2>y1+y2)
                {
                    ++counter;
                    //cout<<ff[i].first<<ff[i].second<<" "<<ff[j].first<<ff[j].second<<endl;
                    ans.len-=(x1+x2);
                    ans.len+=(y1+y2);
                    pr=ff[i].second;
                    int ni;
                    for(k=nx[pr]; k!=ff[j].second;)
                    {
                        ni=nx[k];
                        nx[k]=pr;
                        pr=k;
                        k=ni;
                    }
                    nx[ff[i].first]=ff[j].first;
                    nx[ff[i].second]=ff[j].second;
                    swap(ff[i].second,ff[j].first);
                }
            }
        }
        cnt++;
    }
    cerr<<"cnt: "<<counter<<endl;
    cout<<ans.len<<endl;
    pr=0;

    cout<<ff[0].first+1<<" "<<ff[0].second+1<<" ";
    for(i=0; i<n-1; i++)
    {
        for(j=0; j<n; j++)
        {
            if(j==pr)
                continue;
            if(ff[j].second==ff[pr].second)
            {
                swap(ff[j].first,ff[j].second);
            }
            if(ff[j].first==ff[pr].second)
            {
                cout<<ff[j].second+1<<" ";
                pr=j;
                break;
            }
        }
    }
    cout<<endl;
}


bool u[MAX_N];
double prob[MAX_N];
double tp;
int begSeq[MAX_N];
int bsp=n;
void initBegSeq()
{
    bsp=n;
    for (int i=0; i<n; ++i)
    {
        begSeq[i]=i;
    }
}
int getNext()
{
    /*if (bsp==0) bsp=n;
    int np,next;
    np=rand()%bsp;
    --bsp;
    swap(begSeq[np],begSeq[bsp]);
    return begSeq[bsp];*/
    return rand()%n;
}
void genSol(int p)
{
    for (int i=0; i<n; ++i)
    {
        u[i]=0;
    }
    int beg=getNext();
    double P;
    int q;
    double maxx=0;

    /*tp=0;  //bad idea
    for (int i=0; i<n; ++i)
    {
        tp+=stPher[i];
    }
    if (tp>0)
    {
        for (int i=0; i<n; ++i)
        {
            P=stPher[i]/tp;
            if (rand()<=P*res)
            {
                beg=i;
                break;
            }
            tp-=stPher[i];
        }
    }*/

    int cnt=1;
    auto& cs=s[p];
    int curr=beg;
    int next;
    cs.seq[0]=curr;
    cs.len=0;
    while (cnt<n)
    {
        u[curr]=1;
        if (rand()<=q0*res)
        {
            tp=0;
            for (int i=0; i<n; ++i)
            {
                if (u[i]==1) continue;
                next=i;
                prob[i]=pher[curr][i]/d[curr][i];
                tp+=prob[i];
            }
            if (tp>0)
            {
                for (int i=0; i<n; ++i)
                {
                    if (u[i]==1) continue;
                    P=prob[i]/tp;
                    if (rand()<=P*res)
                    {
                        next=i;
                        break;
                    }
                    tp-=prob[i];
                }
            }
        }
        else
        {
            next=-1;
            maxx=0;
            for (int i=0; i<n; ++i)
            {
                if (u[i]==1) continue;
                prob[i]=pher[curr][i]/d[curr][i];
                if (next==-1 || prob[i]>maxx)
                {
                    next=i;
                    maxx=prob[i];
                }
            }
        }
        /*if (next==-1)
        {
            cerr<<"ERROR"<<endl;
            while (1) {};
        }
        else
        {
            cerr<<"yasha"<<endl;
        }*/
        cs.seq[cnt]=next;
        cs.len+=d[curr][next];
        curr=next;
        ++cnt;
    }
    cs.len+=d[curr][beg];
}
void update()
{
    for (int i=0; i<n; ++i)
    {
        //stPher[i]*=(1-EVAP);
        for (int j=0; j<n; ++j)
        {
            pher[i][j]*=(1-EVAP);
        }
    }
    /*int a;
    for (int i=0;i<n;++i)
    {
        for (int j=0;j<n;++j)
        {
            cerr<<pher[i][j]<<" ";
        }
        cerr<<endl;
    }
    cerr<<endl;
    cin>>a;*/
    double add=Q*ELIT/ans.len;
    //stPher[ans.seq[0]]+=add;
    pher[ans.seq[0]][ans.seq[n-1]]+=add;
    pher[ans.seq[n-1]][ans.seq[0]]+=add;
    for (int j=1; j<n; ++j)
    {
        pher[ans.seq[j-1]][ans.seq[j]]+=add;
        pher[ans.seq[j]][ans.seq[j-1]]+=add;
    }
    for (int i=0; i<POP; ++i)
    {
        add=Q/s[i].len;
        //cerr<<"add: "<<add<<endl;
        //stPher[s[i].seq[0]]+=add;
        pher[s[i].seq[0]][s[i].seq[n-1]]+=add;
        pher[s[i].seq[n-1]][s[i].seq[0]]+=add;
        for (int j=1; j<n; ++j)
        {
            pher[s[i].seq[j-1]][s[i].seq[j]]+=add;
            pher[s[i].seq[j]][s[i].seq[j-1]]+=add;
        }
    }
    /*for (int i=0;i<n;++i)
    {
        for (int j=0;j<n;++j)
        {
            cerr<<pher[i][j]<<" ";
        }
        cerr<<endl;
    }
    cerr<<endl;
    cin>>a;*/
}
high_resolution_clock::time_point start,curr,start2;
void swapEdges(Sol& s, int a, int b)
{
    int k=(b-a)/2;
    for (int i=0; i<k; ++i)
    {
        swap(s.seq[a+1+i],s.seq[b-i]);
    }
}
void opt2(Sol& s, int lim)
{
    if (lim==0) return;
    int minchange,change;
    int ni,nj,mini,minj,ri,rj,rni,rnj;
    int cnt=0;
    do
    {
        minchange=0;
        for (int i=0; i<n; ++i)
        {
            for (int j=i+2; j<n; ++j)
            {
                ni=i+1;
                if (ni==n) ni=0;
                nj=j+1;
                if (nj==n) nj=0;
                ri=s.seq[i];
                rj=s.seq[j];
                rni=s.seq[ni];
                rnj=s.seq[nj];
                change=d[ri][rj]+d[rni][rnj]-d[ri][rni]-d[rj][rnj];
                if (change<minchange)
                {
                    minchange=change;
                    mini=i;
                    minj=j;
                }
            }
        }
        if (minchange<0) swapEdges(s,mini,minj);
        s.len+=minchange;
        ++cnt;
    }
    while (minchange<0 && cnt<lim);
}
void opt2(Sol& s)
{
    int minchange,change;
    int ni,nj,mini,minj,ri,rj,rni,rnj;
    int cnt=0;
    do
    {
        minchange=0;
        for (int i=0; i<n; ++i)
        {
            for (int j=i+2; j<n; ++j)
            {
                ni=i+1;
                if (ni==n) ni=0;
                nj=j+1;
                if (nj==n) nj=0;
                ri=s.seq[i];
                rj=s.seq[j];
                rni=s.seq[ni];
                rnj=s.seq[nj];
                change=d[ri][rj]+d[rni][rnj]-d[ri][rni]-d[rj][rnj];
                if (change<minchange)
                {
                    minchange=change;
                    mini=i;
                    minj=j;
                }
            }
        }
        if (minchange<0) swapEdges(s,mini,minj);
        s.len+=minchange;
        ++cnt;
        curr=high_resolution_clock::now();
    }
    while (minchange<0 && duration_cast<duration<double>>(curr-start).count()<0.98);
}
int arr[MAX_N];
void trapEdges(Sol& s, int a, int b, int c, int mode)
{
    if (mode<0)
    {
        if (mode==-1)
        {
            swapEdges(s,b,c);
        }
        else if (mode==-2)
        {
            swapEdges(s,a,c);
        }
        else if (mode==-3)
        {
            swapEdges(s,b,c);
        }
        return;
    }

    int curr=0;
    int k;
    for (int i=0; i<n; ++i)
    {
        arr[i]=s.seq[i];
    }
    for (int i=b+1; i<=c; ++i)
    {
        s.seq[curr++]=arr[i];
    }
    for (int i=a+1; i<=b; ++i)
    {
        s.seq[curr++]=arr[i];
    }
    for (int i=c+1; i<n; ++i)
    {
        s.seq[curr++]=arr[i];
    }
    for (int i=0; i<=a; ++i)
    {
        s.seq[curr++]=arr[i];
    }
    int d1=c-b;
    int d2=b-a;
    int d3=n-d1-d2;
    if (mode==1)
    {
        k=d1/2;
        for (int i=0; i<k; ++i)
        {
            swap(s.seq[i],s.seq[d1-1-i]);
        }
    }
    else if (mode==2)
    {
        k=d2/2;
        for (int i=0; i<k; ++i)
        {
            swap(s.seq[d1+i],s.seq[d1+d2-1-i]);
        }
    }
    else if (mode==3)
    {
        k=d3/2;
        for (int i=0; i<k; ++i)
        {
            swap(s.seq[d1+d2+i],s.seq[n-1-i]);
        }
    }
}
void opt3(Sol& s)
{
    int minchange,change;
    int ni,nj,nk,mini,minj,mink,minmode,ri,rj,rk,rni,rnj,rnk;
    do
    {
        for (int i=0; i<n && duration_cast<duration<double>>(curr-start).count()<0.98; ++i)
        {
            for (int j=i+2; j<n && duration_cast<duration<double>>(curr-start).count()<0.98; ++j)
            {
                minchange=0;
                for (int k=j+2; k<n; ++k)
                {
                    ni=i+1;
                    if (ni==n) ni=0;
                    nj=j+1;
                    if (nj==n) nj=0;
                    nk=k+1;
                    if (nk==n) nk=0;
                    ri=s.seq[i];
                    rj=s.seq[j];
                    rk=s.seq[k];
                    rni=s.seq[ni];
                    rnj=s.seq[nj];
                    rnk=s.seq[nk];
                    change=d[ri][rnj]+d[rk][rni]+d[rj][rnk]-d[ri][rni]-d[rj][rnj]-d[rk][rnk];
                    if (change<minchange)
                    {
                        minchange=change;
                        mini=i;
                        minj=j;
                        mink=k;
                        minmode=0;
                    }
                    change=d[ri][rk]+d[rnj][rni]+d[rj][rnk]-d[ri][rni]-d[rj][rnj]-d[rk][rnk];
                    if (change<minchange)
                    {
                        minchange=change;
                        mini=i;
                        minj=j;
                        mink=k;
                        minmode=1;
                    }
                    change=d[ri][rnj]+d[rk][rj]+d[rni][rnk]-d[ri][rni]-d[rj][rnj]-d[rk][rnk];
                    if (change<minchange)
                    {
                        minchange=change;
                        mini=i;
                        minj=j;
                        mink=k;
                        minmode=2;
                    }
                    change=d[rnk][rnj]+d[rk][rni]+d[rj][ri]-d[ri][rni]-d[rj][rnj]-d[rk][rnk];
                    if (change<minchange)
                    {
                        minchange=change;
                        mini=i;
                        minj=j;
                        mink=k;
                        minmode=3;
                    }
                    change=d[ri][rni]+d[rj][rk]+d[rnj][rnk]-d[ri][rni]-d[rj][rnj]-d[rk][rnk];
                    if (change<minchange)
                    {
                        minchange=change;
                        mini=i;
                        minj=j;
                        mink=k;
                        minmode=-1;
                    }
                    change=d[rj][rnj]+d[ri][rk]+d[rni][rnk]-d[ri][rni]-d[rj][rnj]-d[rk][rnk];
                    if (change<minchange)
                    {
                        minchange=change;
                        mini=i;
                        minj=j;
                        mink=k;
                        minmode=-2;
                    }
                    change=d[rk][rnk]+d[ri][rj]+d[rni][rnj]-d[ri][rni]-d[rj][rnj]-d[rk][rnk];
                    if (change<minchange)
                    {
                        minchange=change;
                        mini=i;
                        minj=j;
                        mink=k;
                        minmode=-3;
                    }
                }
                if (minchange<0)
                {
                    //cerr<<"mc3: "<<-minchange<<" "<<minmode<<endl;
                    trapEdges(s,mini,minj,mink,minmode);
                    s.len+=minchange;
                }
                curr=high_resolution_clock::now();
            }
            curr=high_resolution_clock::now();
        }
        curr=high_resolution_clock::now();
    }
    while (duration_cast<duration<double>>(curr-start).count()<0.98);
}
void solve()
{
    initBegSeq();
    for (int i=0; i<n; ++i)
    {
        //stPher[i]=Q*POP/n/avg/EVAP/2;
        for (int j=0; j<n; ++j)
        {
            pher[i][j]=Q*POP/n/avg/EVAP/2;
        }
    }

    bool first=1;

    curr=high_resolution_clock::now();
    int reps=REPS_ARR[(n-1)/50];
    while (duration_cast<duration<double>>(curr-start).count()<0.7)
    {
        for (int i=0; i<POP; ++i)
        {
            genSol(i);
            opt2(s[i],reps);
            if (s[i].len<ans.len || first)
            {
                first=0;
                ans=s[i];
            }
        }
        update();
        curr=high_resolution_clock::now();
    }
    /*for (int i=0;i<n;++i)
    {
        for (int j=0;j<n;++j)
        {
            cerr<<pher[i][j]<<" ";
        }
        cerr<<endl;
    }
    cerr<<endl;*/
    opt3(ans);
}
int main()
{
    start=high_resolution_clock::now();
    srand(1);

    ios_base::sync_with_stdio(0);
    cin.tie(0);

    input();
    //start=high_resolution_clock::now();
    solve();
    output();

    return 0;
}
