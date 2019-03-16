#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<conio.h>
#include<algorithm>
#include<chrono>
#define cin in_file
#define cout out_file
using namespace std;
using namespace chrono;
const double INITIAL_T=25;
const int GEN_SIZE=75;
const int KILLED_PER_GEN=50;
double PROBABILITY_CHANGE=0.96;
const int CROSS_RATE=4;
const int MUTATION_TRIES=20;
const int GOOD_MUATIONS=100;
ifstream in_file("keyboard.in");
ofstream out_file("keyboard.out");
struct point
{
    int x,y;
};
struct mut
{
    int a,b;
};
struct mut_by_c
{
    int n;
    int m[8];
};
point pos[32];
int seq[32];
//int cns[32];
mut muts[256];
mut_by_c muts_by_c[32];
int num_muts=0;
int len;
string s;
int i;
void init()
{
    pos[1]={85, 140};
    pos[2]={123, 140};
    pos[3]={161, 140};
    pos[4]={199, 140};
    pos[5]={237, 140};
    pos[6]={275, 140};
    pos[7]={313, 140};
    pos[8]={351, 140};
    pos[9]={389, 140};
    pos[10]={427, 140};
    pos[11]={95, 179};
    pos[12]={133, 179};
    pos[13]={171, 179};
    pos[14]={209, 179};
    pos[15]={247, 179};
    pos[16]={285, 179};
    pos[17]={324, 179};
    pos[18]={362, 179};
    pos[19]={400, 179};
    pos[20]={122, 219};
    pos[21]={160, 219};
    pos[22]={198, 219};
    pos[23]={236, 219};
    pos[24]={275, 219};
    pos[25]={313, 219};
    pos[26]={351, 219};

    seq[0]=1;
    seq[1]=11;
    seq[2]=20;
    seq[3]=12;
    seq[4]=2;
    seq[5]=3;
    seq[6]=13;
    seq[7]=21;
    seq[8]=22;
    seq[9]=14;
    seq[10]=4;
    seq[11]=5;
    seq[12]=15;
    seq[13]=23;
    seq[14]=24;
    seq[15]=16;
    seq[16]=6;
    seq[17]=7;
    seq[18]=17;
    seq[19]=25;
    seq[20]=26;
    seq[21]=18;
    seq[22]=8;
    seq[23]=9;
    seq[24]=19;
    seq[25]=10;

    muts[num_muts++]={1, 2};
    muts[num_muts++]={1, 11};
    muts[num_muts++]={2, 3};
    muts[num_muts++]={2, 12};
    muts[num_muts++]={2, 11};
    muts[num_muts++]={3, 4};
    muts[num_muts++]={3, 13};
    muts[num_muts++]={3, 12};
    muts[num_muts++]={4, 5};
    muts[num_muts++]={4, 14};
    muts[num_muts++]={4, 13};
    muts[num_muts++]={5, 6};
    muts[num_muts++]={5, 15};
    muts[num_muts++]={5, 14};
    muts[num_muts++]={6, 7};
    muts[num_muts++]={6, 16};
    muts[num_muts++]={6, 15};
    muts[num_muts++]={7, 8};
    muts[num_muts++]={7, 17};
    muts[num_muts++]={7, 16};
    muts[num_muts++]={8, 9};
    muts[num_muts++]={8, 18};
    muts[num_muts++]={8, 17};
    muts[num_muts++]={9, 10};
    muts[num_muts++]={9, 19};
    muts[num_muts++]={9, 18};
    muts[num_muts++]={10, 19};
    muts[num_muts++]={11, 12};
    muts[num_muts++]={11, 20};
    muts[num_muts++]={12, 13};
    muts[num_muts++]={12, 21};
    muts[num_muts++]={12, 20};
    muts[num_muts++]={13, 14};
    muts[num_muts++]={13, 22};
    muts[num_muts++]={13, 21};
    muts[num_muts++]={14, 15};
    muts[num_muts++]={14, 23};
    muts[num_muts++]={14, 22};
    muts[num_muts++]={15, 16};
    muts[num_muts++]={15, 24};
    muts[num_muts++]={15, 23};
    muts[num_muts++]={16, 17};
    muts[num_muts++]={16, 25};
    muts[num_muts++]={16, 24};
    muts[num_muts++]={17, 18};
    muts[num_muts++]={17, 26};
    muts[num_muts++]={17, 25};
    muts[num_muts++]={18, 19};
    muts[num_muts++]={18, 26};
    muts[num_muts++]={20, 21};
    muts[num_muts++]={21, 22};
    muts[num_muts++]={22, 23};
    muts[num_muts++]={23, 24};
    muts[num_muts++]={24, 25};
    muts[num_muts++]={25, 26};

    int a,b;
    for (int i=0;i<num_muts;++i)
    {
        a=muts[i].a;
        b=muts[i].b;
        muts_by_c[a].m[muts_by_c[a].n++]=b;
        muts_by_c[b].m[muts_by_c[b].n++]=a;
    }

    /*int curr=1;
    cns['q'-'a']=curr++;
    cns['w'-'a']=curr++;
    cns['e'-'a']=curr++;
    cns['r'-'a']=curr++;
    cns['t'-'a']=curr++;
    cns['y'-'a']=curr++;
    cns['u'-'a']=curr++;
    cns['i'-'a']=curr++;
    cns['o'-'a']=curr++;
    cns['p'-'a']=curr++;
    cns['a'-'a']=curr++;
    cns['s'-'a']=curr++;
    cns['d'-'a']=curr++;
    cns['f'-'a']=curr++;
    cns['g'-'a']=curr++;
    cns['h'-'a']=curr++;
    cns['j'-'a']=curr++;
    cns['k'-'a']=curr++;
    cns['l'-'a']=curr++;
    cns['z'-'a']=curr++;
    cns['x'-'a']=curr++;
    cns['c'-'a']=curr++;
    cns['v'-'a']=curr++;
    cns['b'-'a']=curr++;
    cns['n'-'a']=curr++;
    cns['m'-'a']=curr++;*/
}
struct sol
{
    int ntoc[32];
    int cton[32];
    int fin,fin2;

    point l,r;
    int score;

    void one(char c)
    {
        point& p=pos[cton[c-'a']];
        if (p.x<=l.x)
        {
            score+=round(sqrt((l.x-p.x)*(l.x-p.x)+(l.y-p.y)*(l.y-p.y)));
            l=p;
        }
        else if (p.x>=r.x)
        {
            score+=round(sqrt((r.x-p.x)*(r.x-p.x)+(r.y-p.y)*(r.y-p.y)));
            r=p;
        }
        else
        {
            int ld,rd;
            ld=(l.x-p.x)*(l.x-p.x)+(l.y-p.y)*(l.y-p.y);
            rd=(r.x-p.x)*(r.x-p.x)+(r.y-p.y)*(r.y-p.y);
            if (rd<ld)
            {
                score+=round(sqrt(rd));
                r=p;
            }
            else
            {
                score+=round(sqrt(ld));
                l=p;
            }
        }
    }
    void calc()
    {
        for (int i=0;i<len;++i)
        {
            one(s[i]);
        }
    }
    static int mutC(int c)
    {
        return muts_by_c[c].m[rand()%muts_by_c[c].n];
    }
    void mutate()
    {
        int a,b;
        a=rand()%27+1;
        if (a<=26)
        {
            b=mutC(a);
            swap(ntoc[a],ntoc[b]);
            cton[ntoc[a]]=a;
            cton[ntoc[b]]=b;
        }
        else
        {
            fin=mutC(fin);
        }
    }
    static double P(int e, int en, double T)
    {
        if (en<e) return 1;
        return exp((e-en)/T);
    }
    bool tryToMutate(double T)
    {
        sol sn=*this;
        sn.mutate();
        sn.init();
        sn.calc();

        double r=rand()*RAND_MAX;
        r+=rand();
        r=int(r)%10001;
        r/=10000;
        if (P(score,sn.score,T)>=r)
        {
            *this=sn;
            return 1;
        }
        return 0;
    }
    void init()
    {
        score=0;
        l=pos[fin];
        fin2=cton[s[0]-'a'];
        r=pos[fin2];
        if (l.x>r.x) swap(l,r);
        while (l.x==r.x)
        {
            fin=mutC(fin);
            l=pos[fin];
            if (l.x>r.x) swap(l,r);
        }
    }
    void randomize()
    {
        for (int i=0;i<26;++i)
        {
            /*cton[i]=cns[i];
            ntoc[cns[i]]=i;*/
            ntoc[i+1]=i;
        }
        for (int i=26;i>1;--i)
        {
            swap(ntoc[i],ntoc[rand()%i+1]);
            cton[ntoc[i]]=i;
        }
        fin=rand()%26+1;
    }
    void combine(const sol &p1, const sol &p2)
    {
        //cerr<<endl;
        const sol *p[2]={&p1, &p2};
        fin=p[rand()%2]->fin;
        //cerr<<p1.fin<<" "<<p2.fin<<": "<<fin<<endl;
        int j[2]={0,0};
        bool f;
        for (int i=0;i<26;++i)
        {
            cton[i]=0;
        }
        for (int i=0;i<26;++i)
        {
            f=rand()%2;
            //cerr<<f<<" ";
            while (j[f]<26 && cton[p[f]->ntoc[seq[j[f]]]]!=0) j[f]++;
            if (j[f]<26)
            {
                ntoc[seq[i]]=p[f]->ntoc[seq[j[f]]];
                cton[ntoc[seq[i]]]=seq[i];
            }
        }
        //cerr<<endl;
        //cerr<<endl;
    }
};
bool operator<(sol sol1, sol sol2)
{
    return sol1.score<sol2.score;
}
void input()
{
    cin>>len;
    cin>>s;
}
sol sols[GEN_SIZE],sols2[GEN_SIZE],winner;
void output()
{
    int &fin=winner.fin,&fin2=winner.fin2;
    if (pos[fin].x>pos[fin2].x) swap(fin,fin2);
    while (pos[fin].x==pos[fin2].x)
    {
        cerr<<"FIN ERROR"<<endl;
        fin=sol::mutC(fin);
        if (pos[fin].x>pos[fin2].x) swap(fin,fin2);
    }
    for (int i=1;i<=26;++i)
    {
        cout<<char(winner.ntoc[i]+'a');
    }
    cout<<endl;
    cout<<fin<<" "<<fin2<<endl;
    cerr<<winner.score<<endl;
}
double index_tree[GEN_SIZE];
bool alive[GEN_SIZE];
double get(int p)
{
    //cerr<<"get "<<p<<endl;
    ++p;
    double sum=0;
    while (p>0)
    {
        sum+=index_tree[p-1];
        p-=p&-p;
    }
    return sum;
}
double getOnly(int p)
{
    return get(p)-get(p-1);
}
void add(int p, double a)
{
    ++p;
    while (p<=GEN_SIZE)
    {
        index_tree[p-1]+=a;
        p+=p&-p;
    }
}
void null()
{
    for (int i=0;i<GEN_SIZE;++i)
    {
        index_tree[i]=0;
    }
}
int find(double value)
{
    //cerr<<"";
    int l=0,r=GEN_SIZE-1,m,last;
    while (l<=r)
    {
        m=(l+r)/2;
        //cerr<<l<<" "<<r<<" "<<get(m)<<endl;
        if (get(m)>value) r=m-1;
        else
        {
            last=m;
            l=m+1;
        }
    }
    return last;
}
high_resolution_clock::time_point start,curr,start2;
void solve()
{
    double T=INITIAL_T;
    double reserve=0;
    double prob;
    int kill;
    int j,k,k2;
    int cc;
    bool fail=1;
    for (int i=0;i<GEN_SIZE;++i)
    {
        sol& curr_sol=sols[i];
        curr_sol.randomize();
        curr_sol.init();
        curr_sol.calc();
    }
    sort(sols,sols+GEN_SIZE);
    curr=high_resolution_clock::now();
    while (duration_cast<duration<double>>(curr-start).count()<2.85-reserve)
    {
        start2=high_resolution_clock::now();
        prob=1000;
        null();
        for (int i=0;i<GEN_SIZE;++i)
        {
            alive[i]=1;
            add(i,1000-prob);
            prob*=PROBABILITY_CHANGE;
        }
        for (int i=0;i<KILLED_PER_GEN;++i)
        {
            prob=get(GEN_SIZE-1);
            kill=rand()*RAND_MAX;
            kill+=rand();
            kill=find(kill%int(prob));
            //getch();
            alive[kill]=0;
            add(kill,-getOnly(kill));
        }
        j=0;
        for (int i=0;i<GEN_SIZE;++i)
        {
            if (alive[i]) sols2[j++]=sols[i];
        }
        j=0;
        for (int i=GEN_SIZE-KILLED_PER_GEN;i<GEN_SIZE;++i)
        {
            sol &curr_sol=sols2[i];
            if (rand()%CROSS_RATE==0)
            {
                k=rand()%(GEN_SIZE-KILLED_PER_GEN);
                k2=rand()%(GEN_SIZE-KILLED_PER_GEN);
                curr_sol.combine(sols2[k],sols2[k2]);
            }
            else
            {
                curr_sol=sols2[j++];
                curr_sol.mutate();
            }
            curr_sol.mutate();
            curr_sol.init();
            curr_sol.calc();
            for (int gm=0;gm<GOOD_MUATIONS;++gm)
            {
                cc=0;
                while ((fail=!curr_sol.tryToMutate(T)) && cc<MUTATION_TRIES) ++cc;
                if (fail) break;
            }
            //curr_sol.tryToMutate(T);
        }
        sort(sols2+GEN_SIZE-KILLED_PER_GEN,sols2+GEN_SIZE);
        j=GEN_SIZE-KILLED_PER_GEN;
        i=0;
        k=0;
        while (i<GEN_SIZE-KILLED_PER_GEN && j<GEN_SIZE)
        {
            if (sols2[i].score<sols2[j].score) sols[k++]=sols2[i++];
            else sols[k++]=sols2[j++];
        }
        if (i==GEN_SIZE-KILLED_PER_GEN)
        {
            for (;j<GEN_SIZE;++j)
            {
                sols[k++]=sols2[j];
            }
        }
        else
        {
            for (;i<GEN_SIZE-KILLED_PER_GEN;++i)
            {
                sols[k++]=sols2[i];
            }
        }
        curr=high_resolution_clock::now();
        reserve=duration_cast<duration<double>>(curr-start2).count();

        T-=T*reserve/(3-duration_cast<duration<double>>(curr-start).count());
        //getch();
    }
    cerr<<T<<endl;
    winner=sols[0];
}
int main()
{
    start=high_resolution_clock::now();
    srand(1);

    init();
    input();
    solve();
    output();

    return 0;
}
