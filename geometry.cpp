#include<iostream>
#include<iomanip>
#include<fstream>
#include<algorithm>
#include<math.h>
#include<chrono>
#include<conio.h> ///TODO remove

using namespace std;
using namespace chrono;

#define official

#ifdef official
#define cin inF
#define cout outF
#endif

ifstream inF("geometry.in");
ofstream outF("geometry.out");

const int FAILS_TILL_STOP=1000;

const double PI=3.14159265359;
const double INF=1e9;
const double EPS=1e-9;
const double ZERO_DIST=1.0;

const double E2=2.71/2;
double E2X[32];
void calc_E2X()
{
    E2X[0]=1;
    for (int i=1; i<32; ++i)
    {
        E2X[i]=E2X[i-1]*E2;
    }
}

const int MAX_X=1000;
const int MAX_Y=700;
int MAX_OUT=20;

struct Point
{
    int x,y;
    bool is_outside()
    {
        return x<0 || x>MAX_X || y<0 || y>MAX_Y;
    }
};
bool operator<(const Point& a, const Point& b) //used for sorting and binary search
{
    return a.x<b.x || (a.x==b.x && a.y<b.y);
}
Point operator+(const Point& a, const Point& b)
{
    return {a.x+b.x,a.y+b.y};
}
Point operator-(const Point& a, const Point& b)
{
    return {a.x-b.x,a.y-b.y};
}
Point operator*(const Point& a, int scaler)
{
    return {a.x*scaler,a.y*scaler};
}
Point ortogonal(const Point& a)
{
    return {a.y,-a.x};
}
int cross(const Point& o, const Point& a, const Point& b) //cross section of vector oa and ob
{
    return (a.x-o.x)*(b.y-o.y)-(a.y-o.y)*(b.x-o.x);
}
double dist_p2l(const Point& o, const Point& a, const Point& b) //distance from point to line
{
    double dx,dy;
    dx=b.x-a.x;
    dy=b.y-a.y;
    return fabs(dy*o.x-dx*o.y+b.x*a.y-b.y*a.x)/sqrt(dx*dx+dy*dy);
}
struct State
{
    vector<Point> p; //Points
    //vector<Point> s; //stack
    int out; //number of Point outside the boundaries
    int ops; //number of operations
    int cost; //total cost of operations
    double dist; //distance
    double res; //result
    double alpha; //angle
    vector<Point> tips;

    //functions used for debugging
    void print_points() const
    {
        cerr<<endl;
        cerr<<p.size()<<" points:"<<endl;
        for (int i=0; i<p.size(); ++i)
        {
            cerr<<i<<": "<<p[i].x<<" "<<p[i].y<<endl;
        }
        cerr<<endl;
    }


    void init()
    {
        out=0;
        ops=0;
        cost=0;
    }
    void sort_points() //called only the first time because the Points are kept in sorted order to reduce complexity to N instead of NlogN
    {
        sort(p.begin(),p.end());
    }
    void find_ch(vector<Point>& ch)
    {
        int k=0;
        int n=p.size();
        ch.resize(2*n);
        for (int i=0; i<n; ++i)
        {
            while (k>=2 && ((p[i].x==ch[k-1].x && p[i].y==ch[k-1].y) || cross(ch[k-2],ch[k-1],p[i])<=0)) --k;
            ch[k++]=p[i];
        }
        int t=k+1;
        for (int i=n-2; i>=0; --i)
        {
            while (k>=t && ((p[i].x==ch[k-1].x && p[i].y==ch[k-1].y) || cross(ch[k-2],ch[k-1],p[i])<=0)) --k;
            ch[k++]=p[i];
        }
        ch.resize(k-1);
    }
    void find_dopp(const vector<Point>& ch, vector<double>& dopp, vector<int>& opp) const
    {
        int n=ch.size();
        //cerr<<"Elements in convex hull: "<<n<<endl;
        for (int i=0; i<n; ++i)
        {
            dopp[i]=0;
        }
        if (n<=2) return;
        double cd; //current distance
        dopp[0]=dist_p2l(ch[2],ch[0],ch[1]);
        opp[0]=2;
        for (int i=3; i<n; ++i)
        {
            cd=dist_p2l(ch[i],ch[0],ch[1]);
            if (cd>dopp[0])
            {
                dopp[0]=cd;
                opp[0]=i;
            }
        }
        int j;
        int k;
        for (int i=1; i<n; ++i)
        {
            j=(i+1)%n;
            k=opp[i-1];
            opp[i]=k;
            dopp[i]=dist_p2l(ch[k],ch[i],ch[j]);
            while (true)
            {
                k=(k+1)%n;
                cd=dist_p2l(ch[k],ch[i],ch[j]);
                if (cd<dopp[i]) break;
                dopp[i]=cd;
                opp[i]=k;
            }
        }
        /*for (int i=0;i<n;++i)
        {
            cerr<<i<<": "<<opp[i]<<" "<<dopp[i]<<endl;
        }*/
    }
    void real_calc_res()
    {
        if (out>MAX_OUT || ops>10000) res=INF;
        else res=25*(dist+0.01*E2X[out]*(dist+1))+cost*(dist+1)/50;
        /*else
        {
            proportion=min(1.0,proportion);
            double res0,res1;
            res0=25*(dist+0.01*E2X[out]*(dist+1))+cost*(dist+1)/50;
            res1=100*dist;
            if (out) res1=INF;
            res=res0*proportion+res1*(1-proportion);
            //cerr<<res0<<" "<<res1<<" "<<proportion<<" "<<res<<'\n';
            //getch();
        }*/
    }
    void calc_res(bool printch=false)
    {
        if (p.size()<2)
        {
            alpha=0.0;
            dist=0.0;
            real_calc_res();
            return;
        }
        vector<Point> ch; //convex hull
        find_ch(ch); //find the convex hull
        int n=ch.size();
        if (printch)
        {
            cerr<<"ch: "<<endl;
            for (int i=0;i<n;++i)
            {
                cerr<<ch[i].x<<" "<<ch[i].y<<endl;
            }
            cerr<<endl;
        }
        vector<double> dopp(n); //distance to opposite Point
        vector<int> opp(n); //opposite Point
        find_dopp(ch,dopp,opp);
        dist=-1;
        int j=0;
        for (int i=0; i<n; ++i)
        {
            if (dist<0 || dopp[i]<dist)
            {
                tips.resize(0);
                tips.push_back(ch[opp[i]]);
                tips.push_back(ch[i]);
                tips.push_back(ch[(i+1)%n]);
                dist=dopp[i];
                j=(i+1)%n;
                if (i!=j)
                {
                    alpha=atan2(ch[j].y-ch[i].y,ch[j].x-ch[i].x);
                    if (alpha<0) alpha+=PI;
                }
                else alpha=0;
            }
        }
        real_calc_res();
    }
    int find_point(const Point& a) const
    {
        int l=0,r=p.size()-1,m,last=0;
        while (l<=r)
        {
            m=(l+r)/2;
            if (p[m]<a)
            {
                l=m+1;
            }
            else
            {
                last=m;
                r=m-1;
            }
        }
        //if (p[last].x!=a.x || p[last].y!=a.y) while(1) {};
        return last;
    }
    void homothety(State& new_st, int coeff, int ax, int ay, int bx, int by) const
    {
        Point a= {ax,ay};
        Point b= {bx,by};
        int ai=find_point(a);
        int bi=find_point(b);

        Point b2=(b-a)*coeff+a;
        bool is_out=b2.is_outside();

        int j=0;
        bool b2in=is_out;
        new_st.out=out+is_out;
        new_st.ops=ops+1;
        new_st.cost=cost+(coeff==2?1:6);
        //new_st.s=s;
        new_st.p.resize(p.size()-is_out);
        for (int i=0; i<p.size(); ++i)
        {
            if (!b2in)
            {
                if (!j)
                {
                    if (b2<p[i])
                    {
                        new_st.p[j++]=b2;
                        b2in=true;
                        --i;
                        continue;
                    }
                }
                else
                {
                    if (!(b2<new_st.p[j-1]) && b2<p[i])
                    {
                        new_st.p[j++]=b2;
                        b2in=true;
                        --i;
                        continue;
                    }
                }
            }
            if (i==bi) continue;
            new_st.p[j++]=p[i];
        }
        if (!b2in)
        {
            new_st.p[j++]=b2;
        }
    }
    void ort_homothety(State& new_st, int ax, int ay, int bx, int by) const
    {
        Point a= {ax,ay};
        Point b= {bx,by};
        int ai=find_point(a);
        int bi=find_point(b);

        Point a2=ortogonal(a-b)*2+b;
        Point a3=ortogonal(a-b)*(-2)+b;
        bool is_out2=a2.is_outside();
        bool is_out3=a3.is_outside();
        if (is_out3)
        {
            swap(is_out2,is_out3);
            swap(a2,a3);
        }
        if (!(is_out2 || is_out3) && a3<a2)
        {
            swap(is_out2,is_out3);
            swap(a2,a3);
        }

        int j=0;
        bool a2in=is_out2;
        bool a3in=is_out3;
        new_st.out=out+is_out2+is_out3;
        new_st.ops=ops+1;
        new_st.cost=cost+3;
        //new_st.s=s;
        new_st.p.resize(p.size()+1-is_out2-is_out3);
        for (int i=0; i<p.size(); ++i)
        {
            if (!a2in)
            {
                if (!j)
                {
                    if (a2<p[i])
                    {
                        new_st.p[j++]=a2;
                        a2in=true;
                        --i;
                        continue;
                    }
                }
                else
                {
                    if (!(a2<new_st.p[j-1]) && a2<p[i])
                    {
                        new_st.p[j++]=a2;
                        a2in=true;
                        --i;
                        continue;
                    }
                }
            }
            else if (!a3in)
            {
                if (!j)
                {
                    if (a3<p[i])
                    {
                        new_st.p[j++]=a3;
                        a3in=true;
                        --i;
                        continue;
                    }
                }
                else
                {
                    if (!(a3<new_st.p[j-1]) && a3<p[i])
                    {
                        new_st.p[j++]=a3;
                        a3in=true;
                        --i;
                        continue;
                    }
                }
            }
            if (i==ai) continue;
            new_st.p[j++]=p[i];
        }
        if (!a2in)
        {
            new_st.p[j++]=a2;
        }
        if (!a3in)
        {
            new_st.p[j++]=a3;
        }

    }
    /*void get_from_stack(State& new_st, int ax, int ay) const
    {
        Point a= {ax,ay};
        int ai=find_point(a);

        Point a2=s[s.size()-1];
        bool is_out=a2.is_outside();

        int j=0;
        bool a2in=is_out;
        new_st.out=out+is_out;
        new_st.ops=ops+1;
        new_st.cost=cost+15;
        //new_st.s=s;
        //new_st.s.resize(s.size()-1);
        new_st.p.resize(p.size()-is_out);
        for (int i=0; i<p.size(); ++i)
        {
            if (!a2in)
            {
                if (!j)
                {
                    if (a2<p[i])
                    {
                        new_st.p[j++]=a2;
                        a2in=true;
                        --i;
                        continue;
                    }
                }
                else
                {
                    if (!(a2<new_st.p[j-1]) && a2<p[i])
                    {
                        new_st.p[j++]=a2;
                        a2in=true;
                        --i;
                        continue;
                    }
                }
            }
            if (i==ai) continue;
            new_st.p[j++]=p[i];
        }
        if (!a2in)
        {
            new_st.p[j++]=a2;
        }
    }*/
};
void input(State& s)
{
    int n;
    int m;
    cin>>n;
    s.p.resize(n);
    for (int i=0; i<n; ++i)
    {
        cin>>s.p[i].x>>s.p[i].y;
    }
    cin>>m;
    /*s.s.resize(m);
    for (int i=0; i<n; ++i)
    {
        cin>>s.s[i].x>>s.s[i].y;
    }*/
}
struct Operation
{
    int num;
    vector<int> args;
};
vector<Operation> operations;
double angle;
void output()
{
    cout<<operations.size()<<'\n';
    for (int i=0; i<operations.size(); ++i)
    {
        cout<<operations[i].num<<' ';
        for (int j=0; j<operations[i].args.size(); ++j)
        {
            cout<<(j?" ":"")<<operations[i].args[j];
        }
        cout<<'\n';
    }
    cout<<round(angle*100000)/100000<<endl;
}
void do_operation(const State& s, const Operation& op, State& s2)
{
    if (op.num==1)
    {
        s.homothety(s2,op.args[0],op.args[1],op.args[2],op.args[3],op.args[4]);
    }
    else if (op.num==2)
    {
        s.ort_homothety(s2,op.args[0],op.args[1],op.args[2],op.args[3]);
    }
    /*else
    {
        s.get_from_stack(s2,op.args[0],op.args[1]);
    }*/
}
struct GradedOperation
{
    Operation op;
    double res;
};
bool operator<(const GradedOperation& a, const GradedOperation& b)
{
    return a.res<b.res;
}
double min_possible(const State& s)
{
    //cerr<<"hello"<<endl;
    //return s.res;
    //cerr<<"here\n";
    double min_res=s.res;
    State s2;
    int x1,y1,x2,y2;
    /*if (!s.p.empty() && !s.s.empty()) //operation 3
    {
        for (int i=0; i<s.p.size(); ++i)
        {
            x1=s.p[i].x;
            y1=s.p[i].y;

            s.get_from_stack(s2,x1,y1);
            s2.calc_res();
            if (s2.res<min_res) min_res=s2.res;
        }
    }*/
    if (s.p.size()>1) //operation 2
    {
        for (int i=0; i<s.p.size(); ++i)
        {
            for (int j=0; j<s.p.size(); ++j)
            {
                x1=s.p[i].x;
                y1=s.p[i].y;
                x2=s.p[j].x;
                y2=s.p[j].y;
                if (x1==x2 && y1==y2) continue;

                s.ort_homothety(s2,x1,y1,x2,y2);
                s2.calc_res();
                if (s2.res<min_res) min_res=s2.res;
            }
        }
    }
    if (s.p.size()>1) //operation 1
    {
        for (int i=0; i<s.p.size(); ++i)
        {
            for (int j=0; j<s.p.size(); ++j)
            {
                x1=s.p[i].x;
                y1=s.p[i].y;
                x2=s.p[j].x;
                y2=s.p[j].y;
                if (x1==x2 && y1==y2) continue;

                s.homothety(s2,2,x1,y1,x2,y2);
                s2.calc_res();
                if (s2.res<min_res) min_res=s2.res;

                s.homothety(s2,-2,x1,y1,x2,y2);
                s2.calc_res();
                if (s2.res<min_res) min_res=s2.res;
            }
        }
    }
    return min_res;
}
vector<GradedOperation> grop;
high_resolution_clock::time_point startT,currT,start2T;
vector<int> shuff;
void shuffle(int n)
{
    shuff.resize(n);
    for (int i=0;i<n;++i)
    {
        shuff[i]=i;
    }
    for (int i=n-1;i>0;--i)
    {
        swap(shuff[i],shuff[rand()%(i+1)]);
    }
}
double min_possible_from_tip(const State& s)
{
	double min_res=s.res;
    State s2;
    int x1,y1,x2,y2;
    if (s.p.size()>1) //operation 1
    {
        for (int i=0; i<s.tips.size(); ++i)
        {
            x2=s.tips[i].x;
            y2=s.tips[i].y;
            for (int j=0; j<s.p.size(); ++j)
            {
                x1=s.p[shuff[j]].x;
                y1=s.p[shuff[j]].y;
                if (x1==x2 && y1==y2) continue;

                s.homothety(s2,-2,x1,y1,x2,y2);
                s2.calc_res();
                if (s2.res<min_res) min_res=s2.res;
            }
        }
    }
    return min_res;
}

bool step(State& s, Operation& op, bool force=false)
{
	///returns true if there is no way to improve
    State s2;
    State min_s;
    min_s.res=INF;
    int x1=0,y1=0,x2=0,y2=0;
    if (s.p.size()>1) //operation 1
    {
        shuffle(s.p.size());
        for (int i=0; i<s.tips.size(); ++i)
        {
            x2=s.tips[i].x;
            y2=s.tips[i].y;
            for (int j=0; j<s.p.size(); ++j)
            {
                x1=s.p[shuff[j]].x;
                y1=s.p[shuff[j]].y;
                if (x1==x2 && y1==y2) continue;

                s.homothety(s2,-2,x1,y1,x2,y2);
                s2.calc_res();
                if (s2.res<min_s.res)
                {
                    min_s=s2;
                    op.num=1;
                    op.args.resize(5);
                    op.args[0]=-2;
                    op.args[1]=x1;
                    op.args[2]=y1;
                    op.args[3]=x2;
                    op.args[4]=y2;
                }

                /*s.homothety(s2,2,x1,y1,x2,y2);
                s2.calc_res();
                if (s2.dist<ZERO_DIST && s2.res<s.res)
                {
					cerr<<"PUSH OUTSIDE"<<endl;
                    min_s=s2;
                    op.num=1;
                    op.args.resize(5);
                    op.args[0]=2;
                    op.args[1]=x1;
                    op.args[2]=y1;
                    op.args[3]=x2;
                    op.args[4]=y2;
                    break;
                }*/ ///TODO Try on faster computer
            }
        }
    }
    //cerr<<"    "<<op.args[3]<<" "<<op.args[4]<<endl;
    //getch();
    if (min_s.res>s.res && !force) return true;
    s=min_s;
    return false;
}

bool step_forward(State& s, Operation& op, bool force=false)
{
    ///returns true if there is no way to improve
    State s2;
    State min_s;
    min_s.res=INF;
    int x1=0,y1=0,x2=0,y2=0;
    double min_res=INF;
    double curr_res;
    if (s.p.size()>1) //operation 1
    {
        shuffle(s.p.size());
        for (int i=0; i<s.tips.size(); ++i)
        {
            x2=s.tips[i].x;
            y2=s.tips[i].y;
            for (int j=0; j<s.p.size(); ++j)
            {
                x1=s.p[shuff[j]].x;
                y1=s.p[shuff[j]].y;
                if (x1==x2 && y1==y2) continue;

                s.homothety(s2,-2,x1,y1,x2,y2);
                s2.calc_res();
                if (s2.res<s.res || force)
                {
                    curr_res=min_possible_from_tip(s2);
                    if (curr_res<min_res)
                    {
                        min_s=s2;
                        min_res=curr_res;
                        op.num=1;
                        op.args.resize(5);
                        op.args[0]=-2;
                        op.args[1]=x1;
                        op.args[2]=y1;
                        op.args[3]=x2;
                        op.args[4]=y2;
                    }
                }

                /*s.homothety(s2,2,x1,y1,x2,y2);
                s2.calc_res();
                if (s2.dist<ZERO_DIST && s2.res<s.res)
                {
					cerr<<"PUSH OUTSIDE"<<endl;
                    min_s=s2;
                    op.num=1;
                    op.args.resize(5);
                    op.args[0]=2;
                    op.args[1]=x1;
                    op.args[2]=y1;
                    op.args[3]=x2;
                    op.args[4]=y2;
                    break;
                }*/ ///TODO Try on faster computer
            }
        }
    }
    //cerr<<"    "<<op.args[3]<<" "<<op.args[4]<<endl;
    //getch();
    if (min_res>s.res && !force) return true;
    s=min_s;
    return false;
}
Operation next_op;
bool do_next_op;
bool step_forward2(State& s, Operation& op, bool force=false)
{
    ///returns true if there is no way to improve
    State s2;
    State s3;
    State min_s;
    min_s.res=INF;
    int x1=0,y1=0,x2=0,y2=0;
    Operation curr_next_op;
    bool curr_do_next_op;
    if (s.p.size()>1) //operation 1
    {
        shuffle(s.p.size());
        for (int i=0; i<s.tips.size(); ++i)
        {
            x2=s.tips[i].x;
            y2=s.tips[i].y;
            for (int j=0; j<s.p.size(); ++j)
            {
                x1=s.p[shuff[j]].x;
                y1=s.p[shuff[j]].y;
                if (x1==x2 && y1==y2) continue;

                s.homothety(s2,-2,x1,y1,x2,y2);
                s2.calc_res();
                curr_do_next_op=!step(s2,curr_next_op,false);

                if (s2.res<min_s.res)
                {
                    min_s=s2;
                    next_op=curr_next_op;
                    do_next_op=curr_do_next_op;
                    op.num=1;
                    op.args.resize(5);
                    op.args[0]=-2;
                    op.args[1]=x1;
                    op.args[2]=y1;
                    op.args[3]=x2;
                    op.args[4]=y2;
                }

                /*s.homothety(s2,2,x1,y1,x2,y2);
                s2.calc_res();
                if (s2.dist<ZERO_DIST && s2.res<s.res)
                {
					cerr<<"PUSH OUTSIDE"<<endl;
                    min_s=s2;
                    op.num=1;
                    op.args.resize(5);
                    op.args[0]=2;
                    op.args[1]=x1;
                    op.args[2]=y1;
                    op.args[3]=x2;
                    op.args[4]=y2;
                    break;
                }*/ ///TODO Try on faster computer
            }
        }
    }
    //cerr<<"    "<<op.args[3]<<" "<<op.args[4]<<endl;
    //getch();
    if (min_s.res>s.res && !force) return true;
    s=min_s;
    return false;
}

bool step_forward12(State& s, Operation& op, bool force=false)
{
    ///returns true if there is no way to improve
    State s2;
    State min_s;
    min_s.res=INF;
    int x1=0,y1=0,x2=0,y2=0;
    double min_res=INF;
    double curr_res;
    if (s.p.size()>1) //operation 1
    {
        shuffle(s.p.size());
        for (int i=0; i<s.tips.size(); ++i)
        {
            x2=s.tips[i].x;
            y2=s.tips[i].y;
            for (int j=0; j<s.p.size(); ++j)
            {
                x1=s.p[shuff[j]].x;
                y1=s.p[shuff[j]].y;
                if (x1==x2 && y1==y2) continue;

                s.homothety(s2,-2,x1,y1,x2,y2);
                s2.calc_res();
                if (s2.res<=min_s.res+EPS || force)
                {
                    curr_res=min_possible_from_tip(s2);
                    if (curr_res<min_res)
                    {
                        min_s=s2;
                        min_res=curr_res;
                        op.num=1;
                        op.args.resize(5);
                        op.args[0]=-2;
                        op.args[1]=x1;
                        op.args[2]=y1;
                        op.args[3]=x2;
                        op.args[4]=y2;
                    }
                }

                /*s.homothety(s2,2,x1,y1,x2,y2);
                s2.calc_res();
                if (s2.dist<ZERO_DIST && s2.res<s.res)
                {
					cerr<<"PUSH OUTSIDE"<<endl;
                    min_s=s2;
                    op.num=1;
                    op.args.resize(5);
                    op.args[0]=2;
                    op.args[1]=x1;
                    op.args[2]=y1;
                    op.args[3]=x2;
                    op.args[4]=y2;
                    break;
                }*/ ///TODO Try on faster computer
            }
        }
    }
    //cerr<<"    "<<op.args[3]<<" "<<op.args[4]<<endl;
    //getch();
    if (min_res>s.res && !force) return true;
    s=min_s;
    return false;
}

bool step_greedy(State& s, Operation& op, bool force=false)
{
    ///returns true if there is no way to improve
    State s2;
    State min_s;
    min_s.res=INF;
    int x1,y1,x2,y2;
    /*if (!s.p.empty() && !s.s.empty()) //operation 3
    {
        for (int i=0; i<s.p.size(); ++i)
        {
            x1=s.p[i].x;
            y1=s.p[i].y;

            s.get_from_stack(s2,x1,y1);
            s2.calc_res();
            if (s2.res<min_s.res)
            {
                min_s=s2;
                op.num=3;
                op.args.resize(2);
                op.args[0]=x1;
                op.args[1]=y1;
            }
        }
    }*/
    if (s.p.size()>1) //operation 2
    {
        for (int i=0; i<s.p.size(); ++i)
        {
            for (int j=0; j<s.p.size(); ++j)
            {
                x1=s.p[i].x;
                y1=s.p[i].y;
                x2=s.p[j].x;
                y2=s.p[j].y;
                if (x1==x2 && y1==y2) continue;

                s.ort_homothety(s2,x1,y1,x2,y2);
                s2.calc_res();
                if (s2.res<min_s.res)
                {
                    min_s=s2;
                    op.num=2;
                    op.args.resize(4);
                    op.args[0]=x1;
                    op.args[1]=y1;
                    op.args[2]=x2;
                    op.args[3]=y2;
                }
            }
        }
    }
    if (s.p.size()>1) //operation 1
    {
        for (int i=0; i<s.p.size(); ++i)
        {
            for (int j=0; j<s.p.size(); ++j)
            {
                x1=s.p[i].x;
                y1=s.p[i].y;
                x2=s.p[j].x;
                y2=s.p[j].y;
                if (x1==x2 && y1==y2) continue;

                s.homothety(s2,2,x1,y1,x2,y2);
                s2.calc_res();
                if (s2.res<min_s.res)
                {
                    min_s=s2;
                    op.num=1;
                    op.args.resize(5);
                    op.args[0]=2;
                    op.args[1]=x1;
                    op.args[2]=y1;
                    op.args[3]=x2;
                    op.args[4]=y2;
                }

                s.homothety(s2,-2,x1,y1,x2,y2);
                s2.calc_res();
                if (s2.res<min_s.res)
                {
                    min_s=s2;
                    op.num=1;
                    op.args.resize(5);
                    op.args[0]=-2;
                    op.args[1]=x1;
                    op.args[2]=y1;
                    op.args[3]=x2;
                    op.args[4]=y2;
                }
            }
        }
    }
    //cerr<<"    "<<op.args[0]<<" "<<op.args[1]<<endl;
    if (min_s.res>s.res && !force) return true;
    s=min_s;
    return false;
}
bool step_greedy_from_tip(State& s, Operation& op, bool force=false)
{
    ///returns true if there is no way to improve
    State s2;
    State min_s;
    min_s.res=INF;
    int x1=0,y1=0,x2=0,y2=0;
    /*if (!s.p.empty() && !s.s.empty()) //operation 3
    {
        for (int i=0; i<s.tips.size(); ++i)
        {
            x1=s.tips[i].x;
            y1=s.tips[i].y;
            s.get_from_stack(s2,x1,y1);
            s2.calc_res();
            if (s2.res<min_s.res)
            {
                min_s=s2;
                op.num=3;
                op.args.resize(2);
                op.args[0]=x1;
                op.args[1]=y1;
            }
        }
        //getch();
    }*/
    if (s.p.size()>1) //operation 2
    {
        for (int i=0; i<s.tips.size(); ++i)
        {
            x1=s.tips[i].x;
            y1=s.tips[i].y;
            for (int j=0; j<s.p.size(); ++j)
            {
                x2=s.p[j].x;
                y2=s.p[j].y;
                if (x1==x2 && y1==y2) continue;

                s.ort_homothety(s2,x1,y1,x2,y2);
                s2.calc_res();
                if (s2.res<min_s.res)
                {
                    min_s=s2;
                    op.num=2;
                    op.args.resize(4);
                    op.args[0]=x1;
                    op.args[1]=y1;
                    op.args[2]=x2;
                    op.args[3]=y2;
                }
            }
        }
    }
    if (s.p.size()>1) //operation 1
    {
        for (int i=0; i<s.tips.size(); ++i)
        {
            x2=s.tips[i].x;
            y2=s.tips[i].y;
            for (int j=0; j<s.p.size(); ++j)
            {
                x1=s.p[j].x;
                y1=s.p[j].y;
                if (x1==x2 && y1==y2) continue;

                s.homothety(s2,2,x1,y1,x2,y2);
                s2.calc_res();
                if (s2.res<min_s.res)
                {
                    min_s=s2;
                    op.num=1;
                    op.args.resize(5);
                    op.args[0]=2;
                    op.args[1]=x1;
                    op.args[2]=y1;
                    op.args[3]=x2;
                    op.args[4]=y2;
                }

                s.homothety(s2,-2,x1,y1,x2,y2);
                s2.calc_res();
                if (s2.res<min_s.res)
                {
                    min_s=s2;
                    op.num=1;
                    op.args.resize(5);
                    op.args[0]=-2;
                    op.args[1]=x1;
                    op.args[2]=y1;
                    op.args[3]=x2;
                    op.args[4]=y2;
                }
            }
        }
    }
    //cerr<<"    "<<op.args[3]<<" "<<op.args[4]<<endl;
    //getch();
    if (min_s.res>s.res && !force) return true;
    s=min_s;
    return false;
}


vector<pair<State, int> > potential_answers;
vector<Operation> curr_operations;
void solve_real(State& s)
{
    Operation curr_op;
    bool stop=false;
    s.calc_res();
    int time_to_stop=FAILS_TILL_STOP;
    cerr<<s.p.size()<<endl;
    bool (*step_use)(State&, Operation&, bool);
    bool have_next_op=false;
    if (s.p.size()<=25)
	{
		step_use=&step_forward;
	}
    else if (s.p.size()<=30)
	{
		step_use=&step_forward2;
		have_next_op=true;
		//step_use=&step;
	}
	else
	{
		step_use=&step;
	}
    while (duration_cast<duration<double>>(currT-startT).count()<4.1 && s.dist>EPS && s.p.size()>2 && time_to_stop>0)
    {
        s.calc_res();
        //cerr<<"score (a/d/r): "<<s.alpha*180/PI<<" "<<s.dist<<" "<<s.res<<'\n';

        stop=(*step_use)(s,curr_op,false);
        if (!stop)
        {
            //cerr<<"improve"<<'\n';
            operations.push_back(curr_op);
            if (have_next_op && do_next_op) operations.push_back(next_op);
            time_to_stop=FAILS_TILL_STOP;
        }
        else
        {
            --time_to_stop;
            //cerr<<"not improve at "<<operations.size()<<'\n';
            potential_answers.push_back({s,operations.size()});
            stop=(*step_use)(s,curr_op,true);
            operations.push_back(curr_op);
            if (have_next_op && do_next_op) operations.push_back(next_op);

        }
        //cerr<<"yoyoyoyoyoyoyoy"<<endl;
        currT=high_resolution_clock::now();
        //cerr<<"Time passed:  "<<duration_cast<duration<double>>(currT-startT).count()<<endl;
        //cerr<<s.p.size()<<endl;
    }

    ///TODO: Multiple tries when time is left expected improvement: ~2p.
    ///TODO: Adjust FAILS_TILL_STOP

    s.calc_res();
    //cerr<<"score (a/d/r): "<<s.alpha*180/PI<<" "<<s.dist<<" "<<s.res<<'\n';
    potential_answers.push_back({s,operations.size()});
    int num_ops;
    double min_res=INF;
    for (int i=0; i<potential_answers.size(); ++i)
    {
        //cerr<<potential_answers[i].first.res<<" at "<<potential_answers[i].second<<endl;
        potential_answers[i].first.calc_res();
        if (potential_answers[i].first.res<min_res)
        {
            min_res=potential_answers[i].first.res;
            angle=potential_answers[i].first.alpha;
            num_ops=potential_answers[i].second;
            s=potential_answers[i].first;
        }
    }
    operations.resize(num_ops);
    potential_answers.resize(0);

    MAX_OUT=20;
    time_to_stop=FAILS_TILL_STOP;
    //step_use=&step_greedy_from_tip;
    while (duration_cast<duration<double>>(currT-startT).count()<4.6 && s.dist>EPS && s.p.size()>2 && time_to_stop>0)
    {
        s.calc_res();
        //cerr<<"score (a/d/r): "<<s.alpha*180/PI<<" "<<s.dist<<" "<<s.res<<'\n';

        stop=(*step_use)(s,curr_op,false);
        if (!stop)
        {
            //cerr<<"improve"<<'\n';
            operations.push_back(curr_op);
            if (have_next_op && do_next_op) operations.push_back(next_op);
            time_to_stop=FAILS_TILL_STOP;
        }
        else
        {
            --time_to_stop;
            //cerr<<"not improve at "<<operations.size()<<'\n';
            potential_answers.push_back({s,operations.size()});
            stop=(*step_use)(s,curr_op,true);
            operations.push_back(curr_op);
            if (have_next_op && do_next_op) operations.push_back(next_op);

        }
        //cerr<<"yoyoyoyoyoyoyoy"<<endl;
        currT=high_resolution_clock::now();
        //cerr<<"Time passed:  "<<duration_cast<duration<double>>(currT-startT).count()<<endl;
        //cerr<<s.p.size()<<endl;
    }
    s.calc_res();
    //cerr<<"score (a/d/r): "<<s.alpha*180/PI<<" "<<s.dist<<" "<<s.res<<'\n';
    potential_answers.push_back({s,operations.size()});
    min_res=INF;
    for (int i=0; i<potential_answers.size(); ++i)
    {
        //cerr<<potential_answers[i].first.res<<" at "<<potential_answers[i].second<<endl;
        potential_answers[i].first.calc_res();
        if (potential_answers[i].first.res<min_res)
        {
            min_res=potential_answers[i].first.res;
            angle=potential_answers[i].first.alpha;
            num_ops=potential_answers[i].second;
            s=potential_answers[i].first;
        }
    }
    operations.resize(num_ops);

    cerr<<endl;
    cerr<<"score (a/d/r): "<<s.alpha*180/PI<<" "<<s.dist<<" "<<s.res<<endl;
    cerr<<"Time passed:  "<<duration_cast<duration<double>>(currT-startT).count()<<endl;
}

void step_outside(State& s, Operation& op)
{
    ///returns true if there is no way to improve
    State s2;
    State min_s;
    min_s.res=INF;
    int x1,y1,x2,y2;
    for (int i=0; i<s.p.size(); ++i)
    {
        for (int j=0; j<s.p.size(); ++j)
        {
            x1=s.p[i].x;
            y1=s.p[i].y;
            x2=s.p[j].x;
            y2=s.p[j].y;
            if (x1==x2 && y1==y2) continue;

            s.homothety(s2,2,x1,y1,x2,y2);
            s2.calc_res();
            if (s2.res<min_s.res)
            {
                min_s=s2;
                op.num=1;
                op.args.resize(5);
                op.args[0]=2;
                op.args[1]=x1;
                op.args[2]=y1;
                op.args[3]=x2;
                op.args[4]=y2;
            }

            s.homothety(s2,-2,x1,y1,x2,y2);
            s2.calc_res();
            if (s2.res<min_s.res)
            {
                min_s=s2;
                op.num=1;
                op.args.resize(5);
                op.args[0]=-2;
                op.args[1]=x1;
                op.args[2]=y1;
                op.args[3]=x2;
                op.args[4]=y2;
            }
        }
    }
    s=min_s;
}

void step_outside_rand(State& s, Operation& op, int x1, int y1)
{
    ///returns true if there is no way to improve
    State s2;
    int dist;
    int x2,y2;
    int max_dist=0;
    int mx2,my2;

    for (int i=0; i<s.p.size(); ++i)
    {
        x2=s.p[i].x;
        y2=s.p[i].y;
        if (x2==x1 && y2==y1)
        {
            continue;
        }
        dist=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1);
        if (dist>max_dist)
        {
            max_dist=dist;
            mx2=x2;
            my2=y2;
        }
    }
    op.num=1;
    op.args.resize(5);
    op.args[0]=2;
    op.args[1]=mx2;
    op.args[2]=my2;
    op.args[3]=x1;
    op.args[4]=y1;
    s.homothety(s2,2,mx2,my2,x1,y1);
    s2.calc_res();
    s=s2;
}

void solve_outside(State& s)
{
    State s2;
    Operation curr_op;
    int pn;
    double min_res=INF;
    s2=s;
    curr_operations.resize(0);
    while (s2.p.size()>2 && s2.dist>EPS)
    {
        step_outside(s2,curr_op);
        curr_operations.push_back(curr_op);
    }
    if (s2.res<min_res)
    {
        cerr<<s2.out<<" "<<s2.ops<<" "<<s2.cost<<'\n';
        cerr<<"score (a/d/r): "<<s2.alpha*180/PI<<" "<<s2.dist<<" "<<s2.res<<'\n';
        min_res=s2.res;
        angle=s2.alpha;
        operations=curr_operations;
    }
    //cerr<<endl;
    while (duration_cast<duration<double>>(currT-startT).count()<4.5)
    {
        s2=s;
        curr_operations.resize(0);
        while (s2.p.size()>2 && s2.dist>EPS)
        {
            pn=rand()%s2.p.size();
            step_outside_rand(s2,curr_op,s2.p[pn].x,s2.p[pn].y);
            curr_operations.push_back(curr_op);
        }
        s2.calc_res();
        if (s2.res<min_res)
        {
            cerr<<s2.out<<" "<<s2.ops<<" "<<s2.cost<<'\n';
            cerr<<"score (a/d/r): "<<s2.alpha*180/PI<<" "<<s2.dist<<" "<<s2.res<<'\n';
            min_res=s2.res;
            angle=s2.alpha;
            operations=curr_operations;
        }
        currT=high_resolution_clock::now();
    }
    cerr<<"MIN RESULT: "<<min_res<<endl;
    cerr<<"Time passed:  "<<duration_cast<duration<double>>(currT-startT).count()<<endl;
}
void solve(State& s)
{
    s.init();
    s.sort_points();
    s.calc_res();
    State s2=s;
    int np=s.p.size();
    if (s.p.size()<=22)
    {
    	//MAX_OUT=1;
        //solve_real(s2);
        solve_outside(s2);
    }
    else
    {
    	MAX_OUT=5;
        solve_real(s2);
    }
    /*for (int i=0;i<operations.size();++i)
    {
        //s.print_points();
        //cerr<<operations[i].num<<" "<<operations[i].args[0]<<" "<<operations[i].args[1]<<" "<<operations[i].args[2]<<" "<<operations[i].args[3]<<" "<<operations[i].args[4]<<endl;
        do_operation(s,operations[i],s2);
        s=s2;
        s.calc_res();
    }
    //s.print_points();
    if (fabs(s.alpha-angle)>1000*EPS)
    {
        cerr<<"CRITICAL ERROR "<<angle<<" "<<s.alpha<<endl;
        cerr<<s.dist<<endl;
        //while(np<35) {}; //25 - 35
    }*/
}
int main()
{
    startT=high_resolution_clock::now();
    currT=high_resolution_clock::now();
    srand(3);

    calc_E2X();
    State s;

    input(s);
    solve(s);
    output();

    return 0;
}
