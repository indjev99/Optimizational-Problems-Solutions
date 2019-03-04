#include <iostream>
#include <vector>
#include <unordered_map>
#include <stack>
#include <queue>
#include <unordered_set>
#include <algorithm>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
#include <time.h>
using namespace std;

string name1("a_example.txt");
string name2("b_lovely_landscapes.txt");
string name3("c_memorable_moments.txt");
string name4("d_pet_pictures.txt");
string name5("e_shiny_selfies.txt");

string name=name5;

ifstream inputF(name.c_str());
ofstream outputF(("out3_"+name).c_str());

#define cin inputF
#define cout outputF

unordered_map<string, int> tagMap;
int currTag=0;

struct Image
{
    int num;
    int num2;
    vector<int> tags;
};
const int MAX_N=1e4+10;
vector<Image> hor;
stack<Image> ver;
vector<pair<int,int>> path;
int currScore;
vector<pair<int,int>> maxPath;
int maxScore=-1;
void input()
{
    int n;
    cin>>n;
    char c;
    string s;
    int tgs;
    Image curr;
    for (int i=0;i<n;++i)
    {
        curr.num=i;
        curr.num2=-1;
        cin>>c;
        cin>>tgs;
        curr.tags.resize(tgs);
        for (int j=0;j<tgs;++j)
        {
            cin>>s;
            if (tagMap.find(s)==tagMap.end())
            {
                tagMap[s]=currTag;
                ++currTag;
            }
            curr.tags[j]=tagMap[s];
        }
        sort(curr.tags.begin(),curr.tags.end());
        if (c=='V') ver.push(curr);
        else hor.push_back(curr);
    }
    cerr<<"n: "<<n<<" ver: "<<ver.size()<<" hor: "<<hor.size()<<" tags: "<<currTag<<endl;
}
void output()
{
    cout<<maxPath.size()<<endl;
    for (int i=0;i<maxPath.size();++i)
    {
        cout<<maxPath[i].first;
        if (maxPath[i].second>=0) cout<<' '<<maxPath[i].second;
        cout<<'\n';
    }
}
vector<int> combineTags(const vector<int>& a, const vector<int>& b)
{
    vector<int> c;
    int i=0;
    int j=0;
    while (i<a.size() || j<b.size())
    {
        if (i==a.size()) c.push_back(b[j++]);
        else if (j==b.size()) c.push_back(a[i++]);
        else if (a[i]>b[j]) c.push_back(b[j++]);
        else if (a[i]<b[j]) c.push_back(a[i++]);
        else
        {
            c.push_back(a[i]);
            ++i;
            ++j;
        }
    }
    return c;
}
int commonTags(const vector<int>& a, const vector<int>& b)
{
    int cnt=0;
    int i=0;
    int j=0;
    while (i<a.size() || j<b.size())
    {
        if (i==a.size()) ++j;
        else if (j==b.size()) ++i;
        else if (a[i]>b[j]) ++j;
        else if (a[i]<b[j]) ++i;
        else
        {
            ++cnt;
            ++i;
            ++j;
        }
    }
    return cnt;
}
void match_vert()
{
    Image curr1,curr2,curr3,curr4,curr5,curr6;
    int totalTags=0;
    int currTags;
    int c12,c13,c14,c15,c16,c23,c24,c25,c26,c34,c35,c36,c45,c46,c56;
    while (ver.size()>=2)
    {
        curr1=ver.top();
        ver.pop();
        curr2=ver.top();
        ver.pop();
        if (ver.size()>=4)
        {
            curr3=ver.top();
            ver.pop();
            curr4=ver.top();
            ver.pop();
            curr5=ver.top();
            ver.pop();
            curr6=ver.top();
            ver.pop();
            c12=commonTags(curr1.tags,curr2.tags);
            c13=commonTags(curr1.tags,curr3.tags);
            c14=commonTags(curr1.tags,curr4.tags);
            c15=commonTags(curr1.tags,curr5.tags);
            c16=commonTags(curr1.tags,curr6.tags);
            c23=commonTags(curr2.tags,curr3.tags);
            c24=commonTags(curr2.tags,curr4.tags);
            c25=commonTags(curr2.tags,curr5.tags);
            c26=commonTags(curr2.tags,curr6.tags);
            c34=commonTags(curr3.tags,curr4.tags);
            c35=commonTags(curr3.tags,curr5.tags);
            c36=commonTags(curr3.tags,curr6.tags);
            c45=commonTags(curr4.tags,curr5.tags);
            c46=commonTags(curr4.tags,curr6.tags);
            c56=commonTags(curr5.tags,curr6.tags);
            int minn=min(min(min(min(min(min(min(min(min(min(min(min(min(min(c12,c13),c14),c15),c16),c23),c24),c25),c26),c34),c35),c36),c45),c46),c56);
            if (c13==minn)
            {
                swap(curr2,curr3);
            }
            else if (c14==minn)
            {
                swap(curr2,curr4);
            }
            else if (c15==minn)
            {
                swap(curr2,curr5);
            }
            else if (c16==minn)
            {
                swap(curr2,curr6);
            }
            else if (c23==minn)
            {
                swap(curr1,curr3);
            }
            else if (c24==minn)
            {
                swap(curr1,curr4);
            }
            else if (c25==minn)
            {
                swap(curr1,curr5);
            }
            else if (c26==minn)
            {
                swap(curr1,curr6);
            }
            else if (c34==minn)
            {
                swap(curr1,curr3);
                swap(curr2,curr4);
            }
            else if (c35==minn)
            {
                swap(curr1,curr3);
                swap(curr2,curr5);
            }
            else if (c36==minn)
            {
                swap(curr1,curr3);
                swap(curr2,curr6);
            }
            else if (c45==minn)
            {
                swap(curr1,curr4);
                swap(curr2,curr5);
            }
            else if (c46==minn)
            {
                swap(curr1,curr4);
                swap(curr2,curr6);
            }
            else if (c56==minn)
            {
                swap(curr1,curr5);
                swap(curr2,curr6);
            }
            ver.push(curr3);
            ver.push(curr4);
            ver.push(curr5);
            ver.push(curr6);
        }
        else if (ver.size()>=2)
        {
            curr3=ver.top();
            ver.pop();
            curr4=ver.top();
            ver.pop();
            c12=commonTags(curr1.tags,curr2.tags);
            c13=commonTags(curr1.tags,curr3.tags);
            c14=commonTags(curr1.tags,curr4.tags);
            c23=commonTags(curr2.tags,curr3.tags);
            c24=commonTags(curr2.tags,curr4.tags);
            c34=commonTags(curr3.tags,curr4.tags);
            if (c13<c12 && c13<c14 && c13<c23 && c13<c24 && c13<c34)
            {
                swap(curr2,curr3);
            }
            else if (c14<c12 && c14<c23 && c14<c24 && c14<c34)
            {
                swap(curr2,curr4);
            }
            else if (c23<c12 && c23<c24 && c23<c34)
            {
                swap(curr1,curr3);
            }
            else if (c24<c12 && c24<c34)
            {
                swap(curr1,curr4);
            }
            else if (c34<c12)
            {
                swap(curr1,curr3);
                swap(curr2,curr4);
            }
            ver.push(curr3);
            ver.push(curr4);
        }
        if (!ver.empty())
        {
            curr3=ver.top();
            ver.pop();
            c12=commonTags(curr1.tags,curr2.tags);
            c13=commonTags(curr1.tags,curr3.tags);
            c23=commonTags(curr2.tags,curr3.tags);
            if (c13<c12 && c13<c23)
            {
                swap(curr2,curr3);
            }
            else if (c23<c12)
            {
                swap(curr1,curr3);
            }
            ver.push(curr3);
        }
        hor.push_back({curr1.num,curr2.num,combineTags(curr1.tags,curr2.tags)});
        currTags=hor[hor.size()-1].tags.size();
        totalTags+=currTags;
    }
    cerr<<"ver tags: "<<totalTags<<endl;
}
int score(const vector<int>& a, const vector<int>& b)
{
    int cnt1=0;
    int cnt2=0;
    int cnt3=0;
    int i=0;
    int j=0;
    while (i<a.size() || j<b.size())
    {
        if (i==a.size())
        {
            ++cnt3;
            ++j;
        }
        else if (j==b.size())
        {
            ++cnt1;
            ++i;
        }
        else if (a[i]>b[j])
        {
            ++cnt3;
            ++j;
        }
        else if (a[i]<b[j])
        {
            ++cnt1;
            ++i;
        }
        else
        {
            ++cnt2;
            ++i;
            ++j;
        }
    }
    return min(min(cnt1,cnt2),cnt3);
}
const int MAX_EXPLR=90000;
const int ATTEMPTS=1;
void gen_path()
{
    path.resize(0);
    currScore=0;
    Image last;
    int mxSc;
    int cSc;
    int mxPs;
    for (int i=0;i<hor.size();++i)
    {
        if (i%100==0) cerr<<i<<endl;
        if (i)
        {
            mxSc=-1;
            mxPs=-1;
            for (int j=0;j<MAX_EXPLR && i+j<hor.size();++j)
            {
                cSc=score(last.tags,hor[i+j].tags);
                if (cSc>mxSc)
                {
                    mxSc=cSc;
                    mxPs=i+j;
                }
            }
            //if (i!=mxPs) cerr<<i<<" "<<mxPs<<endl;
            if (i!=mxPs) swap(hor[i],hor[mxPs]);
            currScore+=score(last.tags,hor[i].tags);
        }
        path.push_back({hor[i].num,hor[i].num2});
        last=hor[i];
    }
}
void solve()
{
    match_vert();
    for (int i=0;i<ATTEMPTS;++i)
    {
        random_shuffle(hor.begin(),hor.end());
        gen_path();
        cerr<<i<<" "<<currScore<<endl;
        if (currScore>maxScore)
        {
            maxScore=currScore;
            maxPath=path;
        }
    }
    cerr<<maxScore<<endl;
}
int main()
{
    srand(time(0));
    cerr<<name<<endl;
    input();
    solve();
    output();
    return 0;
}
