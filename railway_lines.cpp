#include <bits/stdc++.h>
using namespace std;

const int MAX_N = 10000;

const long long PRECISION = 1000;

const int RAND_FACT = PRECISION * 61000;

int n, m;
long long l, k;

struct Point
{
    long long x, y;
};

struct Edge
{
    int a, b;
    vector<int> inters;

    long long distCost;
    long long currCost;

    int lastAnsT;
    long long probIn;
};

Point ints[MAX_N];
Edge edges[MAX_N];

long long bestScore;
vector<int> bestAns, ans;

void input()
{
    cin >> n >> m >> k >> l;
    for (int i = 0; i < n; ++i)
    {
        cin >> ints[i].x >> ints[i].y;
    }
    for (int i = 0; i < m; ++i)
    {
        cin >> edges[i].a >> edges[i].b;
        --edges[i].a;
        --edges[i].b;
    }
}

void output()
{
    for (int i = 0; i < n - 1; ++i)
    {
        cout << bestAns[i] + 1 << "\n";
    }
}

int currAnsT;
int currT;
int setT[MAX_N];
int par[MAX_N];
int cmpSz[MAX_N];

int findRootUtil(int a)
{
    if (a == par[a]) return a;
    return par[a] = findRootUtil(par[a]);
}

int findRoot(int a)
{
    if (setT[a] < currT)
    {
        cmpSz[a] = 1;
        setT[a] = currT;
        return par[a] = a;
    }
    return par[a] = findRootUtil(par[a]);
}

bool unite(int a, int b)
{
    a = findRoot(a);
    b = findRoot(b);
    if (a == b) return false;
    if (cmpSz[a] > cmpSz[b])
    {
        par[b] = a;
        cmpSz[a] += cmpSz[b];
    }
    else
    {
        par[a] = b;
        cmpSz[b] += cmpSz[a];
    }
    return true;
}

long long distCost(const Point& a, const Point& b)
{
    long long dstSq = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
    return round(PRECISION * k * sqrt(dstSq));
}

bool ccw(const Point& a, const Point& b, const Point& c)
{
    return (b.y - a.y) * (c.x - b.x) > (b.x - a.x) * (c.y - b.y);
}

bool intersect(const Point& a, const Point& b, const Point& c, const Point& d)
{
    return ccw(a, b, c) != ccw(a, b, d) && ccw(c, d, a) != ccw(c, d, b);
}

priority_queue<pair<long long, int>> pq;

long long randomNumber()
{
    long long rn = rand();
    rn *= RAND_MAX;
    rn += rand();
    rn *= RAND_MAX;
    rn += rand();
    return rn;
}

void init()
{
    ++currT;
    ++currAnsT;
    while (!pq.empty()) pq.pop();
    ans.clear();

    bool doRand = rand() % 2;
    bool doPredict = rand() % 2;

    for (int i = 0; i < m; ++i)
    {
        edges[i].currCost = edges[i].distCost;
        if (doRand) edges[i].currCost += randomNumber() % RAND_FACT;
        edges[i].probIn = PRECISION * (edges[i].lastAnsT == currAnsT - 1);
    }

    for (int i = 0; i < m; ++i)
    {
        for (int other : edges[i].inters)
        {
            edges[other].currCost += edges[i].probIn * l;
        }
    }

    for (int i = 0; i < m; ++i)
    {
        pq.push({-edges[i].currCost, i});
    }
}

void minSpanning()
{
    int cnt = 0;
    while (cnt < n - 1)
    {
        pair<long long, int> currPair = pq.top();
        int edgeInd = currPair.second;
        pq.pop();
        if (edges[edgeInd].currCost > -currPair.first)
        {
            pq.push({-edges[edgeInd].currCost, edgeInd});
            continue;
        }
        if (unite(edges[edgeInd].a, edges[edgeInd].b))
        {
            for (int other : edges[edgeInd].inters)
            {
                edges[other].currCost += PRECISION * l;
            }
            ans.push_back(edgeInd);
            ++cnt;
        }
    }
}

void evalAns()
{
    long long score = 0;
    for (int e : ans)
    {
        score += edges[e].distCost;
        edges[e].lastAnsT = currAnsT;
    }
    for (int e : ans)
    {
        for (int other : edges[e].inters)
        {
            if (other < e) continue;
            if (edges[other].lastAnsT == currAnsT) score += l * PRECISION;
        }
    }
    if (bestScore == -1 || score < bestScore)
    {
        bestScore = score;
        bestAns = ans;
    }
}

void localOpt()
{
    for (int j = 0; j < m; ++j)
    {
        edges[j].currCost = edges[j].distCost;
    }

    for (int e : ans)
    {
        for (int other : edges[e].inters)
        {
            edges[other].currCost += l * PRECISION;
        }
    }

    for (int i = 0; i < n - 1; ++i)
    {
        ++currT;
        for (int j = 0; j < n - 1; ++j)
        {
            if (j == i) continue;
            unite(edges[ans[j]].a, edges[ans[j]].b);
        }

        int bestEdge = ans[i];
        int bestScore = edges[bestEdge].currCost;

        for (int other : edges[bestEdge].inters)
        {
            edges[other].currCost -= l * PRECISION;
        }

        for (int j = 0; j < m; ++j)
        {
            if (findRoot(edges[j].a) != findRoot(edges[j].b))
            {
                int score = edges[j].currCost;
                if (score < bestScore)
                {
                    bestScore = score;
                    bestEdge = j;
                }
            }
        }

        for (int other : edges[bestEdge].inters)
        {
            edges[other].currCost += l * PRECISION;
        }

        ans[i] = bestEdge;
    }
}

void solveIter()
{
    init();
    minSpanning();
    evalAns();
}

void solve()
{
    bestAns.clear();
    bestScore = -1;
    currT += 10;
    currAnsT += 10;

    for (int i = 0; i < m; ++i)
    {
        Edge& e = edges[i];
        e.distCost = distCost(ints[e.a], ints[e.b]);
        e.inters.clear();
    }

    for (int i = 0; i < m; ++i)
    {
        for (int j = i + 1; j < m; ++j)
        {
            if (edges[i].a == edges[j].a || edges[i].a == edges[j].b || edges[i].b == edges[j].a || edges[i].b == edges[j].b) continue;
            if (intersect(ints[edges[i].a], ints[edges[i].b], ints[edges[j].a], ints[edges[j].b]))
            {
                edges[i].inters.push_back(j);
                edges[j].inters.push_back(i);
            }
        }
    }

    for (int i = 0; i < 40; ++i)
    {
        solveIter();
    }

    ++currAnsT;
    ans = bestAns;
    localOpt();
    evalAns();
}

int main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    srand(time(NULL));
    srand(444);

    int t;
    cin >> t;
    for (int i = 0; i < t; ++i)
    {
        input();
        solve();
        cout << "case " << i + 1 << " Y\n";
        output();
    }

    return 0;
}
