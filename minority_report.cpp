#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <algorithm>
#include <vector>
#include <numeric>
#include <math.h>
#include <fstream>
#include <chrono>
#include <random>
using namespace std;
using namespace chrono;

ifstream inF("minority_report.in");
ofstream outF("minority_report.out");

//#define ONLINE_JUDGE
#ifdef ONLINE_JUDGE
#define cout outF
#endif // ONLINE_JUDGE

#define cin inF

std::default_random_engine generator;
time_point<std::chrono::system_clock> startT, currT;

double timeElapsed()
{
    currT = high_resolution_clock::now() ;
    return (duration<double>(currT - startT)).count();
}

long long randNum(long long from, long long to)
{
    if (to < from) return from;
    std::uniform_int_distribution<long long> distribution(from, to);
    return distribution(generator);
}

struct EdgeTo
{
    int to, d;
};

struct Edge
{
    int from, to, d;
};

struct Crime
{
    int x, t, w;
};

const int MAX_C = 10000;
const int MAX_N = 1024;
const int INF = 1e8;

bool allOnes;

vector<Edge> allEs;
vector<EdgeTo> es[MAX_N];
vector<Crime> crs;

int n, e, p, c;

int dists[MAX_N][MAX_N];
int prevs[MAX_N][MAX_N];

struct Police
{
    int loc;
    int enterT;
    int leaveT;

    vector<int> route;
    vector<int> stays;

    void reset()
    {
        loc = -1;
        enterT = 0;
        route.clear();
        stays.clear();
    }

    bool canBeAtCrime(int crimeLoc, int crimeT) const
    {
        if (loc == -1) return true;
        return leaveT + dists[loc][crimeLoc] <= crimeT;
    }

    void beAtCrime(int crimeLoc, int crimeT)
    {
        if (loc != crimeLoc)
        {
            if (loc != -1)
            {
                stays.push_back(leaveT - enterT);
                enterT = leaveT + dists[loc][crimeLoc];
                addRoute(crimeLoc);
            }
            else
            {
                route.push_back(crimeLoc);
                loc = crimeLoc;
            }
        }
        leaveT = crimeT + 1;
    }

    void addRoute(int crimeLoc)
    {
        while (true)
        {
            loc = prevs[crimeLoc][loc];
            route.push_back(loc);
            if (loc != crimeLoc) stays.push_back(0);
            else break;
        }
    }
};

struct Solution
{
    int score;
    vector<Police> pols;
};

void swap(Solution& a, Solution& b)
{
    swap(a.score, b.score);
    swap(a.pols, b.pols);
}

Solution bestSol, sol;

void input()
{
    allOnes = true;
    int a, b, d, x, t, w;
    cin >> n >> e >> p >> c;
    for (int i = 0; i < e; ++i)
    {
        cin >> a >> b >> d;
        allEs.push_back({a, b, d});
        es[a].push_back({b, d});
        es[b].push_back({a, d});
    }
    for (int i = 0; i < c; ++i)
    {
        cin >> x >> t >> w;
        crs.push_back({x, t, w});
        if (w > 1) allOnes = false;
    }
}

void output()
{
    for (const Police& pol : bestSol.pols)
    {
        cout << pol.route.size() << "\n";
        for (int x : pol.route)
        {
            cout << x << " ";
        }
        cout << "\n";
        for (int t : pol.stays)
        {
            cout << t << " ";
        }
        cout << "\n";
    }
}

void findDists()
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            dists[i][j] = INF;
        }
    }

    for (const Edge e: allEs)
    {
        dists[e.from][e.to] = e.d;
        prevs[e.from][e.to] = e.from;
        dists[e.to][e.from] = e.d;
        prevs[e.to][e.from] = e.to;
    }

    for (int i = 0; i < n; ++i)
    {
        dists[i][i] = 0;
        prevs[i][i] = i;
    }

    for (int k = 0; k < n; ++k)
    {
        for (int i = 0; i < n; ++i)
        {
            if (dists[i][k] == INF) continue;
            for (int j = i + 1; j < n; ++j)
            {
                const int newDist = dists[i][k] + dists[k][j];
                if (dists[i][j] > newDist)
                {
                    dists[i][j] = newDist;
                    prevs[i][j] = prevs[k][j];
                    dists[j][i] = newDist;
                    prevs[j][i] = prevs[k][i];
                }
            }
        }
    }

    /*for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            cerr << i << " to " << j << " : " << dists[i][j] << " by " << prevs[i][j] << endl;
        }
    }*/
}

const int skipProb = 30;
const int shuffleProb = 100;
const int wProbMult = 4;
const int wProbUB = 14;
const int costDiv = 20;

void generateSol()
{
    sol.score = 0;
    vector<pair<int, int>> poss;
    for (Police& pol : sol.pols)
    {
        pol.reset();
    }
    bool skip = randNum(0, 1);
    for (const Crime& cr : crs)
    {
        int costType;
        if (randNum(0, shuffleProb - 1) == 0) costType = 0;
        else costType = randNum(1, 2);
        poss.clear();
        for (int i = 0; i < p; ++i)
        {
            const Police& pol = sol.pols[i];
            if (pol.canBeAtCrime(cr.x, cr.t))
            {
                int cost;
                if (costType == 0) cost = randNum(0, 10 * p);
                else if (costType == 1) cost = cr.t + 1 - pol.leaveT;
                else if (costType == 2)
                {
                    if (pol.loc == -1) cost = 0;
                    else cost = dists[pol.loc][cr.x];
                }
                poss.push_back({cost, i});
            }
        }

        if (poss.size() < cr.w) continue;

        int totalCost = 0;
        sort(poss.begin(), poss.end());
        for (int i = 0; i < cr.w; ++i)
        {
            totalCost += cr.t + 1 - sol.pols[poss[i].second].leaveT;
        }

        if (skip && allOnes && randNum(0, skipProb - 1) == 0) continue;
        //if (skip && !allOnes && randNum(0, cr.w * cr.w * wProbMult - 1) < wProbUB) continue;

        if (skip && !allOnes &&
            randNum(0, cr.w * cr.w * wProbMult - totalCost / costDiv - 1) < wProbUB) continue;

        for (int i = 0; i < cr.w; ++i)
        {
            sol.pols[poss[i].second].beAtCrime(cr.x, cr.t);
        }
        sol.score += cr.w * cr.w;

    }
    if (sol.score > bestSol.score)
    {
        swap(sol, bestSol);
    }
}

int dp[MAX_C];
int dpPrev[MAX_C];

void solveSingle()
{
    int best = -1;
    int bestLast = 0;
    for (int i = 0; i < c; ++i)
    {
        if (crs[i].w > 1) continue;

        dp[i] = 1;
        dpPrev[i] = -1;

        for (int j = 0; j < i; ++j)
        {
            if (crs[j].t + 1 + dists[crs[i].x][crs[j].x] <= crs[i].t && dp[j] + 1 > dp[i])
            {
                dp[i] = dp[j] + 1;
                dpPrev[i] = j;
            }
        }

        if (dp[i] > best)
        {
            best = dp[i];
            bestLast = i;
        }
    }

    std::vector<int> crimes;

    while (bestLast != -1)
    {
        crimes.push_back(bestLast);
        bestLast = dpPrev[bestLast];
    }

    reverse(crimes.begin(), crimes.end());

    sol.score = 0;
    Police& pol = sol.pols[0];
    pol.reset();
    for (int ci : crimes)
    {
        const Crime& cr = crs[ci];

        if (!pol.canBeAtCrime(cr.x, cr.t))
        {
            cerr << "warning can't be at crime " << ci << endl;
            continue;
        }
        pol.beAtCrime(cr.x, cr.t);
        ++sol.score;
    }

    swap(sol, bestSol);
}

void solve()
{
    findDists();

    bestSol.score = -1;
    bestSol.pols.resize(p);
    sol.pols.resize(p);

    if (p != 1)
    {
        while (timeElapsed() < 2.0)
        {
            generateSol();
        }
    }
    else
    {
        solveSingle();
    }
}

int main()
{
    startT = high_resolution_clock::now();
    generator.seed(0);

    input();
    solve();
    output();
    cerr << " Score: " << bestSol.score << endl;
    cerr << " Total time: " << timeElapsed() << endl;

    return 0;
}
