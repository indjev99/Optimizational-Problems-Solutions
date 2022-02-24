#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <numeric>
#include <limits>
#include <float.h>
#include <vector>
#include <set>
#include <algorithm>
#include <queue>
#include <iterator>
#include <assert.h>

const bool DEBUG = true;

std::ifstream inF("newyear.in");
std::ofstream outF("newyear.out");

std::mt19937 generator;
int randNum(int lb, int ub)
{
    std::uniform_int_distribution<int> distribution(lb, ub - 1);
    return distribution(generator);
}

double maxTimeLimit = 4.5;
std::chrono::high_resolution_clock::time_point startT, currT;

double timePassed()
{
    using namespace std::chrono;
    currT = high_resolution_clock::now();
    double time = duration_cast<duration<double>>(currT - startT).count();
    return time;
}

typedef unsigned long long ull;

const ull INF = ULLONG_MAX;

const int MAX_N = 5e3;
const int MAX_ML = 2e4;
const int MAX_K = 200;

struct GraphEdge
{
    int to;
    ull cost;
};

struct MarkovEdge
{
    int to;
    double prob;
};

struct Assignment
{
    int house = -1;
    double cost = INF;
};

int n, m, l, k;
std::vector<GraphEdge> graphEdges[MAX_N];
std::vector<MarkovEdge> markovEdges[MAX_N];

int houses[MAX_N];
Assignment assignments[MAX_N];

void input()
{
    inF >> n >> m >> l >> k;

    int from, to;
    ull cost;
    double prob;

    for (int i = 0; i < m; ++i)
    {
        inF >> from >> to >> cost;
        graphEdges[from - 1].push_back({to - 1, cost});
        graphEdges[to - 1].push_back({from - 1, cost});
    }

    for (int i = 0; i < l; ++i)
    {
        inF >> from >> to >> prob;
        markovEdges[from - 1].push_back({to - 1, prob});
    }
}


void output()
{
    for (int i = 0; i < k; ++i)
    {
        outF << houses[i] + 1 << " ";
    }
    outF << std::endl;

    for (int i = 0; i < n; ++i)
    {
        outF << assignments[i].house + 1 << " ";
    }
    outF << std::endl;
}

ull dists[MAX_N][MAX_N];
double costs[MAX_N][MAX_N];

typedef std::pair<ull, int> pqPair;

void findDists(int idx)
{
    std::priority_queue<pqPair, std::vector<pqPair>, std::greater<pqPair>> pq;
    std::fill(dists[idx], dists[idx] + n, INF);

    dists[idx][idx] = 0;
    pq.push({0, idx});

    while (!pq.empty())
    {
        pqPair curr = pq.top();
        pq.pop();

        if (curr.first > dists[idx][curr.second]) continue;

        for (const GraphEdge& edge : graphEdges[curr.second])
        {
            if (curr.first + edge.cost < dists[idx][edge.to])
            {
                dists[idx][edge.to] = curr.first + edge.cost;
                pq.push({curr.first + edge.cost, edge.to});
            }
        }
    }
}

int CONV_ITERS_1 = 50;
int CONV_ITERS_2 = 10;

double weights[MAX_N];
double currWeights[2][MAX_N];

void findWeights()
{
    bool curr = false;

    std::fill(currWeights[curr], currWeights[curr] + n, 0);
    currWeights[curr][0] = 1;

    if (n <= 500 && l <= 5e3)
    {
        CONV_ITERS_1 *= 20;
        CONV_ITERS_2 *= 20;
    }

    for (int i = 0; i < CONV_ITERS_1 + CONV_ITERS_2; ++i)
    {
        std::fill(currWeights[!curr], currWeights[!curr] + n, 0);

        for (int j = 0; j < n; ++j)
        {
            if (i >= CONV_ITERS_1) weights[j] += currWeights[curr][j];

            for (const MarkovEdge& edge : markovEdges[j])
            {
                currWeights[!curr][edge.to] += currWeights[curr][j] * edge.prob;
            }
        }

        curr = !curr;
    }
}

const double PROB_LAST = 0.001;

void findCosts(int idx)
{
    for (int i = 0; i < n; ++i)
    {
        costs[idx][i] = dists[idx][i];
        for (const MarkovEdge& edge : markovEdges[i])
        {
            costs[idx][i] += (1 - PROB_LAST) * edge.prob * dists[idx][edge.to];
        }
        costs[idx][i] *= weights[i];
    }
}

bool trySwap(int oldIdx, int newIdx)
{
    std::swap(houses[oldIdx], houses[newIdx]);

    int newHouse = houses[oldIdx];
    int oldHouse = houses[newIdx];

    std::vector<std::pair<int, int>> newAssignments;

    double delta = 0;
    for (int j = 0; j < n; ++j)
    {
        if (costs[newHouse][j] < assignments[j].cost)
        {
            delta += costs[newHouse][j] - assignments[j].cost;
            newAssignments.push_back({j, newHouse});
        }
        else if (assignments[j].house == oldHouse)
        {
            double best = -1;
            double bestCost = INF;
            for (int i = 0; i < k; ++i)
            {
                if (costs[houses[i]][j] < bestCost)
                {
                    best = houses[i];
                    bestCost = costs[houses[i]][j];
                }
            }
            delta += bestCost - assignments[j].cost;
            newAssignments.push_back({j, best});
        }
    }

    if (delta >= 0)
    {
        std::swap(houses[oldIdx], houses[newIdx]);
        return false;
    }

    std::cerr << "Improvement: " << delta / CONV_ITERS_2 << std::endl;

    for (const std::pair<int, int>& newAss : newAssignments)
    {
        assignments[newAss.first].house = newAss.second;
        assignments[newAss.first].cost = costs[newAss.second][newAss.first];
    }

    return true;
}

bool tryAdd(int newIdx)
{
    int newHouse = houses[newIdx];

    double baseDelta = 0;
    std::vector<double> deltas(n, 0);
    for (int j = 0; j < n; ++j)
    {
        if (costs[newHouse][j] < assignments[j].cost)
        {
            baseDelta += costs[newHouse][j] - assignments[j].cost;
        }
        else
        {
            double bestCost = costs[newHouse][j];
            for (int i = 0; i < k; ++i)
            {
                if (houses[i] == assignments[j].house) continue;
                if (costs[houses[i]][j] < bestCost)
                {
                    bestCost = costs[houses[i]][j];
                }
            }
            deltas[assignments[j].house] += bestCost - assignments[j].cost;
        }
    }

    int best = -1;
    double bestDelta = INF;
    for (int i = 0; i < k; ++i)
    {
        if (baseDelta + deltas[houses[i]] < bestDelta)
        {
            best = i;
            bestDelta = baseDelta + deltas[houses[i]];
        }
    }

    if (bestDelta >= 0)
    {
        return false;
    }

    return trySwap(best, newIdx);
}

double bestTotalCost = INF;
int bestHouses[MAX_N];

void generateSolution(int numCands, int t)
{
    for (int i = 0; i < n; ++i)
    {
        assignments[i].house = -1;
        assignments[i].cost = INF;
    }

    int s = 0;

    while (timePassed() < 4.6 && s < std::min(t, k))
    {
        if (numCands > 5 && timePassed() > 4.4) numCands = 5;
        if (numCands > 1 && timePassed() > 4.5) numCands = 1;

        std::vector<int> cands;

        int best = -1;
        double bestTotalCost = INF;

        int p = s;

        if (numCands == 0) goto doAssignment;

        while (p - s < numCands && p < t)
        {
            int selected = randNum(p, t);
            std::swap(houses[p], houses[selected]);
            ++p;
        }

        if (numCands == 1) goto doAssignment;

        for (int i = s; i < p; ++i)
        {
            double totalCost = 0;
            int house = houses[i];

            for (int j = 0; j < n; ++j)
            {
                totalCost += std::min(costs[house][j], assignments[j].cost);
            }

            if (totalCost < bestTotalCost)
            {
                best = i;
                bestTotalCost = totalCost;
            }
        }

        std::swap(houses[s], houses[best]);

        doAssignment:

        int house = houses[s];

        for (int j = 0; j < n; ++j)
        {
            if (costs[house][j] < assignments[j].cost)
            {
                assignments[j].cost = costs[house][j];
                assignments[j].house = house;
            }
        }

        ++s;
    }

    double totalCost = 0;
    for (int i = 0; i < n; ++i)
    {
        totalCost += assignments[i].cost;
    }

    if (totalCost < bestTotalCost)
    {
        bestTotalCost = totalCost;
        std::copy(houses, houses + k, bestHouses);

        std::cerr << "New best edge cost: " << bestTotalCost / CONV_ITERS_2 << std::endl;
    }
}

int BASE_NUM_CANDS = 150;

double connectedness[MAX_N];

bool cmpHouses(int a, int b)
{
    return connectedness[a] > connectedness[b];
}

void solve()
{
    findWeights();

    std::iota(houses, houses + n, 0);

    double avgEdgeCost = 0.0;
    for (int i = 0; i < n; ++i)
    {
        for (auto edge : graphEdges[i])
        {
            avgEdgeCost += edge.cost;
        }
    }
    avgEdgeCost /= m * 2;

    for (int i = 0; i < n; ++i)
    {
        connectedness[i] = weights[i] * n / CONV_ITERS_2 * 3;
        for (auto edge : graphEdges[i])
        {
            connectedness[i] += 1.0 + avgEdgeCost / edge.cost / 3;
        }
    }

    std::sort(houses, houses + n, cmpHouses);

    int t = 0;

    while (timePassed() < 3.8 && t < n)
    {
        int house = houses[t];
        findDists(house);
        findCosts(house);
        ++t;
    }

    for (int i = t; i < k; ++i)
    {
        assignments[houses[i]].cost = 0;
        assignments[houses[i]].house = houses[i];
    } 

    int numCands;

    if (t <= k) numCands = 0;
    else numCands = std::min((ull) BASE_NUM_CANDS * MAX_N * MAX_K / n / k, (ull) t);

    generateSolution(numCands, t);

    if (n <= 20 && t == n)
    {
        while (timePassed() < 4.5)
        {
            generateSolution(1, t);
        }

        std::copy(bestHouses, bestHouses + k, houses);
        generateSolution(0, k);
    }

    while (timePassed() < 4.7 && t > k)
    {
        int newIdx = randNum(k, t);
        tryAdd(newIdx);
    }

    for (int i = 0; i < k; ++i)
    {
        assignments[houses[i]].house = houses[i];
    }

    double totalCost = 0;
    for (int i = 0; i < n; ++i)
    {
        totalCost += assignments[i].cost;
    }

    double edgeCost = totalCost / CONV_ITERS_2;

    std::cerr << "Final edge cost: " << edgeCost << std::endl;
}

int main()
{
    generator.seed(0);
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    startT = std::chrono::high_resolution_clock::now();

    input();
    solve();
    output();

    return 0;
}
