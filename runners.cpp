#pragma GCC optimize ("Ofast")

#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <numeric>
#include <float.h>
#include <vector>
#include <set>
#include <algorithm>
#include <queue>
#include <tuple>
#include <iterator>
#include <assert.h>

#define FORCE_INLINE __attribute__((always_inline)) inline

constexpr bool DEBUG = true;
constexpr double TIME_MULT = 1;

std::ifstream inF("runners.in");
std::ofstream outF("runners.out");

std::mt19937 generator;
int randNum(int lb, int ub)
{
    std::uniform_int_distribution<int> distribution(lb, ub - 1);
    return distribution(generator);
}

std::chrono::high_resolution_clock::time_point startT, currT;

double timePassed()
{
    using namespace std::chrono;
    currT = high_resolution_clock::now();
    double time = duration_cast<duration<double>>(currT - startT).count();
    return time * TIME_MULT;
}

using Real = double;
const Real INF = 1e20;

const int MAX_N = 1e5;
const int MAX_K = 1e5;
const int MAX_XY = 1e9 - 1;

const int MAX_K_1 = 10;
const int MAX_K_2 = 1e2;
const int MAX_K_3 = 1e3;
const int MAX_K_4 = 1e4;
const int MAX_K_5 = 1e5;

struct Point
{
    int x, y;
};

int n, k;
Point comp[MAX_N];
Real speedUnsort[MAX_K];

const int MAX_DIGS = 6;

const Real EXP[MAX_DIGS + 1] = {
    1,
    1e-1,
    1e-2,
    1e-3,
    1e-4,
    1e-5,
    1e-6,
};

FORCE_INLINE
Real readSpeed()
{
    static std::string str;

    inF >> str;

    int size = str.size();

    if (size <= 2) return std::stoi(str);
    else if (str[1] != '.') return 10;

    int frac = std::stoi(str.substr(2));

    return str[0] - '0' + frac * EXP[size - 2];
}

void input()
{
    inF >> n >> k;

    for (int r = 0; r < k; ++r)
    {
        speedUnsort[r] = readSpeed();
    }

    for (int c = 0; c < n; ++c)
    {
        inF >> comp[c].x >> comp[c].y;
        --comp[c].x;
        --comp[c].y;
    }
}

Real speed[MAX_K];
int perm[MAX_K];
int invPerm[MAX_K];

bool cmpSpeed(int i, int j)
{
    return speedUnsort[i] < speedUnsort[j];
}

void sortSpeeds()
{
    std::iota(perm, perm + k, 0);
    std::sort(perm, perm + k, cmpSpeed);
    for (int r = 0; r < k; ++r)
    {
        speed[r] = speedUnsort[perm[r]];
        invPerm[perm[r]] = r;
    }
}

Real bestCost;
int bestStart[MAX_K];
int bestRunner[MAX_N];

void output()
{
    for (int r = 0; r < k; ++r)
    {
        int start = bestStart[invPerm[r]];
        outF << comp[start].x + 1 << " " << comp[start].y + 1 << "\n";
    }

    for (int c = 0; c < n; ++c)
    {
        outF << perm[bestRunner[c]] + 1 << "\n";
    }
}

Real gridSize;
int maxCellX;
int numCells;

FORCE_INLINE
int getCellX(int x)
{
    return x / gridSize;
}

const Real EPS = 1e-6;

void setGridSize(int runnersPerCell, bool useN)
{
    int numPoints = useN ? n : k;

    int freq = std::max((int) sqrt(numPoints / runnersPerCell), 2);

    gridSize = MAX_XY / (freq - EPS);
    maxCellX = getCellX(MAX_XY); 
    numCells = (maxCellX + 1) * (maxCellX + 1);
}

FORCE_INLINE
int getCell(int cellX, int cellY)
{
    return cellY + cellX * (maxCellX + 1);
}

FORCE_INLINE
int getClosestCellX(int x)
{
    int cellX = getCellX(x);
    if (cellX == 0) ++cellX;
    else if (cellX == maxCellX) --cellX;
    else if (x - cellX * gridSize < gridSize / 2) --cellX;
    else ++cellX;
    return cellX;
}

FORCE_INLINE
std::vector<int> getClosestCells(const Point& locP)
{
    int cellX = getCellX(locP.x);
    int cellY = getCellX(locP.y);
    int closestCellX = getClosestCellX(locP.x);
    int closestCellY = getClosestCellX(locP.y);
    return {
        getCell(cellX, cellY),
        getCell(closestCellX, cellY),
        getCell(cellX, closestCellY),
        getCell(closestCellX, closestCellY)
    };
}

Real getSpeed(int idx)
{
    return speed[idx];
}

const int NONE = -1;
const int UNKNOWN = -2;

FORCE_INLINE
bool isNone(int val)
{
    return val == NONE;
}

FORCE_INLINE
bool notNone(int val)
{
    return val != NONE;
}

struct Runner
{
    int idx;
    Real speed;
    int start;
    int loc;
    Point locP;
    int cell;
    int cellIdx;
    Real dist;
    Real totalDist;

    void reset(int newIdx)
    {
        idx = newIdx;
        speed = getSpeed(newIdx);
        start = NONE;
        loc = NONE;
        cell = NONE;
        cellIdx = NONE;
        dist = 0;
    }
};

struct Node
{
    Point locP;
    std::vector<int> closestCells;
    int cell;
    int next;
    int planRun;

    void reset(int idx)
    {
        locP = comp[idx];
        closestCells = getClosestCells(locP);
        cell = closestCells.front();
        next = NONE;
        planRun = NONE;
    }
};

const int MAX_CELLS = 2e4;

Runner runner[MAX_K];
Node node[MAX_N];
std::vector<int> runInCell[MAX_CELLS];

void updateSol()
{
    Real cost = 0;
    for (int r = 0; r < k; ++r)
    {
        cost += runner[r].totalDist * speed[r];
    }

    if (cost >= bestCost && bestCost > 0) return;

    bestCost = cost;

    std::cerr << " New best: " << bestCost << std::endl;

    for (int r = 0; r < k; ++r)
    {
        bestStart[r] = runner[r].start;
        bestRunner[bestStart[r]] = r;
    }

    for (int c = 0; c < n; ++c)
    {
        if (notNone(node[c].next)) bestRunner[node[c].next] = bestRunner[c];
    }
}

void resetCells()
{
    for (int cell = 0; cell < numCells; ++cell)
    {
        runInCell[cell].clear();
    }
}

void reset()
{
    resetCells();

    for (int c = 0; c < n; ++c)
    {
        node[c].reset(c);
    }

    for (int r = 0; r < k; ++r)
    {
        runner[r].reset(r);
    }
}

void initRun(Runner& curr);

void softReset()
{
    resetCells();

    for (int c = 0; c < n; ++c)
    {
        node[c].planRun = NONE;
    }

    for (int r = 0; r < k; ++r)
    {
        initRun(runner[r]);
    }
}

void finalize()
{
    for (int r = 0; r < k; ++r)
    {
        runner[r].totalDist = runner[r].dist;
    }
}

FORCE_INLINE
Real dist(const Point& a, const Point& b)
{
    Real dx = a.x - b.x;
    Real dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy);
}

FORCE_INLINE
Real dist(int loc1, int loc2)
{
    return notNone(loc1) && notNone(loc2) ? dist(node[loc1].locP, node[loc2].locP) : 0;
}

FORCE_INLINE
Real dist(const Runner& curr, int loc)
{
    return notNone(curr.loc) && notNone(loc) ? dist(curr.locP, node[loc].locP) : 0;
}

FORCE_INLINE
int getLocNext(int loc)
{
    return notNone(loc) ? node[loc].next : NONE;
}

FORCE_INLINE
int& getNext(Runner& curr)
{
    return notNone(curr.loc) ? node[curr.loc].next : curr.start;
}

FORCE_INLINE
void remFromCell(int cell, int cellIdx)
{
    if (isNone(cell) || isNone(cellIdx)) return;
    std::vector<int>& currIn = runInCell[cell];
    std::swap(currIn[cellIdx], currIn.back());
    runner[currIn[cellIdx]].cellIdx = cellIdx;
    currIn.pop_back();
}

FORCE_INLINE
int addToCell(int cell, int r)
{
    if (isNone(cell)) return NONE;
    std::vector<int>& currIn = runInCell[cell];
    currIn.push_back(r);
    return currIn.size() - 1;
}

FORCE_INLINE
void addNodeToCell(int c)
{
    runInCell[node[c].cell].push_back(c);
}

FORCE_INLINE
void updateRunCell(Runner& curr, int cell)
{
    remFromCell(curr.cell, curr.cellIdx);
    curr.cell = cell;
    curr.cellIdx = addToCell(curr.cell, curr.idx);
}

FORCE_INLINE
void updateRunCellLoc(Runner& curr, int loc)
{
    updateRunCell(curr, notNone(loc) ? node[loc].cell : NONE);
}

FORCE_INLINE
void initRun(Runner& curr)
{
    curr.dist = 0;
    curr.loc = NONE;
    curr.cell = NONE;
    curr.cellIdx = NONE;
    updateRunCellLoc(curr, curr.start);
    node[curr.start].planRun = curr.idx;
}

FORCE_INLINE
void pivotRun(Runner& curr, const Point& locP)
{
    curr.loc = UNKNOWN;
    curr.locP = locP;
    curr.cell = NONE;
    curr.cellIdx = NONE;
    int cell = getClosestCells(locP).front();
    updateRunCell(curr, cell);
}

FORCE_INLINE
void advanceRun(Runner& curr)
{
    int next = getNext(curr);
    if (isNone(next))
    {
        std::cerr << "ERROR: Advancing past end" << std::endl;
        return;
    }
    curr.dist += dist(curr, next);
    curr.loc = next;
    curr.locP = node[next].locP;
    updateRunCellLoc(curr, next);
}

FORCE_INLINE
void extendRun(Runner& curr, int loc)
{
    getNext(curr) = loc;
    advanceRun(curr);
}

FORCE_INLINE
void updatePlan(Runner& curr)
{
    int next = getNext(curr);
    if (notNone(next)) node[next].planRun = curr.idx;
}

FORCE_INLINE
void advancePlanRun(Runner& curr)
{
    advanceRun(curr);
    updatePlan(curr);
}

FORCE_INLINE
void swapRuns(Runner& left, Runner& right)
{
    int& leftNext = getNext(left);
    int& rightNext = getNext(right);

    Real leftNextD = dist(left, leftNext);
    Real rightNextD = dist(right, rightNext);

    Real leftCrossD = dist(left, rightNext);
    Real rightCrossD = dist(right, leftNext);

    Real leftRemD = left.totalDist - left.dist - leftNextD;
    Real rightRemD = right.totalDist - right.dist - rightNextD;

    left.totalDist = left.dist + leftCrossD + rightRemD;
    right.totalDist = right.dist + rightCrossD + leftRemD;

    std::swap(leftNext, rightNext);

    if (isNone(left.loc)) updateRunCellLoc(left, left.start);
    if (isNone(right.loc)) updateRunCellLoc(right, right.start);

    updatePlan(left);
    updatePlan(right);
}

#define FOR_CLOSEST(r, c) for (int cell : node[c].closestCells) for (int r : runInCell[cell])

int startPerm[MAX_N];

bool firstGen = true;

void genNN()
{
    std::iota(startPerm, startPerm + k, 0);
    if (!firstGen) std::shuffle(startPerm, startPerm + k, generator);
    firstGen = false;

    for (int r = 0; r < k; ++r)
    {
        extendRun(runner[r], startPerm[r]);
    }

    for (int c = k; c < n; ++c)
    {
        int minR = NONE;
        Real minT = INF;

        FOR_CLOSEST(r, c)
        {
            Real t = dist(runner[r], c) * runner[r].speed;
            if (t < minT)
            {
                minR = r;
                minT = t;
            }
        }

        if (isNone(minR))
        {
            std::cerr << "ERROR: No runners nearby" << std::endl;
            minR = 0;
        }

        extendRun(runner[minR], c);
    }
}

struct State
{
    Real totalTime;
    std::vector<int> pos;
    int lastRunner;
    int prevState;
};

namespace std
{
    FORCE_INLINE
    void swap(State& left, State& right)
    {
        std::swap(left.totalTime, right.totalTime);
        std::swap(left.pos, right.pos);
        std::swap(left.lastRunner, right.lastRunner);
        std::swap(left.prevState, right.prevState);
    }
}

FORCE_INLINE
bool operator<(State& left, State& right)
{
    return left.totalTime < right.totalTime;
}

const int MAX_NUM_STATES = 25;
const int MAX_STATES_K = 50;

int numStates;

State newStates[MAX_NUM_STATES * MAX_STATES_K];
std::vector<State> states[MAX_N];

int assignment[MAX_N];

void genMultState()
{
    for (int c = 0; c < n; ++c)
    {
        states[c].resize(numStates);
    }

    std::iota(startPerm, startPerm + k, 0);
    for (int s = 0; s < numStates; ++s)
    {
        auto& state = states[k - 1][s];

        state.totalTime = 0;
        state.pos.resize(k);
        std::copy(startPerm, startPerm + k, state.pos.begin());

        std::shuffle(startPerm, startPerm + k, generator);
    }

    for (int c = k; c < n; ++c)
    {
        int numNewStates = 0;
        auto& oldStates = states[c - 1];

        for (int r = 0; r < k; ++r)
        {
            for (int s = 0; s < numStates; ++s)
            {
                auto& state = oldStates[s];
                Real newTotalTime = state.totalTime + runner[r].speed * dist(state.pos[r], c);
                ++numNewStates;
                auto& newState = newStates[numNewStates - 1];
                newState.totalTime = newTotalTime;
                newState.pos.resize(k);
                std::copy(state.pos.begin(), state.pos.end(), newState.pos.begin());
                newState.pos[r] = c;
                newState.lastRunner = r;
                newState.prevState = s;
                int i = numNewStates - 1;
                while (i > 0 && newStates[i].totalTime < newStates[i - 1].totalTime)
                {
                    std::swap(newStates[i], newStates[i - 1]);
                    --i;
                }
            }
        }

        std::set<long long> usedHashes;

        int numNewFinalStates = 0;
        auto& newFinalStates = states[c];

        for (int s = 0; s < numNewStates && numNewFinalStates < numStates; ++s)
        {
            auto& newState = newStates[s];
            long long hash = 0;
            for (int r = 0; r < k; ++r)
            {
                hash = hash * n + newState.pos[r];
            }
            bool used = !usedHashes.insert(hash).second;
            if (used && numNewStates - s > numStates - numNewFinalStates) continue;
            std::swap(newFinalStates[numNewFinalStates++], newStates[s]);
        }
    }

    int s = 0;
    int c = n - 1;

    while (c >= k)
    {
        auto& state = states[c][s];
        assignment[c] = state.lastRunner;
        s = state.prevState;
        --c;
    }

    for (int r = 0; r < k; ++r)
    {
        extendRun(runner[r], states[c][s].pos[r]);
    }

    for (int c = k; c < n; ++c)
    {
        extendRun(runner[assignment[c]], c);
    }
}

struct MeanTracker
{
    long long x;
    long long y;
    long long cnt;

    void reset()
    {
        x = 0;
        y = 0;
        cnt = 0;
    }

    void add(const Point& locP)
    {
        x += locP.x;
        y += locP.y;
        cnt += 1;
    }

    Point get() const
    {
        if (cnt == 0) return node[randNum(0, n)].locP;
        int px = (x + cnt / 2) / cnt;
        int py = (y + cnt / 2) / cnt;
        return {px, py};
    }
};

MeanTracker meanTracker[MAX_K];

double kMeansTime;

void genKMeans()
{
    std::iota(startPerm, startPerm + k, 0);
    std::shuffle(startPerm, startPerm + k, generator);

    for (int r = 0; r < k; ++r)
    {
        meanTracker[r].reset();
        meanTracker[r].add(node[startPerm[r]].locP);
    }

    do
    {
        softReset();

        for (int r = 0; r < k; ++r)
        {
            pivotRun(runner[r], meanTracker[r].get());
            meanTracker[r].reset();
        }

        for (int c = 0; c < n; ++c)
        {
            int minR = NONE;
            Real minT = INF;

            FOR_CLOSEST(r, c)
            {
                Real t = dist(runner[r], c) * runner[r].speed;
                if (t < minT)
                {
                    minR = r;
                    minT = t;
                }
            }

            if (isNone(minR))
            {
                std::cerr << "ERROR: No runners nearby" << std::endl;
                minR = 0;
            }

            assignment[c] = minR;
            meanTracker[minR].add(node[c].locP);
        }
    }
    while (timePassed() < kMeansTime);

    for (int r = 0; r < k; ++r)
    {
        if (meanTracker[r].cnt > 0) continue;
        int c;
        do
        {
            c = randNum(0, n);
        }
        while (meanTracker[assignment[c]].cnt == 1);
        --meanTracker[assignment[c]].cnt;
        --meanTracker[r].cnt;
        assignment[c] = r;
    }

    softReset();

    for (int c = 0; c < n; ++c)
    {
        extendRun(runner[assignment[c]], c);
    }
}

const int MIN_MATCHING_K = 36000;

std::vector<std::tuple<Real, int, int>> edges;

void genEdges()
{
    for (int c = 0; c < n; ++c)
    {
        addNodeToCell(c);
    }

    for (int c = 0; c < n; ++c)
    {
        FOR_CLOSEST(c2, c)
        {
            if (c2 <= c) continue;
            Real d = dist(c, c2);
            edges.emplace_back(d, c, c2);
        }
    }

    std::sort(edges.begin(), edges.end());
}

bool hasPrev[MAX_N];

void genMatching()
{
    std::fill(hasPrev, hasPrev + n, false);

    int comps = n;

    for (auto& [d, from, to] : edges)
    {
        if (comps == k) break;

        if (hasPrev[to] || notNone(node[from].next)) continue;

        node[from].next = to;
        hasPrev[to]= true;
        --comps;
    }

    if (comps > k)
    {
        reset();
        genNN();
        return;
    }

    int nextRunner = 0;
    for (int c = 0; c < n; ++c)
    {
        if (!hasPrev[c]) assignment[c] = nextRunner++;
        if (notNone(node[c].next)) assignment[node[c].next] = assignment[c];
        extendRun(runner[assignment[c]], c);
    }
}

bool cmpDist(const Runner& left, const Runner& right)
{
    return left.dist > right.dist;
}

void optSortRunners()
{
    std::sort(runner, runner + k, cmpDist);
    for (int r = 0; r < k; ++r)
    {
        runner[r].idx = r;
        runner[r].speed = speed[r];
    }
}

bool opt2Opt()
{
    softReset();

    bool succ = false;

    for (int c = 0; c < n && timePassed() < 4.25; ++c)
    {
        assert(node[c].planRun >= 0 && node[c].planRun < k);

        int planIdx = node[c].planRun;
        Runner& plan = runner[planIdx];

        int minR = planIdx;
        bool minSwapSpeeds = false;
        bool minSwapAgain = false;
        Real minDt = 0;

        Real planNextD = dist(plan, c);
        Real planRemD = plan.totalDist - plan.dist - planNextD;

        Real planSkipD = dist(plan, getLocNext(c));
        Real planSkipRemD = planRemD - dist(c, getLocNext(c));

        FOR_CLOSEST(r, c)
        {
            if (r == planIdx) continue;

            assert(r >= 0 && r < k);

            Runner& curr = runner[r];

            int currNext = getNext(curr);

            Real currNextD = dist(curr, currNext);
            Real currRemD = curr.totalDist - curr.dist - currNextD;

            Real planCrossD = dist(plan, currNext);
            Real currCrossD = dist(curr, c);

            Real currRetD = dist(c, currNext);

            Real baseT = plan.totalDist * plan.speed + curr.totalDist * curr.speed;

            Real planD = plan.dist + planCrossD + currRemD;
            Real currD = curr.dist + currCrossD + planRemD;

            Real newT = planD * plan.speed + currD * curr.speed;

            if (newT - baseT < minDt)
            {
                minR = r;
                minSwapSpeeds = false;
                minSwapAgain = false;
                minDt = newT - baseT;
            }

            newT = planD * curr.speed + currD * plan.speed;

            if (newT - baseT < minDt)
            {
                minR = r;
                minSwapSpeeds = true;
                minSwapAgain = false;
                minDt = newT - baseT;
            }

            planD = plan.dist + planSkipD + planSkipRemD;
            currD = curr.dist + currCrossD + currRetD + currRemD;

            newT = planD * plan.speed + currD * curr.speed;

            if (newT - baseT < minDt)
            {
                minR = r;
                minSwapSpeeds = false;
                minSwapAgain = true;
                minDt = newT - baseT;
            }

            newT = planD * curr.speed + currD * plan.speed;

            if (newT - baseT < minDt)
            {
                minR = r;
                minSwapSpeeds = true;
                minSwapAgain = true;
                minDt = newT - baseT;
            }
        }

        Runner& newPlan = runner[minR];

        if (minR != planIdx)
        {
            succ = true;
            if (minSwapSpeeds) std::swap(plan.speed, newPlan.speed);
            swapRuns(plan, newPlan);
        }

        advancePlanRun(newPlan);

        if (minR != planIdx && minSwapAgain) swapRuns(plan, newPlan);
    }

    return succ;
}

template <class F>
void gen(F genFunc)
{
    reset();
    genFunc();
    finalize();
    optSortRunners();
    updateSol();
}

template <class F>
bool opt(F optFunc)
{
    if (!optFunc()) return false;
    optSortRunners();
    updateSol();
    return true;
}

void solve()
{
    bool doMatching = k >= MIN_MATCHING_K;

    if (doMatching)
    {
        setGridSize(5, true);
        reset();
        genEdges();
    }

    int runnersPerCell;
    if (k >= MIN_MATCHING_K) runnersPerCell = 10;
    else if (k >= 700) runnersPerCell = 15;
    else if (k >= 100) runnersPerCell = 10;
    else runnersPerCell = 5;

    setGridSize(runnersPerCell, false);

    numStates = std::min(100 / k, MAX_NUM_STATES);
    bool doMultState = k <= MAX_STATES_K;

    kMeansTime = 2.0;
    bool doKMeans = k > MAX_STATES_K && k <= 400;

    while (timePassed() < 3.75)
    {
        if (doMultState)
        {
            gen(genMultState); 
            doMultState = false;
        }
        else if (doKMeans)
        {
            gen(genKMeans);
            doKMeans = false;
        }
        else if (doMatching)
        {
            gen(genMatching);
            doMatching = false;
        }
        else gen(genNN);

        while (timePassed() < 4.25)
        {
            if (!opt(opt2Opt)) break;
        }
    }
}

int main()
{
    generator.seed(0);
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    startT = std::chrono::high_resolution_clock::now();

    input();
    sortSpeeds();
    solve();
    output();

    std::cerr << " Time passed: " << timePassed() << std::endl;

    return 0;
}
