#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <numeric>
#include <limits>
#include <vector>
#include <set>
#include <algorithm>
#include <queue>
#include <iterator>
#include <assert.h>

std::ifstream inF("cheating.in");
std::ofstream outF("cheating.out");

std::mt19937 generator;
int randNum(int lb, int ub)
{
    std::uniform_int_distribution<int> distribution(lb, ub);
    return distribution(generator);
}

typedef unsigned long long ull;

const int SMALL_N = 200;
const int MID_N = 2000;
const int MAX_N = 100000;
const ull INF = ULLONG_MAX;

int n, m;
ull pCost, qCost;
int page[MAX_N + 1];

struct Solution
{
    ull cost = INF;
    int path[MAX_N + 1];

    Solution& operator=(const Solution& other)
    {
        cost = other.cost;
        std::copy(other.path + 1, other.path + m + 1, path + 1);
        return *this;
    }

    bool operator<(const Solution& other) const
    {
        return cost < other.cost;
    }
};

Solution bestSol;

void input()
{
    inF >> n >> m >> pCost >> qCost;
    page[0] = 0;
    for (int i = 1; i <= n; ++i)
    {
        inF >> page[i];
    }
}

void output()
{
    std::cerr << "Best: " << bestSol.cost << std::endl;

    for (int i = 1; i <= m; ++i)
    {
        if (i > 1) outF << " ";
        outF << bestSol.path[i];
    }
    outF << std::endl;
}

double timeLimit = 4.8;
double maxTime;
std::chrono::high_resolution_clock::time_point startT, currT;
bool haveTime(bool urgent=false)
{
    using namespace std::chrono;
    currT = high_resolution_clock::now();
    double time = duration_cast<duration<double>>(currT - startT).count();
    if (urgent) return time < timeLimit;
    else return time < maxTime;
}

inline ull dist(int a, int b)
{
    return qCost * std::abs(a - b) + pCost * std::abs(page[a] - page[b]);
}

inline ull dist(const int qp1[2], const int qp2[2])
{
    return qCost * std::abs(qp1[0] - qp2[0]) + pCost * std::abs(qp1[1] - qp2[1]);
}

struct Node
{
    int qp[2];
    Node* left;
    Node* right;

    Node(const int qp[2]):
        qp{qp[0], qp[1]},
        left(nullptr),
        right(nullptr) {}
};

Node* newNode(const int qp[2])
{
    return new Node(qp);
}

inline bool cmpQP(const int qp1[2], const int qp2[2], bool dir)
{
    return qp1[dir] < qp2[dir] || (dir && qp1[dir] == qp2[dir] && qp1[0] < qp2[0]);
}

void insertNode(Node*& curr, const int qp[2], bool dir=false)
{
    if (curr == nullptr)
    {
        curr = newNode(qp);
        return;
    }

    if (cmpQP(qp, curr->qp, dir)) insertNode(curr->left, qp, !dir);
    else insertNode(curr->right, qp, !dir);
}

Node* findMinNode(Node* curr, bool targetDir, bool dir=false)
{
    if (curr == nullptr) return nullptr;

    if (dir == targetDir)
    {
        if (curr->left == nullptr) return curr;
        return findMinNode(curr->left, targetDir, !dir);
    }

    Node* min = curr;
    Node* minL = findMinNode(curr->left, targetDir, dir);
    Node* minR = findMinNode(curr->right, targetDir, dir);
    if (minL != nullptr && cmpQP(minL->qp, min->qp, targetDir)) min = minL;
    if (minR != nullptr && cmpQP(minR->qp, min->qp, targetDir)) min = minR;
    return min;
}

Node* findMaxNode(Node* curr, bool targetDir, bool dir=false)
{
    if (curr == nullptr) return nullptr;

    if (dir == targetDir)
    {
        if (curr->right == nullptr) return curr;
        return findMaxNode(curr->right, targetDir, !dir);
    }

    Node* max = curr;
    Node* maxL = findMaxNode(curr->left, targetDir, dir);
    Node* maxR = findMaxNode(curr->right, targetDir, dir);
    if (maxL != nullptr && cmpQP(max->qp, maxL->qp, targetDir)) max = maxL;
    if (maxR != nullptr && cmpQP(max->qp, maxR->qp, targetDir)) max = maxR;
    return max;
}

void removeNode(Node*& curr, const int qp[2], bool dir=false)
{
    if (curr == nullptr)
    {
        assert(false);
        return;
    }

    if (qp[0] == curr->qp[0])
    {
        if (curr->right != nullptr)
        {
            Node* min = findMinNode(curr->right, dir, !dir);
            curr->qp[0] = min->qp[0]; 
            curr->qp[1] = min->qp[1];
            removeNode(curr->right, min->qp, !dir);
        }
        else if (curr->left != nullptr)
        {
            Node* max = findMaxNode(curr->left, dir, !dir);
            curr->qp[0] = max->qp[0]; 
            curr->qp[1] = max->qp[1];
            removeNode(curr->left, max->qp, !dir);
        }
        else curr = nullptr;
        return;
    }

    if (cmpQP(qp, curr->qp, dir)) removeNode(curr->left, qp, !dir);
    else removeNode(curr->right, qp, !dir);
}

std::pair<ull, int> closestNode(Node*& curr, const int qp[2], ull minDist=INF, bool dir=false)
{
    if (curr == nullptr) return {0, 0};

    Node* close;
    Node* far;
    if (qp[dir] < curr->qp[dir])
    {
        close = curr->left;
        far = curr->right;
    }
    else
    {
        close = curr->right;
        far = curr->left;
    }

    std::pair<ull, int> min = {dist(qp, curr->qp), curr->qp[0]};
    minDist = std::min(minDist, min.first);
    std::pair<ull, int> minC = closestNode(close, qp, minDist, !dir);
    if (minC.second != 0 && minC < min) min = minC;
    minDist = std::min(minDist, min.first);
    if ((dir ? pCost : qCost) * std::abs(qp[dir] - curr->qp[dir]) <= minDist)
    {
        std::pair<ull, int> minF = closestNode(far, qp, minDist, !dir);
        if (minF.second != 0 && minF < min) min = minF;
    }
    return min;
}

int countNodes(Node*& curr)
{
    if (curr == nullptr) return 0;
    return 1 + countNodes(curr->left) + countNodes(curr->right);
}

void printNodes(Node*& curr, int depth=0)
{
    for (int i = 0; i < depth; ++i)
    {
        std::cerr << " ";
    }
    if (curr == nullptr)
    {
        std::cerr << ":" << std::endl;
        return;
    }
    std::cerr << ": " << curr->qp[0] << " " << curr->qp[1] << std::endl;
    if (curr->left != nullptr || curr->right != nullptr)
    {
        printNodes(curr->left, depth + 1);
        printNodes(curr->right, depth + 1);
    }
    if (depth == 0) std::cerr << std::endl;
}

Node* root = nullptr;

void insertQ(int q)
{
    int qp[2] = {q, page[q]};
    insertNode(root, qp);
}

void removeQ(int q)
{
    int qp[2] = {q, page[q]};
    removeNode(root, qp);
}

int closestQ(int q)
{
    int qp[2] = {q, page[q]};
    std::pair<ull, int> res = closestNode(root, qp);
    assert(res.second != 0);
    assert(res.first == dist(q, res.second));
    return res.second;
}

int qs[MAX_N + 1];

void init()
{
    std::iota(qs, qs + n + 1, 0);
    std::shuffle(qs + 1, qs + n + 1, generator);
    for (int i = 1; i <= n; ++i)
    {
        insertQ(qs[i]);
    }
    std::iota(qs, qs + n + 1, 0);
}

int genQGreedy(int pos)
{
    int prevQ = qs[0];
    qs[0] = closestQ(prevQ);
    removeQ(qs[0]);
    return qs[0];
}

const double RES = 1e8;

int optIters;

double explorRate;

const double alpha = 1;
const double beta = 1;

double baseCost;

const double minPher = 1e-6;

double initPherMult;
double evapRate;

double regAddRate;
double eliteAddRate;

double pher[MID_N + 1][MID_N + 1];

void initPher(double initPher)
{
    for (int i = 0; i <= n; ++i)
    {
        for (int j = 0; j <= n; ++j)
        {
            pher[i][j] = initPher;
        }
    }
}

void addPher(const Solution& sol, double amnt)
{
    for (int i = 0; i < m; ++i)
    {
        pher[sol.path[i]][sol.path[i + 1]] += amnt;
        pher[sol.path[i + 1]][sol.path[i]] += amnt;
    }
}

void updatePher(const Solution& sol)
{
    for (int i = 0; i <= n; ++i)
    {
        for (int j = 0; j <= n; ++j)
        {
            pher[i][j] = std::min(minPher, (1 - evapRate));
        }
    }

    addPher(bestSol, eliteAddRate * baseCost / bestSol.cost);
    addPher(sol, regAddRate * baseCost / sol.cost);
}

double weight[MID_N + 1];

int genQ(int pos)
{
    int prevQ = qs[pos - 1];
    int next = 0;
    if (randNum(0, RES - 1) < explorRate * RES)
    {
        double sum = 0;
        for (int i = pos; i <= n; ++i)
        {
            ull currDist = dist(prevQ, qs[i]);
            weight[i] = sum + pher[prevQ][qs[i]] / currDist;
            sum = weight[i];
        }
        assert(sum > 0);
        double mult = RES / sum;
        int roll = randNum(0, RES - 1);
        for (int i = pos; i <= n; ++i)
        {
            if (roll < weight[i] * mult)
            {
                next = i;
                break;
            }
        }
    }
    else
    {
        ull bestDist = INF;
        for (int i = pos; i <= n; ++i)
        {
            ull currDist = dist(prevQ, qs[i]);
            if (currDist < bestDist)
            {
                bestDist = currDist;
                next = i;
            }
        }
    }
    assert(next != 0);
    std::swap(qs[pos], qs[next]);
    return qs[pos];
}

void genSol(Solution& sol, bool greedy)
{
    qs[0] = 0;
    sol.cost = 0;
    for (int i = 1; i <= m; ++i)
    {
        sol.path[i] = greedy ? genQGreedy(i) : genQ(i);
        sol.cost += dist(sol.path[i - 1], sol.path[i]);
    }
    sol.cost -= qCost;
}

inline ull swapChange(const int path[], int i, int j)
{
    ull curr = dist(path[i], path[i + 1]) + dist(path[j], path[j + 1]);
    ull alt = dist(path[i], path[j]) + dist(path[i + 1], path[j + 1]);
    if (curr <= alt) return 0;
    return curr - alt;
}

void swapEdges(int path[], int i, int j)
{
    int k = (j - i) / 2;
    for (int t = 0; t < k; ++t)
    {
        std::swap(path[i + t + 1], path[j - t]);
    }
}

#define SLOW 0
#define NORMAL 1
#define FAST 2

int optMode;

int bestI, bestJ;
ull bestChange;

void one2OptSlow(const Solution& sol)
{
    for (int i = 0; i < m; ++i)
    {
        for (int j = i + 2; j < m; ++j)
        {
            ull change = swapChange(sol.path, i, j);
            if (change > bestChange)
            {
                bestChange = change;
                bestI = i;
                bestJ = j;
            }
        }
    }
}

#define END 0
#define START 1

struct Event
{
    int type;
    int line;
    int line2;

    Event(int type, int line):
        type(type),
        line(line) {}
};

int pointIdxs[MAX_N + 1];
std::vector<Event> events[MAX_N + 1];
std::vector<int> currLines;

void addLineNormal(const Solution& sol, int line)
{
    for (int i = 0; i < (int) currLines.size(); ++i)
    {
        int cl = currLines[i];
        assert(cl != line);
        if (cl == line - 1 || cl == line + 1) continue;

        int currI = std::min(line, cl);
        int currJ = std::max(line, cl);
        ull change = swapChange(sol.path, currI, currJ);
        if (change > bestChange)
        {
            bestChange = change;
            bestI = currI;
            bestJ = currJ;
        }
    }
    currLines.push_back(line);
}

void remLineNormal(int line)
{
    for (int i = 0; i < (int) currLines.size(); ++i)
    {
        if (currLines[i] == line)
        {
            std::swap(currLines[i], currLines.back());
            currLines.pop_back();
            return;
        }
    }
    assert(false);
}

bool careful;

void iterateEvents(const Solution& sol)
{
    for (int i = 0; i <= m; ++i)
    {
        for (Event ev : events[pointIdxs[i]])
        {
            if (ev.type == START) addLineNormal(sol, ev.line);
            else if (ev.type == END) remLineNormal(ev.line);
        }
    }
}

void carefullyIterateEvents(const Solution& sol)
{
    for (int i = 0; i <= m; ++i)
    {
        if (!haveTime(true)) break;
        for (Event ev : events[pointIdxs[i]])
        {
            if (ev.type == START) addLineNormal(sol, ev.line);
            else if (ev.type == END) remLineNormal(ev.line);
        }
    }
}

void one2OptNormal(const Solution& sol)
{
    events[sol.path[0]] = {};
    for (int i = 0; i < m; ++i)
    {
        if (sol.path[i] < sol.path[i + 1])
        {
            events[sol.path[i]].push_back(Event(START, i));
            events[sol.path[i + 1]] = {Event(END, i)};
        }
        else
        {
            events[sol.path[i]].push_back(Event(END, i));
            events[sol.path[i + 1]] = {Event(START, i)};
        }
    }

    assert(currLines.empty());
    if (careful) carefullyIterateEvents(sol);
    else iterateEvents(sol);
}

void full2Opt(Solution& sol, int iters=-1)
{
    if (optMode == NORMAL)
    {
        for (int i = 0; i <= m; ++i)
        {
            pointIdxs[i] = sol.path[i];
        }
        std::sort(pointIdxs, pointIdxs + m + 1);
    }
    while (iters-- && haveTime())
    {
        bestChange = 0;
        if (optMode == SLOW) one2OptSlow(sol);
        else if (optMode == NORMAL) one2OptNormal(sol);

        if (bestChange == 0) break;
        sol.cost -= bestChange;
        swapEdges(sol.path, bestI, bestJ);
    }
}

inline ull trapChange(const int path[], int i, int j, int k, int mode)
{
    ull curr = dist(path[i], path[i + 1]) + dist(path[j], path[j + 1]) + dist(path[k], path[k + 1]);
    ull alt = curr;
    if (k == m)
    {
        if (mode == 0)
        {
            alt = dist(path[i], path[j]) + dist(path[i + 1], path[j + 1]) + dist(path[k], path[k + 1]);
        }
    }
    else if (mode == 0)
    {
        alt = dist(path[i], path[j + 1]) + dist(path[k], path[i + 1]) + dist(path[j], path[k + 1]);
    }
    else if (mode == 1)
    {
        alt = dist(path[i], path[k]) + dist(path[j + 1], path[i + 1]) + dist(path[j], path[k + 1]);
    }
    else if (mode == 2)
    {
        alt = dist(path[i], path[j + 1]) + dist(path[k], path[j]) + dist(path[i + 1], path[k + 1]);
    }
    else if (mode == 3)
    {
        alt = dist(path[i], path[j]) + dist(path[i + 1], path[k]) + dist(path[j + 1], path[k + 1]);
    }
    if (curr <= alt) return 0;
    return curr - alt;
}

int temp[MAX_N + 1];

void trapEdges(int path[], int i, int j, int k, int mode)
{
    if (k == m)
    {
        assert(mode == 0);
        swapEdges(path, i, j);
        return;
    }

    std::copy(path, path + m + 1, temp);
    int curr = 0;
    for (int t = 0; t <= i; ++t)
    {
        path[curr++] = temp[t];
    }
    if (mode == 0)
    {
        for (int t = j + 1; t <= k; ++t)
        {
            path[curr++] = temp[t];
        }
        for (int t = i + 1; t <= j; ++t)
        {
            path[curr++] = temp[t];
        }
    }
    else if (mode == 1)
    {
        for (int t = k; t >= j + 1; --t)
        {
            path[curr++] = temp[t];
        }
        for (int t = i + 1; t <= j; ++t)
        {
            path[curr++] = temp[t];
        }
    }
    else if (mode == 2)
    {
        for (int t = j + 1; t <= k; ++t)
        {
            path[curr++] = temp[t];
        }
        for (int t = j; t >= i + 1; --t)
        {
            path[curr++] = temp[t];
        }
    }
    else if (mode == 3)
    {
        for (int t = j; t >= i + 1; --t)
        {
            path[curr++] = temp[t];
        }
        for (int t = k; t >= j + 1; --t)
        {
            path[curr++] = temp[t];
        }
    }
    for (int t = k + 1; t <= m; ++t)
    {
        path[curr++] = temp[t];
    }
}

int bestK, bestMode;

void one3Opt(const Solution& sol)
{
    for (int i = 0; i < m && haveTime(true); ++i)
    {
        for (int j = i + 2; j < m; ++j)
        {
            for (int k = j + 2; k <= m; ++k)
            {
                for (int mode = 0; mode <= 3; ++mode)
                {
                    ull change = trapChange(sol.path, i, j, k, mode);
                    if (change > bestChange)
                    {
                        bestChange = change;
                        bestI = i;
                        bestJ = j;
                        bestK = k;
                        bestMode = mode;
                    }
                }
            }
        }
    }
}

void full3Opt(Solution& sol, int iters=-1)
{
    while (iters-- && haveTime())
    {
        bestChange = 0;
        one3Opt(sol);
        if (bestChange == 0) break;
        sol.cost -= bestChange;
        trapEdges(sol.path, bestI, bestJ, bestK, bestMode);
    }
}

void updateBestSol(const Solution& sol)
{
    if (sol < bestSol)
    {
        std::cerr << "New best: " << sol.cost << "\n";
        bestSol = sol;
    }
}

int numNotInSol;
int notInSol[MAX_N + 1];

inline ull repChange(const int path[], int i, int j)
{
    ull curr = dist(path[i - 1], path[i]) + (i < m ? dist(path[i], path[i + 1]) : 0);
    ull alt = dist(path[i - 1], notInSol[j]) + (i < m ? dist(notInSol[j], path[i + 1]) : 0);
    if (curr <= alt) return 0;
    return curr - alt;
}

void oneRepOpt(const Solution& sol)
{
    for (int i = 1; i <= m; ++i)
    {
        for (int j = 0; j < numNotInSol; ++j)
        {
            ull change = repChange(sol.path, i, j);
            if (change > bestChange)
            {
                bestChange = change;
                bestI = i;
                bestJ = j;
            }
        }
    }
}

void fullRepOpt(Solution& sol, int iters=-1)
{
    for (int i = 0; i <= m; ++i)
    {
        pointIdxs[i] = sol.path[i];
    }
    std::sort(pointIdxs, pointIdxs + m + 1);
    int j = 0;
    numNotInSol = 0;
    for (int i = 1; i <= n; ++i)
    {
        while (pointIdxs[j] < i) ++j;
        if (i == pointIdxs[j]) continue;
        notInSol[numNotInSol++] = i;
    }
    while (iters-- && haveTime())
    {
        bestChange = 0;
        oneRepOpt(sol);
        if (bestChange == 0) break;
        sol.cost -= bestChange;
        std::swap(sol.path[bestI], notInSol[bestJ]);
    }
}

Solution sol;

void solve()
{
    init();

    optMode = NORMAL;
    int itersDone = 0;

    if (n <= SMALL_N)
    {
        maxTime = 4.0;

        optIters = -1;
        explorRate = 0.35;

        initPherMult = 5;
        evapRate = 0.002;

        regAddRate = 0.2;
        eliteAddRate = 0.2;
    }
    else if (n <= MID_N)
    {
        maxTime = 3.5;

        optIters = 75;
        explorRate = 0.01;

        initPherMult = 5;
        evapRate = 0.0035;

        regAddRate = 0.35;
        eliteAddRate = 0.35;
    }
    else
    {
        maxTime = 3.8;
        timeLimit = 4.6;
        careful = true;

        optIters = -1;
    }

    genSol(sol, true);
    full2Opt(sol, optIters);
    updateBestSol(sol);

    if (n > MID_N) return;

    baseCost = sol.cost;
    initPher(sol.cost * initPherMult);

    do
    {
        genSol(sol, false);
        full2Opt(sol, optIters);
        updateBestSol(sol);
        updatePher(sol);

        ++itersDone;
    }
    while (haveTime());

    if (n <= SMALL_N)
    {
        maxTime = 4.3;
        fullRepOpt(bestSol);
        maxTime = 4.6;
        full3Opt(bestSol);
        optMode = SLOW;
    }
    else if (n <= MID_N)
    {
        maxTime = 4.0;
        fullRepOpt(bestSol);
        maxTime = 4.4;
        optMode = SLOW;
        full2Opt(bestSol);
        maxTime = 4.6;
        careful = true;
        optMode = NORMAL;
        full2Opt(bestSol);
    }

    std::cerr << "Iters: " << itersDone << std::endl;
}

int main()
{
    generator.seed(1);
    startT = std::chrono::high_resolution_clock::now();

    input();
    solve();
    output();

    return 0;
}
