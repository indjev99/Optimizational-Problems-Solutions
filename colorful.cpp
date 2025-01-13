#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <cmath>
#include <random>
#include <queue>
#include <unordered_set>
#include <map>
#include <chrono>
#include <cassert>

std::mt19937 generator;

int randInt(int lb, int rb)
{
    std::uniform_int_distribution<int> distribution(lb, rb - 1);
    return distribution(generator);
}

double randReal(double lb, double rb)
{
    std::uniform_real_distribution<double> distribution(lb, rb);
    return distribution(generator);
}

std::chrono::high_resolution_clock::time_point startT;

double timePassed()
{
    std::chrono::high_resolution_clock::time_point currT = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::duration<double>>(currT - startT).count();
}

struct Point
{
    int x;
    int y;
};

int pointDistSq(const Point& l, const Point& r)
{
    return
        (l.x - r.x) * (l.x - r.x) +
        (l.y - r.y) * (l.y - r.y);
}

struct Color
{
    int r;
    int g;
    int b;
};

int colDistSq(const Color& l, const Color& r)
{
    return
        (l.r - r.r) * (l.r - r.r) +
        (l.g - r.g) * (l.g - r.g) +
        (l.b - r.b) * (l.b - r.b);
}

struct QAns
{
    Point loc;
    int k;

    Color col;
    std::vector<Point> ps;
};

constexpr int MAX_COL = 255;
constexpr int MAX_K = 5000;
constexpr int MAX_Q = 5000;

struct PAns
{
    PAns() {}

    int k = 10 * MAX_K;
    Color col = {MAX_COL, MAX_COL, MAX_COL};

    Color extraCol;

    std::vector<int> qIs;
    std::vector<int> rejIs;
};

constexpr int MAX_NM = 512;

int n, m, maxK, maxQs;

PAns ans[MAX_NM][MAX_NM];

void input()
{
    std::cin >> n >> m >> maxK >> maxQs;
}

Color real[MAX_NM][MAX_NM];

void inputReal()
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            std::cin >> real[i][j].r >> real[i][j].g >> real[i][j].b;
        }
    }
}

void output()
{
    std::cout << "Ready" << "\n";

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            std::cout << ans[i][j].col.r << " " << ans[i][j].col.g << " " << ans[i][j].col.b << " ";
        }
        std::cout << "\n";
    }

    std::cout << std::flush;
}

void outputScore()
{
    long long totalScore = 0;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            totalScore += colDistSq(ans[i][j].col, real[i][j]);
        }
    }

    std::cout << std::fixed << std::setprecision(1) << (double) totalScore / n / m << std::endl;
    // std::cout << totalScore << std::endl;
}

int numQs = 0;

std::vector<QAns> qs;

void sendQOnline(QAns& q)
{
    std::cout << q.loc.x << " " << q.loc.y << " " << q.k << "\n" << std::flush;

    std::cin >> q.col.r >> q.col.g >> q.col.b;

    int s;
    std::cin >> s;

    q.ps.resize(s);

    for (Point& p : q.ps)
    {
        std::cin >> p.x >> p.y;
    }
}

bool locVis[MAX_NM][MAX_NM];

void sendQLocal(QAns& q)
{
    std::queue<Point> qu;

    q.col = real[q.loc.x][q.loc.y];

    locVis[q.loc.x][q.loc.y] = true;
    qu.push(q.loc);

    while (!qu.empty())
    {
        Point p = qu.front();
        qu.pop();

        q.ps.push_back(p);

        std::vector<Point> adj = {
            {p.x, p.y - 1},
            {p.x, p.y + 1},
            {p.x - 1, p.y},
            {p.x + 1, p.y},
        };

        for (Point next : adj)
        {
            if (next.x < 0 || next.x >= n || next.y < 0 || next.y >= m) continue;

            if (locVis[next.x][next.y]) continue;

            if (colDistSq(real[next.x][next.y], q.col) > q.k) continue;

            locVis[next.x][next.y] = true;
            qu.push(next);
        }
    }

    for (Point p : q.ps)
    {
        locVis[p.x][p.y] = false;
    }

    q.ps.shrink_to_fit();
}

void sendQ(QAns& q)
{

#ifndef LOCAL
    sendQOnline(q);
#endif

#ifdef LOCAL
    sendQLocal(q);
#endif

    numQs++;

    qs.push_back(q);
}

const int MAX_DIST = MAX_NM * 2 + 1;
const int MAX_ADJ = 4;

int dists[MAX_NM][MAX_NM];
int adjs[MAX_NM][MAX_NM];

std::vector<Point> distPq[MAX_DIST + 1];
std::vector<Point> adjsPq[MAX_ADJ + 1];

void init()
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            dists[i][j] = MAX_DIST;
            distPq[MAX_DIST].push_back({i, j});
            adjs[i][j] = 4 - (i == 0) - (i == n - 1) - (j == 0) - (j == m - 1);
        }
    }
}

void bfs(std::vector<Point>& origins)
{
    std::queue<Point> q;

    for (Point p : origins)
    {
        if (dists[p.x][p.y] == 0) continue;
        dists[p.x][p.y] = 0;
        q.push(p);
    }

    while (!q.empty())
    {
        Point p = q.front();
        q.pop();

        int nextD = dists[p.x][p.y] + 1;

        std::vector<Point> adj = {
            {p.x, p.y - 1},
            {p.x, p.y + 1},
            {p.x - 1, p.y},
            {p.x + 1, p.y},
        };

        for (Point next : adj)
        {
            if (next.x < 0 || next.x >= n || next.y < 0 || next.y >= m) continue;

            if (nextD == 1)
            {
                adjs[next.x][next.y]--;
                adjsPq[adjs[next.x][next.y]].push_back(next);
            }

            if (dists[next.x][next.y] <= nextD) continue;

            dists[next.x][next.y] = nextD;
            distPq[nextD].push_back(next);
            q.push(next);
        }
    }
}

const int NUM_Q_SMPL_CANDS = 7;
const double Q_SMPL_SZ_POW = 0.35;

double SCALE_K = 0.8;

const double SCALE_SZ_TARGET_MAX = 0.875;
const double SCALE_SZ_TARGET_MIN = 0.125;
const double SCALE_ADJ_DOWN = 1 / 1.02;
const double SCALE_ADJ_UP = 1.02;
const double SCALE_MAX = 0.8;
const double SCALE_MIN = 0.2;

const double COL_K_POW = 1.5;
const double COL_D_POW = 0.75;
const double COL_MC_SIM_POW = 1.5;
int FIX_MAX_D = 20;
const double FIX_D_POW = 4.0;

const double TARGET_COLOR_TIME_END = 8.0;

int MAX_VALIDS_ITERS = 750;

void fixOne(int x, int y, bool extra = false);
void colorOne(int x, int y);

void colorOrFixOne(int x, int y)
{
    if (ans[x][y].k > maxK) fixOne(x, y);
    else colorOne(x, y);
}

bool measureSurpr = false;

double makeQ(int x, int y)
{
    QAns q;

    q.loc.x = x;
    q.loc.y = y;

    PAns& pa = ans[q.loc.x][q.loc.y];

    if (pa.k == 0) return 0;

    q.k = pa.k <= maxK ? (int) std::round(pa.k * SCALE_K) : maxK;
    q.k = std::max(q.k, 1);

    if (measureSurpr)
    {
        colorOrFixOne(x, y);
    }

    sendQ(q);

    double surprize = colDistSq(pa.col, q.col);
    double totalSurprize = surprize;

    pa.k = 0;
    pa.col = q.col;

    int qNum = numQs - 1;

    if (measureSurpr)
    {
        std::shuffle(q.ps.begin(), q.ps.end(), generator);
    }

    for (auto p : q.ps)
    {
        PAns& pa = ans[p.x][p.y];

        if (measureSurpr && (p.x != x || p.y != y))
        {
            totalSurprize += surprize * std::min((double) pa.k / q.k, 1.0) / std::pow(pointDistSq(p, q.loc), COL_D_POW / 2);
        }

        pa.qIs.push_back(qNum);

        if (pa.k < q.k) continue;

        pa.k = q.k;
        pa.col = q.col;
    }

    for (auto p : q.ps)
    {
        std::vector<Point> adj = {
            {p.x, p.y - 1},
            {p.x, p.y + 1},
            {p.x - 1, p.y},
            {p.x + 1, p.y},
        };

        for (Point next : adj)
        {
            if (next.x < 0 || next.x >= n || next.y < 0 || next.y >= m) continue;

            PAns& pa = ans[next.x][next.y];

            if ((pa.qIs.empty() || pa.qIs.back() != qNum) &&
                (pa.rejIs.empty() || pa.rejIs.back() != qNum))
            {
                pa.rejIs.push_back(qNum);
            }
        }
    }

    bfs(q.ps);

    return totalSurprize;
}

int firstD1 = 0;

void makeQs()
{
    firstD1 = maxQs;

    int dist = MAX_DIST + 1;
    int numAdj = MAX_ADJ + 1;

    std::vector<Point> ps;

    const int NUM_RETRY_ITERS = 200;

    auto samplePoint = [&](){
        double bestScore = 0;
        double bestAvgSize = 0;
        Point bestPoint = {0, 0};

        for (int iter = 0, iter2 = 0; iter < NUM_Q_SMPL_CANDS && iter2 < NUM_RETRY_ITERS; iter++, iter2++)
        {
            int x = randInt(0, n);
            int y = randInt(0, m);

            PAns pans = ans[x][y];

            if (pans.k == 0 || pans.k > maxK)
            {
                --iter;
                continue;
            }

            double avgSize = 0;
            int avgSizeCnt = 0;

            for (int i : pans.qIs)
            {
                if (qs[i].k <= pans.k * 1.1)
                {
                    avgSize += qs[i].ps.size();
                    avgSizeCnt++;
                }
            }

            avgSize /= avgSizeCnt;

            double score = std::pow(avgSize, Q_SMPL_SZ_POW) * pans.k;

            if (score > bestScore)
            {
                bestScore = score;
                bestAvgSize = avgSize;
                bestPoint = {x, y};
            }
        }

        double surpr = makeQ(bestPoint.x, bestPoint.y);

        double szRatio = qs.back().ps.size() / bestAvgSize;

        // std::cerr << "  " << std::setprecision(2) << SCALE_K << " : " << std::setprecision(2) << szRatio << std::endl;

        if (szRatio > SCALE_SZ_TARGET_MAX) SCALE_K *= SCALE_ADJ_DOWN; 
        else if (szRatio < SCALE_SZ_TARGET_MIN) SCALE_K *= SCALE_ADJ_UP; 

        SCALE_K = std::min(std::max(SCALE_K, SCALE_MIN), SCALE_MAX);

        return surpr;
    };

    const double EXP_DEC = 0.925;
    const int LAST_USE_MAX = 20;

    double expDecSurpr[2] = {0, 0};
    double expDecW[2] = {0, 0};
    int lastUsed[2] = {LAST_USE_MAX, LAST_USE_MAX};

    while (numQs < maxQs)
    {
        while (numAdj >= 0 && ps.empty())
        {
            if (dist > 1)
            {
                dist--;
                if (dist == 1) firstD1 = numQs;
            }
            if (dist == 1 && numAdj >= 0) numAdj--;
            if (numAdj >= 0)
            {
                if (dist > 1) ps = distPq[dist];
                else ps = adjsPq[numAdj];
                std::shuffle(ps.begin(), ps.end(), generator);
            }
        }

        if (ps.empty()) break;

        while (!ps.empty())
        {
            Point p = ps.back();

            if (dists[p.x][p.y] == dist && (dist != 1 || adjs[p.x][p.y] == numAdj)) break;

            ps.pop_back();
        }

        if (ps.empty()) continue;

        measureSurpr = dist <= 2;

        if (!measureSurpr)
        {
            Point p = ps.back();
            ps.pop_back();

            makeQ(p.x, p.y);

            continue;
        }

        int choice;

        if (lastUsed[0] == LAST_USE_MAX) choice = 0;
        else if (lastUsed[1] == LAST_USE_MAX) choice = 1;
        else if (expDecSurpr[0] / expDecW[0] >= expDecSurpr[1] / expDecW[1]) choice = 0;
        else choice = 1;

        lastUsed[choice] = 0;

        lastUsed[0]++;
        lastUsed[1]++;

        double surpr;

        if (choice == 0)
        {
            Point p = ps.back();
            ps.pop_back();
            surpr = makeQ(p.x, p.y);
        }
        else
        {
            surpr = samplePoint();
        }

        expDecSurpr[choice] = expDecSurpr[choice] * EXP_DEC + surpr;
        expDecW[choice] = expDecW[choice] * EXP_DEC + 1;

        // std::cerr << "qs: " << numQs << " choice: " << choice << " k: " << qs.back().k << " hits: " << qs.back().ps.size() << " sprz: " << surpr << std::endl; 
    }

    while (numQs < maxQs)
    {
        samplePoint();
    }
}

void fixOne(int x, int y, bool extra)
{
    double totalR = 0;
    double totalG = 0;
    double totalB = 0;
    double totalW = 0;

    for (int i = std::max(0, x - FIX_MAX_D); i < std::min(n, x + FIX_MAX_D + 1); i++)
    {
        for (int j = std::max(0, y - FIX_MAX_D); j < std::min(m, y + FIX_MAX_D + 1); j++)
        {
            if (ans[i][j].k > maxK) continue;
            if (i == x && j == y) continue;

            int distSq = pointDistSq({x, y}, {i, j});

            if (distSq > FIX_MAX_D * FIX_MAX_D) continue;

            double w = 1.0 / std::pow(distSq, FIX_D_POW / 2);

            totalW += w;
            totalR += w * ans[i][j].col.r;
            totalG += w * ans[i][j].col.g;
            totalB += w * ans[i][j].col.b;
        }
    }

    if (!extra && numQs == maxQs)
    {
        double totalW2 = 0;
        double totalR2 = 0;
        double totalG2 = 0;
        double totalB2 = 0;

        const int FIX_NB_MAX_D = 35;
        const double MIN_W2 = 2;
        const double BASE_MULT = 1;
        const double BIAS = 0.8;

        int start = firstD1;

        for (int i = start; i < numQs; i++)
        {
            if (qs[i].k != maxK) continue;

            int distSq = pointDistSq({x, y}, qs[i].loc);

            if (distSq > FIX_NB_MAX_D * FIX_NB_MAX_D) continue;

            double w = 1.0;

            totalW2 += w;
            totalR2 += w * qs[i].col.r;
            totalG2 += w * qs[i].col.g;
            totalB2 += w * qs[i].col.b;
        }

        if (totalW2 >= MIN_W2)
        {
            Color newCol = {
                (int) std::round(totalR2 / totalW2),
                (int) std::round(totalG2 / totalW2),
                (int) std::round(totalB2 / totalW2)
            };

            long long totalDSqOld = 0;
            long long totalDSqNew = 0;

            for (int i = start; i < numQs; i++)
            {
                if (qs[i].k != maxK) continue;

                int distSq = pointDistSq({x, y}, qs[i].loc);

                if (distSq > FIX_NB_MAX_D * FIX_NB_MAX_D) continue;

                totalDSqOld += colDistSq(qs[i].col, ans[qs[i].loc.x][qs[i].loc.y].extraCol);
                totalDSqNew += colDistSq(qs[i].col, newCol);
            }

            if (totalDSqNew * BIAS < totalDSqOld)
            {
                double mult = BASE_MULT * totalW / totalW2;

                totalW2 *= mult;
                totalR2 *= mult;
                totalG2 *= mult;
                totalB2 *= mult;

                totalW += totalW2;
                totalR += totalR2;
                totalG += totalG2;
                totalB += totalB2;
            }
        }
    }

    Color& col = !extra ? ans[x][y].col : ans[x][y].extraCol;

    if (totalW == 0)
    {
        col = {128, 128, 128};
        return;
    }

    totalR /= totalW;
    totalG /= totalW;
    totalB /= totalW;

    int r = std::max(std::min((int) std::round(totalR), MAX_COL), 0);
    int g = std::max(std::min((int) std::round(totalG), MAX_COL), 0);
    int b = std::max(std::min((int) std::round(totalB), MAX_COL), 0);

    col = {r, g, b};
}

int failCnt = 0;

void fixRIs(int x, int y)
{
    PAns& pans = ans[x][y];

    if (pans.rejIs.empty()) return;

    Color& col = pans.col;
    std::vector<int> rejIs = pans.rejIs;

    const int MAX_ITERS = 10;

    double minScore = 1e9;

    for (int iter = 0; iter < MAX_ITERS; iter++)
    {
        std::vector<Color> cands = {
            {col.r, col.g, col.b},
            {col.r - 1, col.g, col.b},
            {col.r + 1, col.g, col.b},
            {col.r, col.g - 1, col.b},
            {col.r, col.g + 1, col.b},
            {col.r, col.g, col.b - 1},
            {col.r, col.g, col.b + 1},
            {col.r - 1, col.g - 1, col.b},
            {col.r + 1, col.g - 1, col.b},
            {col.r, col.g - 1, col.b - 1},
            {col.r, col.g - 1, col.b + 1},
            {col.r - 1, col.g + 1, col.b},
            {col.r + 1, col.g + 1, col.b},
            {col.r, col.g + 1, col.b - 1},
            {col.r, col.g + 1, col.b + 1},
            {col.r - 1, col.g, col.b - 1},
            {col.r + 1, col.g, col.b - 1},
            {col.r - 1, col.g, col.b + 1},
            {col.r + 1, col.g, col.b + 1}
        };

        minScore = 1e9;
        Color bestCand = col;

        for (Color cand : cands)
        {
            if (cand.r < 0 || cand.r > MAX_COL || cand.g < 0 || cand.g > MAX_COL || cand.b < 0 || cand.b > MAX_COL)
            {
                continue;
            }

            double score = 0;

            for (int i : rejIs)
            {
                Color qCol = qs[i].col;

                int distSq = colDistSq(qCol, cand);

                if (distSq > qs[i].k) continue;

                score += 1 + std::sqrt(qs[i].k) - std::sqrt(distSq);
            }

            if (score < minScore)
            {
                minScore = score;
                bestCand = cand;
            }
        }

        if (bestCand.r == col.r && bestCand.g == col.g && bestCand.b == col.b) break;

        col = bestCand;
    }

    if (minScore > 0)
    {
        failCnt++;
    }
}

void fixAll()
{
    for (int i = 0; i < numQs; i++)
    {
        if (qs[i].k == maxK) fixOne(qs[i].loc.x, qs[i].loc.y, true);
    }

    int totalCnt = 0;
    for (int i = 0; i < n; i++)
    {
        if (timePassed() > 9.3) FIX_MAX_D = 8;
        for (int j = 0; j < m; ++j)
        {
            if (ans[i][j].k > maxK)
            {
                fixOne(i, j);
                fixRIs(i, j);
                totalCnt++;
            }
        }
    }

    // std::cerr << "fail rate: " << (double) 100 * failCnt / totalCnt << std::endl;
}

int currComp;
int srComp[MAX_NM][MAX_NM];
std::vector<std::unordered_set<int>> compRejs = {{}};

void spreadRejsBfs(int x, int y)
{
    currComp++;
    compRejs.push_back({});

    std::queue<Point> q;

    q.push({x, y});
    srComp[x][y] = currComp;

    while (!q.empty())
    {
        Point p = q.front();
        q.pop();

        for (int i : ans[p.x][p.y].rejIs)
        {
            compRejs[currComp].insert(i);
        }

        std::vector<Point> adj = {
            {p.x, p.y - 1},
            {p.x, p.y + 1},
            {p.x - 1, p.y},
            {p.x + 1, p.y},
        };

        for (Point next : adj)
        {
            if (next.x < 0 || next.x >= n || next.y < 0 || next.y >= m) continue;

            if (srComp[next.x][y]) continue;

            if (ans[next.x][next.y].qIs != ans[x][y].qIs) continue;

            srComp[next.x][next.y] = currComp;
            q.push(next);
        }
    }
}

void spreadRejs()
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            if (!srComp[i][j]) spreadRejsBfs(i, j);
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            ans[i][j].rejIs = std::vector<int>(compRejs[srComp[i][j]].begin(), compRejs[srComp[i][j]].end());
            std::sort(ans[i][j].rejIs.begin(), ans[i][j].rejIs.end());
        }
    }
}

void colorOne(int x, int y)
{
    PAns& pans = ans[x][y];

    if (pans.k == 0) return;
    if (pans.qIs.empty()) return;

    double totalR = 0;
    double totalG = 0;
    double totalB = 0;
    double totalW = 0;

    for (int i : pans.qIs)
    {
        int distSq = pointDistSq(qs[i].loc, {x, y});

        double w = 1.0 / std::pow(qs[i].k, COL_K_POW) / std::pow(distSq, COL_D_POW / 2.0);

        totalW += w;
        totalR += w * qs[i].col.r;
        totalG += w * qs[i].col.g;
        totalB += w * qs[i].col.b;
    }

    totalR /= totalW;
    totalG /= totalW;
    totalB /= totalW;

    int r = std::max(std::min((int) std::round(totalR), MAX_COL), 0);
    int g = std::max(std::min((int) std::round(totalG), MAX_COL), 0);
    int b = std::max(std::min((int) std::round(totalB), MAX_COL), 0);

    pans.col = {r, g, b};
}

bool isValid(const Color& col, const std::vector<int>& qIs, const std::vector<int>& rejIs)
{
    for (int i : qIs)
    {
        if (colDistSq(col, qs[i].col) > qs[i].k) return false;
    }

    for (int i : rejIs)
    {
        if (colDistSq(col, qs[i].col) <= qs[i].k) return false;
    }

    return true;
}

std::map<std::vector<int>, std::vector<Color>> qIsToValidsBase;

std::vector<Color> findValidsBase(const std::vector<int>& qIs)
{
    auto it = qIsToValidsBase.find(qIs);
    if (it != qIsToValidsBase.end()) return it->second;

    int minR = 0;
    int maxR = MAX_COL;
    int minG = 0;
    int maxG = MAX_COL;
    int minB = 0;
    int maxB = MAX_COL;

    for (int i : qIs)
    {
        int k = qs[i].k;
        Color qCol = qs[i].col;

        int dist = std::sqrt(k);

        minR = std::max(minR, qCol.r - dist);
        maxR = std::min(maxR, qCol.r + dist);
        minG = std::max(minG, qCol.g - dist);
        maxG = std::min(maxG, qCol.g + dist);
        minB = std::max(minB, qCol.b - dist);
        maxB = std::min(maxB, qCol.b + dist);
    }

    std::vector<Color> valids;

    for (int iter = 0; iter < MAX_VALIDS_ITERS; iter++)
    {
        Color col;

        col.r = randInt(minR, maxR + 1);
        col.g = randInt(minG, maxG + 1);
        col.b = randInt(minB, maxB + 1);

        if (!isValid(col, qIs, {})) continue;

        valids.push_back(col);
    }

    valids.shrink_to_fit();

    qIsToValidsBase.emplace(qIs, valids);

    return valids;
}

std::vector<Color> findValids(const std::vector<int>& qIs, const std::vector<int>& rejIs)
{
    std::vector<Color> validsBase = findValidsBase(qIs);
    std::vector<Color> valids;

    for (Color col : validsBase)
    {
        if (!isValid(col, {}, rejIs)) continue;

        valids.push_back(col);
    }

    return valids;
}

void cleanCache(int x0)
{
    std::map<std::vector<int>, std::vector<Color>> oldQIsToValidsBase;

    oldQIsToValidsBase = std::move(qIsToValidsBase);

    for (int x = x0; x < std::min(n, x0 + 2); x++)
    {
        for (int y = 0; y < m; y++)
        {
            if (ans[x][y].k == 0) continue;

            auto it = oldQIsToValidsBase.find(ans[x][y].qIs);
            if (it == oldQIsToValidsBase.end()) continue;
        
            qIsToValidsBase.emplace(it->first, std::move(it->second));

            oldQIsToValidsBase.erase(it);
        }
    }
}

void colorOneMC(int x, int y)
{
    if (y == 0 && x >= 0 && x % 50 == 0) cleanCache(x);

    PAns& pans = ans[x][y];

    if (pans.k == 0) return;
    if (pans.qIs.empty()) return;

    std::vector<Color> valids = findValids(pans.qIs, pans.rejIs);

    if (valids.empty())
    {
        return;
    }

    double totalR = 0;
    double totalG = 0;
    double totalB = 0;
    double totalW = 0;

    for (Color col : valids)
    {
        double w = 1.0 / (1.0 + std::pow(colDistSq(col, pans.col), COL_MC_SIM_POW));

        totalW += w;
        totalR += w * col.r;
        totalG += w * col.g;
        totalB += w * col.b;
    }

    totalR /= totalW;
    totalG /= totalW;
    totalB /= totalW;

    int r = std::max(std::min((int) std::round(totalR), MAX_COL), 0);
    int g = std::max(std::min((int) std::round(totalG), MAX_COL), 0);
    int b = std::max(std::min((int) std::round(totalB), MAX_COL), 0);

    pans.col = {r, g, b};
}

void colorAll()
{
    double colorTimePrev = timePassed();

    const int REC_FREQ = 20;

    for (int i = 0; i < n; i++)
    {
        if (i > 0 && i % REC_FREQ == 0)
        {
            double colorTimeNow = timePassed();

            double delta = colorTimeNow - colorTimePrev;

            double remTimeProj = delta * (n - i) / REC_FREQ;

            double targetRemTime = std::max(TARGET_COLOR_TIME_END - colorTimeNow, 0.0);

            if (i > REC_FREQ && std::abs(remTimeProj - targetRemTime) > 0.2)
            {
                MAX_VALIDS_ITERS = std::round(MAX_VALIDS_ITERS * targetRemTime / remTimeProj);
            }

            colorTimePrev = colorTimeNow;
        }

        for (int j = 0; j < m; ++j)
        {
            colorOne(i, j);
            colorOneMC(i, j);
        }
    }

    // std::cerr << "fail rate 2: " << (double) 100 * failCnt2 / n / m << std::endl;
}

void solve()
{
    init();

    makeQs();

    spreadRejs();

    colorAll();

    fixAll();


    // int maxDist = -1;

    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < m; ++j)
    //     {
    //         maxDist = std::max(dists[i][j], maxDist);
    //     }
    // }

    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < m; ++j)
    //     {
    //         int comp = MAX_COL - MAX_COL * dists[i][j] / maxDist;
    //         ans[i][j].col = {comp, comp, comp};
    //     }
    // }
}

int main()
{
    startT = std::chrono::high_resolution_clock::now();

    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    generator.seed(0);

    input();

#ifdef LOCAL
    inputReal();
#endif

    solve();

#ifndef LOCAL
    output();
#endif

#ifdef LOCAL
    outputScore();
#endif

    // std::cerr << "Total time: " << timePassed() << std::endl;

    return 0;
}
