#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <queue>
#include <chrono>

std::mt19937 rng;

int randInt(int lb, int ub)
{
    std::uniform_int_distribution<int> distr(lb, ub - 1);
    return distr(rng);
}

std::chrono::high_resolution_clock::time_point startT;

double timePassed()
{
    std::chrono::high_resolution_clock::time_point currT = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::duration<double>>(currT - startT).count();
}

const double TL = 7.5;

const int MAX_NM = 512;
const int MAX_SHOTS = 10000;
const int MAX_D = 100;
const int MAX_RGB = 255;

struct Col
{
    int r;
    int g;
    int b;
};

struct Point
{
    int x;
    int y;
};

struct Shot
{
    int x;
    int y;
    int d;
    int r;
    int g;
    int b;
};

int n;
int m;
Col img[MAX_NM][MAX_NM];

std::vector<Shot> shots;

void input()
{
    int seed;

    std::cin >> n;
    std::cin >> m;
    std::cin >> seed;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            std::cin >> img[i][j].r >> img[i][j].g >> img[i][j].b;
        }
    }
}

void output()
{
    std::cout << shots.size() << std::endl;

    for (auto& s : shots)
    {
        s.r = std::max(std::min(s.r, MAX_RGB), 0);
        s.g = std::max(std::min(s.g, MAX_RGB), 0);
        s.b = std::max(std::min(s.b, MAX_RGB), 0);

        std::cout << s.x << " " << s.y << " " << s.d << " " << s.r << " " << s.g << " " << s.b << "\n";
    }
}

double WS_EPS = 0;
double TARGET_RAT = 0.6;

int maxYOffs[MAX_D + 1][MAX_D];

void precompMaxYOffs()
{
    for (int d = 1; d <= MAX_D; d++)
    {
        for (int x = 0; x < d; x++)
        {
            for (int y = 0; y < d; y++)
            {
                if (x * x + y * y >= d * d) break;

                maxYOffs[d][x] = y;
            }
        }
    }
}

const int MAX_DIST_SEEN = MAX_NM * MAX_NM * 2;
double sqrts[MAX_DIST_SEEN + 1];

void precompSqrt()
{
    for (int i = 0; i <= MAX_DIST_SEEN; i++)
    {
        sqrts[i] = std::sqrt(i);
    }
}

bool useInitShots;
std::vector<Shot> initShots;

void generateInitShotsGrid()
{
    int side = 1;

    while (true)
    {
        side += 1;
        int maxOff = side / 2;
        int d = std::ceil(std::sqrt(2 * maxOff * maxOff) / 0.6);

        if (d > MAX_D)
        {
            side -= 1;
            break;
        }
    }

    int maxOff = side / 2;
    int d = std::ceil(std::sqrt(2 * maxOff * maxOff) / 0.6);

    for (int i = side / 2; i <= n + side / 2 - 1; i += side)
    {
        for (int j = side / 2; j <= m + side / 2 - 1; j += side)
        {
            int x = std::min(i, n - 1);
            int y = std::min(j, m - 1);
            initShots.push_back({x, y, d});
        }
    }

    std::cerr << side << " : " << d << " -> " << initShots.size() << std::endl;
}

std::vector<int> wastedShots;
double weight[MAX_NM][MAX_NM];

void colorShots()
{
    wastedShots.clear();
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            weight[i][j] = 1;
        }
    }

    for (int i = shots.size() - 1; i >= 0; i--)
    {
        Shot& s = shots[i];

        s.r = img[s.x][s.y].r;
        s.g = img[s.x][s.y].g;
        s.b = img[s.x][s.y].b;

        double tw = 0;
        double tr = 0;
        double tg = 0;
        double tb = 0;

        for (int x = std::max(s.x - s.d + 1, 0); x <= std::min(s.x + s.d - 1, n - 1); x++)
        {
            int myo = maxYOffs[s.d][std::abs(s.x - x)];
            for (int y = std::max(s.y - myo, 0); y <= std::min(s.y + myo, m - 1); y++)
            {
                int dx = x - s.x;
                int dy = y - s.y;
                int dist2 = dx * dx + dy * dy;

                if (weight[x][y] == 0) continue;

                double dist = sqrts[dist2];

                double cw = (dist / s.d - 0.6) * 2.5;
                cw = std::min(std::max(cw, 0.0), 1.0);
                cw = 1 - cw * cw;

                double w = weight[x][y] * cw;
                weight[x][y] *= 1 - cw;

                tw += w;
                tr += w * img[x][y].r;
                tg += w * img[x][y].g;
                tb += w * img[x][y].b;
            }
        }

        if (tw <= WS_EPS && (!useInitShots || i >= (int) initShots.size()))
        {
            wastedShots.push_back(i);
        }

        if (tw < 1e-6) continue;

        s.r = std::round(tr / tw);
        s.g = std::round(tg / tw);
        s.b = std::round(tb / tw);
    }

    std::reverse(wastedShots.begin(), wastedShots.end());
}

void removeWastedShots()
{
    std::vector<Shot> goodShots;
    int j = 0;
    for (int i : wastedShots)
    {
        while (j < i && j < (int) shots.size())
        {
            goodShots.push_back(shots[j]);
            j++;
        }
        j++;
    }
    while (j < (int) shots.size())
    {
        goodShots.push_back(shots[j]);
        j++;
    }
    std::swap(shots, goodShots);
}

double evalSol()
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            weight[i][j] = (MAX_RGB - img[i][j].r) * (MAX_RGB - img[i][j].r) + (MAX_RGB - img[i][j].g) * (MAX_RGB - img[i][j].g) + (MAX_RGB - img[i][j].b) * (MAX_RGB - img[i][j].b);
        }
    }

    for (const Shot& s : shots)
    {
        for (int x = std::max(s.x - s.d + 1, 0); x <= std::min(s.x + s.d - 1, n - 1); x++)
        {
            int myo = maxYOffs[s.d][std::abs(s.x - x)];
            for (int y = std::max(s.y - myo, 0); y <= std::min(s.y + myo, m - 1); y++)
            {
                int dx = x - s.x;
                int dy = y - s.y;
                int dist2 = dx * dx + dy * dy;

                double dist = sqrts[dist2];

                double cw = (dist / s.d - 0.6) * 2.5;
                cw = std::min(std::max(cw, 0.0), 1.0);
                cw = 1 - cw * cw;

                weight[x][y] *= 1 - cw;

                weight[x][y] += cw * (
                    (s.r - img[x][y].r) * (s.r - img[x][y].r) + (s.g - img[x][y].g) * (s.g - img[x][y].g) + (s.b - img[x][y].b) * (s.b - img[x][y].b)
                );
            }
        }
    }

    double totalErr = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            totalErr += weight[i][j];
        }
    }

    return totalErr;
}

struct Reg
{
    int x;
    int y;

    int sx;
    int sy;

    int d;

    double var;
};

bool operator<(const Reg& a, const Reg& b)
{
    if (a.d > 100 || b.d > 100) return a.d < b.d;
    return a.var < b.var;
}

int getD(int sx, int sy)
{
    int x = sx / 2;
    int y = sy / 2;
    int d = std::ceil(sqrts[x * x + y * y] / TARGET_RAT);
    d = std::max(d, 1);
    return d;
}

double getVar(const Reg& reg)
{
    int x0 = std::min(reg.x + reg.sx / 2, n - 1);
    int y0 = std::min(reg.y + reg.sy / 2, m - 1);
    int d = reg.d;

    double tw = 0;
    double tr = 0;
    double tg = 0;
    double tb = 0;

    for (int x = std::max(x0 - d + 1, 0); x <= std::min(x0 + d - 1, n - 1); x++)
    {
        int myo = maxYOffs[d][std::abs(x0 - x)];
        for (int y = std::max(y0 - myo, 0); y <= std::min(y0 + myo, m - 1); y++)
        {
            int dx = x - x0;
            int dy = y - y0;
            int dist2 = dx * dx + dy * dy;

            double dist = sqrts[dist2];

            double cw = (dist / d - 0.6) * 2.5;
            cw = std::min(std::max(cw, 0.0), 1.0);
            cw = 1 - cw * cw;

            tw += cw;
            tr += cw * img[x][y].r;
            tg += cw * img[x][y].g;
            tb += cw * img[x][y].b;
        }
    }

    int r = std::round(tr / tw);
    int g = std::round(tg / tw);
    int b = std::round(tb / tw);

    double tv = 0;

    for (int x = std::max(x0 - d + 1, 0); x <= std::min(x0 + d - 1, n - 1); x++)
    {
        int myo = maxYOffs[d][std::abs(x0 - x)];
        for (int y = std::max(y0 - myo, 0); y <= std::min(y0 + myo, m - 1); y++)
        {
            double mul = 1.0;

            if (x < reg.x) mul *= 0.5;
            if (x >= reg.x + reg.sx) mul *= 0.5;

            if (y < reg.y) mul *= 0.5;
            if (y >= reg.y + reg.sy) mul *= 0.5;

            int dx = x - x0;
            int dy = y - y0;
            int dist2 = dx * dx + dy * dy;

            double dist = sqrts[dist2];

            double cw = (dist / d - 0.6) * 2.5;
            cw = std::min(std::max(cw, 0.0), 1.0);
            cw = 1 - cw * cw;

            cw *= mul;

            tv += cw * ((img[x][y].r - r) * (img[x][y].r - r) + (img[x][y].g - g) * (img[x][y].g - g) + (img[x][y].b - b) * (img[x][y].b - b));
        }
    }

    return tv;
}

void finalizeReg(Reg& reg)
{
    reg.d = getD(reg.sx, reg.sy);
    reg.var = getVar(reg);
}

std::priority_queue<Reg> cachePq;

void generateShotsVar(int numShots = MAX_SHOTS, bool resetCache = false)
{
    shots.clear();

    std::priority_queue<Reg> pq;

    if (resetCache)
    {
        Reg init = {0, 0, n, m};
        finalizeReg(init);
        pq.push(init);
    }
    else
    {
        pq = cachePq;
    }

    if (useInitShots) numShots -= initShots.size();

    while ((int) pq.size() + 3 <= numShots)
    {
        Reg curr = pq.top();

        if (curr.sx == 1 && curr.sy == 1) break;

        pq.pop();

        if (curr.sx > curr.sy)
        {
            Reg u = {curr.x, curr.y, curr.sx / 2, curr.sy};
            Reg b = {curr.x + curr.sx / 2, curr.y, (curr.sx + 1) / 2, curr.sy};

            finalizeReg(u);
            finalizeReg(b);

            pq.push(u);
            pq.push(b);
        }
        else // if (curr.sy > curr.sx)
        {
            Reg l = {curr.x, curr.y, curr.sx, curr.sy / 2};
            Reg r = {curr.x, curr.y + curr.sy / 2, curr.sx, (curr.sy + 1) / 2};

            finalizeReg(l);
            finalizeReg(r);

            pq.push(l);
            pq.push(r);
        }
        // else
        // {
        //     Reg ul = {curr.x, curr.y, curr.sx / 2, curr.sy / 2};
        //     Reg ur = {curr.x, curr.y + curr.sy / 2, curr.sx / 2, (curr.sy + 1) / 2};
        //     Reg bl = {curr.x + curr.sx / 2, curr.y, (curr.sx + 1) / 2, curr.sy / 2};
        //     Reg br = {curr.x + curr.sx / 2, curr.y + curr.sy / 2, (curr.sx + 1) / 2, (curr.sy + 1) / 2};

        //     finalizeReg(ul);
        //     finalizeReg(ur);
        //     finalizeReg(bl);
        //     finalizeReg(br);

        //     pq.push(ul);
        //     pq.push(ur);
        //     pq.push(bl);
        //     pq.push(br);
        // }
    }

    if (resetCache)
    {
        cachePq = pq;
    }

    shots = initShots;

    while (!pq.empty())
    {
        Reg curr = pq.top();
        pq.pop();

        shots.push_back({curr.x + curr.sx / 2, curr.y + curr.sy / 2, curr.d});
    }
}

double bestErr = -1;
std::vector<Shot> bestShots;
double bestWsEps;
double bestTargetRat;

void solveOne()
{
    if (timePassed() > TL) return;

    useInitShots = WS_EPS > 0 || TARGET_RAT > 0.6;

    generateShotsVar(MAX_SHOTS, true);

    int lb = MAX_SHOTS;
    int rb = MAX_SHOTS * 1.2;

    while (rb - lb > 1)
    {
        int mid = (lb + rb) / 2;

        generateShotsVar(mid);
        colorShots();
        removeWastedShots();

        // std::cerr << mid << ": " << shots.size() << " " << wastedShots.size() << std::endl;

        if ((int) shots.size() <= MAX_SHOTS) lb = mid;
        else rb = mid;
    }

    generateShotsVar(lb);
    colorShots();
    removeWastedShots();

    // std::cerr << lb << ": " << shots.size() << " " << wastedShots.size() << std::endl;

    colorShots();

    double err = evalSol();

    // std::cerr << "  " << WS_EPS << " " << TARGET_RAT << ": " << err << std::endl;

    if (bestErr == -1 || err < bestErr)
    {
        bestErr = err;
        bestShots = shots;
        bestWsEps = WS_EPS;
        bestTargetRat = TARGET_RAT;
    }
}

void solve()
{
    precompMaxYOffs();
    precompSqrt();

    generateInitShotsGrid();

    WS_EPS = 0;
    TARGET_RAT = 0.6;
    solveOne();

    WS_EPS = 0.06;
    TARGET_RAT = 0.7;
    solveOne();

    WS_EPS = 0.2;
    TARGET_RAT = 0.7;
    solveOne();

    WS_EPS = 0.06;
    TARGET_RAT = 0.75;
    solveOne();

    WS_EPS = 0.2;
    TARGET_RAT = 0.75;
    solveOne();

    WS_EPS = 0.6;
    TARGET_RAT = 0.75;
    solveOne();

    WS_EPS = 0.2;
    TARGET_RAT = 0.8;
    solveOne();

    WS_EPS = 0.6;
    TARGET_RAT = 0.8;
    solveOne();

    WS_EPS = 2.0;
    TARGET_RAT = 0.8;
    solveOne();

    WS_EPS = 0.6;
    TARGET_RAT = 0.85;
    solveOne();

    WS_EPS = 2.0;
    TARGET_RAT = 0.85;
    solveOne();

    WS_EPS = 0.6;
    TARGET_RAT = 0.9;
    solveOne();

    WS_EPS = 2.0;
    TARGET_RAT = 0.9;
    solveOne();

    std::cerr << bestErr << std::endl;
    std::cerr << bestWsEps << std::endl;
    std::cerr << bestTargetRat << std::endl;

    shots = bestShots;
}

int main()
{
    startT = std::chrono::high_resolution_clock::now();

    std::ios::sync_with_stdio(false);

    input();
    solve();
    output();

    return 0;
}
