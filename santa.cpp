#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <map>
#include <stdlib.h>
#include <math.h>
#include <unordered_map>
using namespace std;

//#define DEBUG

int test;

const string tests[] = {
    "a",
    "b",
    "c",
    "d",
    "e",
    "f"
};

ifstream inFile;
ofstream outFile;

const int MAX_T = 1e4;
const int MAX_D = 1e2;
const int MAX_W = 10;
const int MAX_G = 1e4;

struct AccRange
{
    int maxW;
    int acc;
};

struct Gift
{
    int num;
    int score;
    int weight;
    int x, y;

    bool used;
};

int t, w, g;
long long d;
AccRange accRanges[MAX_W];
std::string giftNames[MAX_G];
Gift gifts[MAX_G];

#define AccUp 0
#define AccDown 1
#define AccLeft 2
#define AccRight 3
#define Float 4
#define LoadCarrots 5
#define LoadGift 6
#define DeliverGift 7

struct Action
{
    int code;
    int arg;

    Action() {};
    Action(int code, int arg): code(code), arg(arg) {};
};

struct Solution
{
    std::vector<Action> acts;
    int score;
};

Solution bestSol, sol;

void input()
{
    inFile >> t >> d >> w >> g;
    for (int i = 0; i < w; ++i)
    {
        inFile >> accRanges[i].maxW >> accRanges[i].acc;
    }
    for (int i = 0; i < g; ++i)
    {
        gifts[i].num = i;
        inFile >> giftNames[i] >> gifts[i].score >> gifts[i].weight >> gifts[i].x >> gifts[i].y;
    }
}

void output()
{
    outFile << bestSol.acts.size() << "\n";
    for (int i = 0; i < (int) bestSol.acts.size(); ++i)
    {
        switch (bestSol.acts[i].code)
        {
        case AccUp:
            outFile << "AccUp " << bestSol.acts[i].arg << "\n";
            break;
        case AccDown:
            outFile << "AccDown " << bestSol.acts[i].arg << "\n";
            break;
        case AccRight:
            outFile << "AccRight " << bestSol.acts[i].arg << "\n";
            break;
        case AccLeft:
            outFile << "AccLeft " << bestSol.acts[i].arg << "\n";
            break;
        case Float:
            outFile << "Float " << bestSol.acts[i].arg << "\n";
            break;
        case LoadCarrots:
            outFile << "LoadCarrots " << bestSol.acts[i].arg << "\n";
            break;
        case LoadGift:
            outFile << "LoadGift " << giftNames[bestSol.acts[i].arg] << "\n";
            break;
        case DeliverGift:
            outFile << "DeliverGift " << giftNames[bestSol.acts[i].arg] << "\n";
            break;
        }
    }
}

int getMaxAcc(int weight)
{
    for (int i = 0; i < w; ++i)
    {
        if (weight <= accRanges[i].maxW) return accRanges[i].acc;
    }
    return 0;
}

std::vector<int> getClose(int x, int y)
{
    std::vector<int> close;
    for (int i = 0; i < g; ++i)
    {
        if (gifts[i].used) continue;
        long long dx = x - gifts[i].x;
        long long dy = y - gifts[i].y;
        if (dx * dx + dy * dy <= d * d) close.push_back(i);
    }
    return close;
}

std::pair<int, int> pickAcc(int x, int y, int vx, int vy, int weight, int maxAcc)
{
    double best = 0;
    int bestGift = -1;

    if (maxAcc == 0) return {0, 0};

    for (int i = 0; i < g; ++i)
    {
        if (gifts[i].used) continue;
        if (gifts[i].weight > weight) continue;
        long long dx = x + vx - gifts[i].x;
        long long dy = y + vy - gifts[i].y;
        double dist = sqrt(dx * dx + dy * dy);
        double pot = gifts[i].score / dist;
        if (pot > best)
        {
            best = pot;
            bestGift = i;
        }
    }

    int tx = 0;
    int ty = 0;

    if (bestGift != -1)
    {
        tx = gifts[bestGift].x;
        ty = gifts[bestGift].y;
    }

    double tvx = (tx - x);
    double tvy = (ty - y);
    double dist = sqrt(tvx * tvx + tvy * tvy);

    // 3, 1.5 -- C
    // 5, 2.5 -- E, F

    if (dist / d > 12)
    {
        tvx *= (d * 3.5 + maxAcc * 8 + 1) / std::max<double>(dist, 1);
        tvy *= (d * 3.5 + maxAcc * 8 + 1) / std::max<double>(dist, 1);
    }
    else
    {
        tvx *= std::max<double>(d, 1.0) / std::max<double>(dist, 1);
        tvy *= std::max<double>(d, 1.0) / std::max<double>(dist, 1);
    }

    // std::cerr << "Tv magnitude: " << sqrt(tvx * tvx + tvy * tvy) << " : " << d << std::endl;

    if (abs(tvx - vx) > abs(tvy - vy))
    {
        if (tvx > vx) return {AccRight, std::min<int>(round(tvx - vx), maxAcc)};
        else return {AccLeft, std::min<int>(round(vx - tvx), maxAcc)};
    }
    if (abs(tvy - vy) > 0)
    {
        if (tvy > vy) return {AccUp, std::min<int>(round(tvy - vy), maxAcc)};
        else return {AccDown, std::min<int>(round(vy - tvy), maxAcc)};
    }
    else return {0, 0};
}

int minX, maxX, minY, maxY;

void findBB()
{
    for (int i = 0; i < g; ++i)
    {
        if (i == 0 || gifts[i].x < minX) minX = gifts[i].x;
        if (i == 0 || gifts[i].x > maxX) maxX = gifts[i].x;
        if (i == 0 || gifts[i].y < minY) minY = gifts[i].y;
        if (i == 0 || gifts[i].y > maxY) maxY = gifts[i].y;
    }

    std::cout << "Bounding box: " << minX << " " << maxX << " and " << minY << " " << maxY << std::endl;
}

void solveOne(int wInit, int cInit)
{
    sol.score = 0;
    sol.acts.clear();
    for (int i = 0; i < g; ++i)
    {
        gifts[i].used = false;
    }

    bool canDeliver = true;
    bool canAcc = true;

    int x = 0;
    int y = 0;
    int vx = 0;
    int vy = 0;
    int time = 0;
    int weight = wInit;
    int carrots = cInit;

    std::vector<int> close = getClose(x, y);
    for (int i : close)
    {
        sol.acts.emplace_back(LoadGift, gifts[i].num);
        sol.acts.emplace_back(DeliverGift, gifts[i].num);
        sol.score += gifts[i].score;
        gifts[i].used = true;
    }

    std::vector<std::pair<int, int>> delivered;
    sol.acts.emplace_back(LoadCarrots, cInit);

    int giftLoadIndex = sol.acts.size();

    while (time < t)
    {
        if ((long long) (x) * x + (long long) (y) * y <= (long long) (d) * d)
        {
            if (cInit > carrots)
            {
                sol.acts.emplace_back(LoadCarrots, cInit - carrots);
            }
            carrots = cInit;
            weight = wInit;
            giftLoadIndex = sol.acts.size();
        }
        if (canDeliver)
        {
            canDeliver = false;
            std::vector<int> close = getClose(x, y);
            for (int i : close)
            {
                if (gifts[i].weight > weight) continue;
                sol.acts.emplace_back(DeliverGift, gifts[i].num);
                sol.score += gifts[i].score;
                weight -= gifts[i].weight;
                delivered.push_back({gifts[i].num, giftLoadIndex});
                gifts[i].used = true;
            }
        }
        else if (canAcc && carrots > 0)
        {
            canAcc = false;
            int maxAcc = getMaxAcc(weight + carrots);

            auto [dir, a] = pickAcc(x, y, vx, vy, weight, maxAcc);
            if (a > 0)
            {
                --carrots;
                sol.acts.emplace_back(dir, a);
                switch (dir)
                {
                case AccUp:
                    vy += a;
                    break;
                case AccDown:
                    vy -= a;
                    break;
                case AccRight:
                    vx += a;
                    break;
                case AccLeft:
                    vx -= a;
                    break;
                }
            }
        }
        else
        {
            canDeliver = true;
            canAcc = true;

            sol.acts.emplace_back(Float, 1);

            x += vx;
            y += vy;
            time += 1;
        }
    }

    for (int i = delivered.size() - 1; i >= 0; --i)
    {
        Action act = {LoadGift, delivered[i].first};
        sol.acts.insert(sol.acts.begin() + delivered[i].second, act);
    }    

    std::cerr << "Score: " << sol.score << " with " << wInit << " " << cInit << std::endl;

    if (sol.score > bestSol.score)
    {
        // std::cerr << "New score: " << sol.score << " with " << wInit << " " << cInit << std::endl;
        bestSol = sol;
    }
}

bool cmpScore(const Gift& left, const Gift& right)
{
    return left.score > right.score;
}

void solve()
{
    std::sort(gifts, gifts + g, cmpScore);

    bestSol.score = -1;

    findBB();

    for (double weight = accRanges[w - 1].maxW; weight >= accRanges[w - 1].maxW / 2; weight /= 1.25)
    {
        for (double carrots = 1; carrots < weight; carrots *= 1.25)
        {
            solveOne(weight - carrots, carrots);
        }
    }

    std::cerr << "Best score: " << bestSol.score << std::endl;
}

int main()
{
    cin >> test;
    inFile.open((tests[test] + ".txt").c_str());
    outFile.open(("sol_" + tests[test] + ".txt").c_str());

    srand(42);

    input();
    solve();
    output();

    return 0;
}
