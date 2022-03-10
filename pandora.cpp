#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <map>
#include <stdlib.h>
#include <math.h>
#include <unordered_map>
#include <random>
#include <numeric>
#include <assert.h>

//#define DEBUG

int test;

const std::string tests[] = {
    "00",
    "01",
    "02",
    "03",
    "04",
    "05"
};

const int expTaken[] = {
    4,
    687,
    1489,
    953,
    1339,
    9810
};

std::ifstream inFile;
std::ofstream outFile;
// ifstream bestFile;

long long bestTotalScore = -1;
long long totalScore;

std::string outputStr;

std::mt19937 generator;
int randNum(int lb, int ub)
{
    std::uniform_int_distribution<int> distribution(lb, ub - 1);
    return distribution(generator);
}

const int MAX = 1e6;
const int INF = 1e9;

int stamInit, maxStam, maxTime, D;

struct Demon
{
    int stamCost;
    int recTime;
    int recStam;
    int numFragments;
    std::vector<int> fragments;
    std::vector<int> cumFrags;
    int totalFrags;

    inline int score(int time) const
    {
        int left = maxTime - time - 1;

        if (left >= numFragments) return totalFrags;
        else return cumFrags[left];
    }
};

double avgScore;
Demon demons[MAX];

void input()
{
    inFile >> stamInit >> maxStam >> maxTime >> D;
    
    for (int i = 0; i < D; ++i)
    {
        inFile >> demons[i].stamCost >> demons[i].recTime >> demons[i].recStam >> demons[i].numFragments;

        assert(demons[i].stamCost > 0);
        assert(demons[i].recTime > 0);

        demons[i].fragments.resize(demons[i].numFragments);
        for (int j = 0; j < demons[i].numFragments; ++j)
        {
            inFile >> demons[i].fragments[j];
        }
        demons[i].cumFrags.resize(demons[i].numFragments);
        int prev = 0;
        for (int j = 0; j < demons[i].numFragments; ++j)
        {
            demons[i].cumFrags[j] = prev + demons[i].fragments[j];
            prev = demons[i].cumFrags[j];
        }
        demons[i].totalFrags = prev;

        avgScore += demons[i].totalFrags;
    }

    avgScore /= D;

    std::cerr << "Input done" << std::endl;
}

int score;
int perm[MAX];

int bestScore = -1;
int bestPerm[MAX];

void output()
{
    for (int i = 0; i < D; ++i)
    {
        outFile << bestPerm[i] << "\n";
    }
}


void tryUpdate(int taken, int numLook, int mode, double mult)
{
    if (score <= bestScore) return;

    bestScore = score;
    std::copy(perm, perm + D, bestPerm);

    std::cerr << "New best score: " << score << "  taken: " << taken << " / (" << maxTime << ", " << D << ") with " << numLook << " " << mode << " " << mult << std::endl;
}

int stamGets[MAX];

void solveOne(int numLook, int mode, double mult)
{
    std::shuffle(perm, perm + D, generator);

    score = 0;

    int time = 0;
    int stam = stamInit;

    std::fill(stamGets, stamGets + maxTime, 0);
    
    int i;
    for (i = 0; i < D && time < maxTime; ++i, ++time)
    {
        stam += stamGets[time];
        stam = std::min(stam, maxStam);

        // std::clog << time << " : " << i << " with " << stam << " @ " << score << "\n";

        for (int s = 0; s < numLook && i + s < D; ++s)
        {
            std::swap(perm[i + s], perm[randNum(i + s, D)]);
        }

        int best = -1;
        double bestHeur = 0;

        for (int s = 0; s < numLook && i + s < D; ++s)
        {
            const Demon& dem = demons[perm[i + s]];

            if (dem.stamCost > stam) continue;

            double heur = 0;

            int nextStam = std::min(stam - dem.stamCost + stamGets[time + 1] + (dem.recTime == 1 ? dem.recStam : 0), maxStam);

            int rec = dem.recTime <= maxTime - time - 1 ? dem.recStam : 0;
            rec = std::min(stamGets[time + dem.recTime] + rec, maxStam) - rec;
            if (dem.recTime == 1) rec = 0;

            switch (mode)
            {
            case 0:
            {
                heur = (double) dem.score(time) / avgScore + (double) nextStam / maxStam * mult;
                break;
            }
            
            case 1:
            {
                heur = (double) dem.score(time) / avgScore + (double) (nextStam + (double) rec / dem.recTime) / maxStam * mult;
                break;
            }

            default:
                std::cerr << "No heur mode" << std::endl;
            }
    
            if (time == maxTime - 1) heur = dem.score(time);

            if (best == -1 || heur > bestHeur)
            {
                best = i + s;
                bestHeur = heur;
            }
        }

        if (best == -1)
        {
            --i;
            continue;
        }

        std::swap(perm[i], perm[best]);

        const Demon& dem = demons[perm[i]];

        score += dem.score(time);
        stam -= dem.stamCost;
        stamGets[time + dem.recTime] += dem.recStam;

        assert(stam >= 0);
    }

    tryUpdate(i, numLook, mode, mult);
}

int TRIALS = 300;

void solve()
{
    std::iota(perm, perm + D, 0);

    for (int i = 0; i < TRIALS; ++i)
    {
        int num = randNum(15, 50);
        num *= num;
        solveOne(num, randNum(0, 2), (double) randNum(8, 25) / 4);
    }
}

int main()
{
    std::cin >> test;
    inFile.open((tests[test] + ".txt").c_str());
    outFile.open(("sol_" + tests[test] + ".txt").c_str());
    // bestFile.open(("best_" + tests[test] + ".txt").c_str());

    generator.seed(0);

    input();
    solve();
    output();

    return 0;
}
