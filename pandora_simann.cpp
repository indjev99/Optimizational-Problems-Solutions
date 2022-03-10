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
std::ifstream bestFile;

long long bestTotalScore = -1;
long long totalScore;

std::string outputStr;

std::mt19937 generator;
int randNum(int lb, int ub)
{
    std::uniform_int_distribution<int> distribution(lb, ub - 1);
    return distribution(generator);
}

double randReal(double lb, double ub)
{
    std::uniform_real_distribution<double> distribution(lb, ub);
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

void inputBest()
{
    for (int i = 0; i < D; ++i)
    {
        bestFile >> perm[i];
    }
}

bool tryUpdate(int taken)
{
    if (score <= bestScore) return false;

    bestScore = score;
    std::copy(perm, perm + D, bestPerm);

    // std::cerr << "New best score: " << score << "  taken: " << taken << " / (" << maxTime << ", " << D << ")" << std::endl;

    return true;
}

int stamGets[MAX];

bool tryOne()
{
    score = 0;

    int time = 0;
    int stam = stamInit;

    std::fill(stamGets, stamGets + maxTime, 0);
    
    int i;
    for (i = 0; i < D && time < maxTime; ++i, ++time)
    {
        stam += stamGets[time];
        stam = std::min(stam, maxStam);

        const Demon& dem = demons[perm[i]];

        if (dem.stamCost > stam)
        {
            --i;
            continue;
        }

        score += dem.score(time);
        stam -= dem.stamCost;
        stamGets[time + dem.recTime] += dem.recStam;
    }

    return tryUpdate(i);
}

int TRIALS;
int ratio;

double transProb(int newScore, int oldScore, double temp)
{
    if (newScore > oldScore) return 1;

    return exp((newScore - oldScore) / temp);
}

void solve()
{
    tryOne();

    std::cerr << "Initial score: " << bestScore << std::endl;

    double maxTemp = bestScore / ratio;

    for (int k = 0; k < TRIALS; ++k)
    {
        int i, j;

        i = randNum(0, std::min(maxTime, D));

        if (maxTime >= D / 2) j = randNum(0, D);
        else
        {
            if (randNum(0, 2) == 0) j = randNum(0, maxTime);
            else j = randNum(maxTime, D);
        }

        if (i == j)
        {
            --k;
            continue;
        }

        int oldScore = score;

        std::swap(perm[i], perm[j]);
        tryOne();

        int newScore = score;
        double temp = maxTemp * (TRIALS - k) / TRIALS;
        double prob = transProb(newScore, oldScore, temp);
        double sample = randReal(0, 1);

        if (sample > prob)
        {
            std::swap(perm[i], perm[j]);
            score = oldScore;
        }

        if (k % (TRIALS / 100) == 0) std::cerr << "Curr SA score: " << score << " " << (k * 100 + TRIALS / 2) / TRIALS << "%" << std::endl;
    }

    std::copy(perm, perm + D, bestPerm);

    std::cerr << "Final score: " << score << std::endl;
}

int main()
{
    std::cin >> test >> TRIALS >> ratio;
    inFile.open((tests[test] + ".txt").c_str());
    outFile.open(("sa_" + tests[test] + ".txt").c_str());
    bestFile.open(("best_" + tests[test] + ".txt").c_str());

    generator.seed(time(nullptr));

    input();
    inputBest();
    solve();
    output();

    return 0;
}
