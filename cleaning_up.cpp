#include <iostream>
#include <numeric>
#include <algorithm>
#include <random>
#include <chrono>

std::mt19937 generator(0);

int randInt(int lb, int ub)
{
    std::uniform_int_distribution<int> distr(lb, ub);
    return distr(generator);
}

std::chrono::high_resolution_clock::time_point startTChrono, currTChrono;

double timePassed()
{
    using namespace std::chrono;
    currTChrono = high_resolution_clock::now();
    double time = duration_cast<duration<double>>(currTChrono - startTChrono).count();
    return time;
}

const int MAX_N = 100;
const int MAX_DAYS = 500000;

int n, days;
int t[MAX_N];

void input()
{
    std::cin >> n;
    std::cin >> days;

    for (int i = 0; i < n; i++)
    {
        std::cin >> t[i];
    }
}

int bestScore;
int bestOddNext[MAX_N];
int bestEvenNext[MAX_N];

void output()
{
    for (int i = 0; i < n; i++)
    {
        std::cout << bestOddNext[i] << " " << bestEvenNext[i] << "\n";
    }
}

int oddNext[MAX_N];
int evenNext[MAX_N];

int currT[MAX_N];

int perm[MAX_N];

int noise[MAX_N];

void solveOne(int noiseLevel = 0)
{
    std::fill(oddNext, oddNext + n, 0);
    std::fill(evenNext, evenNext + n, 0);

    std::fill(currT, currT + n, 0);
    currT[0] = 1;

    for (int i = 0; i < n; i++)
    {
        noise[i] = randInt(0, noiseLevel);
    }

    std::iota(perm, perm + n, 0);
    std::sort(perm, perm + n, [](int i, int j)
    {
        return t[i] + noise[i] > t[j] + noise[j];
    });

    for (int ii = 0; ii < n; ii++)
    {
        int i = perm[ii];

        for (int rep = 1; rep >= 0; rep--)
        {
            int toGive = (t[i] + rep) / 2;

            double maxScore = -1e6;
            int bestJ = -1;

            for (int jj = 0; jj < n; jj++)
            {
                int j = perm[jj];

                if (j == i) continue;

                double score = t[j] - currT[j];

                if (score - toGive >= 0)
                {
                    double ratio = (double) (score - toGive) / t[j];

                    if (ratio <= 0.0033) score *= 7;
                    else if (ratio <= 0.0067) score *= 5;
                    else if (ratio <= 0.01) score *= 3;

                    // if (ratio <= 0.01)
                    // {
                    //     score *= 3 + (0.01 - ratio) * 550;
                    // }

                    score *= std::pow(t[j], 0.6);
                }

                if (score > maxScore)
                {
                    maxScore = score;
                    bestJ = j;
                }            
            }

            currT[bestJ] += toGive;
            if (rep == 1) oddNext[i] = bestJ;
            else evenNext[i] = bestJ;
        }
    }
}

int actTs[MAX_N];

void simulate()
{
    std::fill(actTs, actTs + n, 0);

    int curr = 0;

    for (int i = 0; i < days; i++)
    {
        actTs[curr]++;

        if (actTs[curr] % 2 == 1) curr = oddNext[curr];
        else curr = evenNext[curr];
    }
}

int score;

void calcScore()
{
    score = 0;
    for (int i = 0; i < n; i++)
    {
        score += std::abs(actTs[i] - t[i]);
    }
}

void maybeUpdate()
{
    simulate();
    calcScore();

    if (score < bestScore)
    {
        // std::cerr << "New best score: " << score << std::endl;
    
        bestScore = score;
        std::copy(evenNext, evenNext + n, bestEvenNext);
        std::copy(oddNext, oddNext + n, bestOddNext);
    }
}

int prevActTs[MAX_N];

void optimize(bool persist = false)
{
    int bestDScore = -1e6;

    int bestX = 0;
    int bestY = 0;
    int bestSx = 0;
    int bestSy = 0;

    for (int trials = 0; trials < 2000; trials++)
    {
        int x = randInt(0, n - 1);
        int y = randInt(0, n - 1);

        if (x == y) continue;

        int sx = randInt(0, 1);
        int sy = randInt(0, 1);

        int fx = sx ? oddNext[x] : evenNext[x];
        int fy = sy ? oddNext[y] : evenNext[y];

        if (fx == x || fy == y || fx == y || fy == x || fx == fy) continue;

        int dx = (actTs[x] + sx) / 2;
        int dy = (actTs[y] + sy) / 2;

        int oldScore = std::abs(t[fx] - actTs[fx]) + std::abs(t[fy] - actTs[fy]);
        int newScore = std::abs(t[fx] - (actTs[fx] - dx + dy)) + std::abs(t[fy] - (actTs[fy] - dy + dx));

        int dScore = oldScore - newScore;

        if (dScore > bestDScore)
        {
            bestDScore = dScore;

            bestX = x;
            bestY = y;
            bestSx = sx;
            bestSy = sy;
        }
    }

    if (bestDScore <= 0) return;

    int& fx = bestSx ? oddNext[bestX] : evenNext[bestX];
    int& fy = bestSy ? oddNext[bestY] : evenNext[bestY];

    int oldScore = score;
    std::copy(actTs, actTs + n, prevActTs);

    std::swap(fx, fy);

    maybeUpdate();

    if (score > oldScore && !persist)
    {
        std::swap(fx, fy);
        std::copy(prevActTs, prevActTs + n, actTs);
        score = oldScore;
    }
}

void solve()
{
    bestScore = 1e6;

    solveOne(0);
    maybeUpdate();

    // std::cerr << "Trying noise:" << std::endl;

    // while (timePassed() < 0.1)
    // {
    //     int noiseLevel = randInt(50, 200);
    //     solveOne(noiseLevel);
    //     maybeUpdate();
    // }

    // std::cerr << "Local optimization:" << std::endl;

    score = bestScore;
    std::copy(bestEvenNext, bestEvenNext + n, evenNext);
    std::copy(bestOddNext, bestOddNext + n, oddNext);

    while (timePassed() < 1.0)
    {
        optimize(true);
    }

    score = bestScore;
    std::copy(bestEvenNext, bestEvenNext + n, evenNext);
    std::copy(bestOddNext, bestOddNext + n, oddNext);
    maybeUpdate();

    while (timePassed() < 1.95)
    {
        optimize();
    }

    std::cerr << bestScore << std::endl;
}

int main()
{
    startTChrono = std::chrono::high_resolution_clock::now();

    input();
    solve();
    output();

    return 0;
}
