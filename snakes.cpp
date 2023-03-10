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
#include <array>
#include <assert.h>

//#define DEBUG

#define FOR_S(i) for (int i = 0; i < s; ++i)
#define FOR_NM(i, j) for (int i = 0; i < n; ++i) for (int j = 0; j < m; ++j)

int test;

const std::string tests[] = {
    "00",
    "01",
    "02",
    "03",
    "04",
    "05",
    "06"
};

std::ifstream inFile;
std::ofstream outFile;

std::mt19937 generator;
int randNum(int lb, int ub)
{
    std::uniform_int_distribution<int> distribution(lb, ub - 1);
    return distribution(generator);
}

const int INF = 1e9;

const int MAX = 5e3;

const int WORMHOLE = -1e6;

struct Point
{
    int x;
    int y;
};

int n, m, s;
int lens[MAX];
int vals[MAX][MAX];
std::vector<Point> wormholes;

void input()
{
    inFile >> m >> n >> s;
    for (int i = 0; i < s; ++i)
    {
        inFile >> lens[i];
    }
    std::string val;
    FOR_NM(i, j)
    {
        inFile >> val;
        if (val == "*")
        {
            vals[i][j] = WORMHOLE;
            wormholes.push_back({i, j});
        }
        else vals[i][j] = std::stoi(val);
    }
}

int score;

int bestScore = -INF;
std::array<std::vector<Point>, MAX> paths;

void output()
{
    FOR_S(i)
    {
        auto& path = paths[i];
        if (path.empty())
        {
            outFile << "\n";
            continue;
        }
        Point prev = {INF, INF};
        bool prevWormhole = false;
        bool start = true;
        int len = 0;
        for (Point p : path)
        {
            if (p.x < 0 || p.y < 0 || p.x >= n || p.y >= m)
            {
                std::cerr << "ERROR: Point outside grid: " << i << " : " << p.y << " " << p.x << std::endl;
                exit(1);
            }

            if (!start) outFile << " ";
            if (start)
            {
                outFile << p.y << " " << p.x;

                if (vals[p.x][p.y] == WORMHOLE)
                {
                    std::cerr << "ERROR: Started on wormhole: " << i << std::endl;
                    exit(1);
                }

                start = false;
            }
            else if (prevWormhole) outFile << p.y << " " << p.x;
            else if (p.x == (prev.x - 1 + n) % n && p.y == prev.y) outFile << "U";
            else if (p.x == (prev.x + 1) % n && p.y == prev.y) outFile << "D";
            else if (p.x == prev.x && p.y == (prev.y - 1 + m) % m) outFile << "L";
            else if (p.x == prev.x && p.y == (prev.y + 1) % m) outFile << "R";
            else
            {
                std::cerr << "ERROR: Invalid next segment: " << i << std::endl;
                exit(1);
            }
            if (vals[p.x][p.y] == WORMHOLE && !prevWormhole) prevWormhole = true;
            else
            {
                prevWormhole = false;
                ++len;
            }
            prev = p;
        }
        if (prevWormhole)
        {
            std::cerr << "ERROR: Ended on wormhole: " << i << std::endl;
            exit(1);
        }
        if (len != lens[i])
        {
            std::cerr << "ERROR: Incorrect length: " << i << " : " << len << " / " << lens[i] << std::endl;
            exit(1);
        }
        outFile << "\n";
    }
}

bool taken[MAX][MAX];

int MAX_TRIALS = 100;

int TOUCH_SCORE = 5;
int UNTOUCH_SCORE = -5;

int touchScore(const Point& curr)
{
    int num = 0;

    Point a;

    a = {(curr.x - 1 + n) % n, curr.y};
    num += taken[a.x][a.y];

    a = {(curr.x + 1) % n, curr.y};
    num += taken[a.x][a.y];

    a = {curr.x, (curr.y - 1 + m) % m};
    num += taken[a.x][a.y];

    a = {curr.x, (curr.y + 1) % m};
    num += taken[a.x][a.y];

    return num * TOUCH_SCORE + (4 - num) * UNTOUCH_SCORE;
}

void solve()
{
    for (Point p : wormholes)
    {
        taken[p.x][p.y] = true;
    }

    int snakesDone = 0;
    score = 0;

    std::vector<int> perm(s);
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(), [](int i, int j){ return lens[i] > lens[j]; });

    for (int i : perm)
    {
        auto& bestPath = paths[i];
        bestPath.clear();

        int bestPathScore = -INF;
        int bestPathScoreExt = -INF;

        for (int trials = 0; trials < MAX_TRIALS; ++trials)
        {
            Point curr;

            do
            {
                curr.x = randNum(0, n);
                curr.y = randNum(0, m);
            }
            while (taken[curr.x][curr.y]);

            std::vector<Point> path;
        
            int len = 0;
            int pathScore = 0;
            int pathScoreExt = 0;

            path.push_back(curr);
            pathScore += vals[curr.x][curr.y];
            pathScoreExt += touchScore(curr) + vals[curr.x][curr.y];
            taken[curr.x][curr.y] = true;
            ++len;

            while (len < lens[i])
            {
                Point next = {INF, INF};
                int nextScore = -INF;
                int nextScoreExt = -INF;

                Point a;

                a = {(curr.x - 1 + n) % n, curr.y};
                if (!taken[a.x][a.y])
                {
                    int valExt = touchScore(a) + vals[a.x][a.y];
                    if (valExt > nextScoreExt)
                    {
                        next = a;
                        nextScore = vals[a.x][a.y];
                        nextScoreExt = valExt;
                    }
                }

                a = {(curr.x + 1) % n, curr.y};
                if (!taken[a.x][a.y])
                {
                    int valExt = touchScore(a) + vals[a.x][a.y];
                    if (valExt > nextScoreExt)
                    {
                        next = a;
                        nextScore = vals[a.x][a.y];
                        nextScoreExt = valExt;
                    }
                }

                a = {curr.x, (curr.y - 1 + m) % m};
                if (!taken[a.x][a.y])
                {
                    int valExt = touchScore(a) + vals[a.x][a.y];
                    if (valExt > nextScoreExt)
                    {
                        next = a;
                        nextScore = vals[a.x][a.y];
                        nextScoreExt = valExt;
                    }
                }

                a = {curr.x, (curr.y + 1) % m};
                if (!taken[a.x][a.y])
                {
                    int valExt = touchScore(a) + vals[a.x][a.y];
                    if (valExt > nextScoreExt)
                    {
                        next = a;
                        nextScore = vals[a.x][a.y];
                        nextScoreExt = valExt;
                    }
                }

                if (nextScoreExt == -INF) break;

                curr = next;

                path.push_back(curr);
                taken[curr.x][curr.y] = true;
                pathScore += nextScore;
                pathScoreExt += nextScoreExt;
                ++len;
            }

            for (Point p : path)
            {
                taken[p.x][p.y] = false;
            }

            if (len == lens[i] && pathScore > 0 && pathScoreExt > bestPathScoreExt)
            {
                bestPath = path;
                bestPathScore = pathScore;
                bestPathScoreExt = pathScoreExt;
            }
        }

        for (Point p : bestPath)
        {
            taken[p.x][p.y] = true;
        }

        if (!bestPath.empty())
        {
            score += bestPathScore;
            ++snakesDone;
        }

        // if (bestPath.empty()) std::cerr << "No path for " << i << " ; len: " << lens[i] << std::endl;
        // else std::cerr << "Path for " << i << " ; len: " << lens[i] << " ; score: " << bestPathScore << " ; scoreExt: " << bestPathScoreExt << std::endl;
    }

    std::cerr << "Score: " << score << " ; Snakes done: " << snakesDone << " / " << s << std::endl;
}

int main()
{
    generator.seed(0);

    std::cin >> test;

    std::cin >> MAX_TRIALS;
    std::cin >> TOUCH_SCORE;
    std::cin >> UNTOUCH_SCORE;

    inFile.open((tests[test] + ".txt").c_str());
    outFile.open(("sol_" + tests[test] + ".txt").c_str());

    generator.seed(0);

    input();
    solve();
    output();

    return 0;
}
