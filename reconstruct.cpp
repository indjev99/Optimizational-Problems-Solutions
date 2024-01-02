#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <random>
#include <cmath>
#include <string>
#include <array>

#define _USE_MATH_DEFINES

std::mt19937 generator;
int randInt(int lb, int rb)
{
    std::uniform_int_distribution<int> distribution(lb, rb - 1);
    return distribution(generator);
}

constexpr int MAX_N = 512;

constexpr int KBUCK_SIZE = 16;
constexpr int MAX_NUM_KBUCKS = MAX_N / KBUCK_SIZE;

int n = MAX_N;
int k;

int NUM_BUCKS;
int MAX_KBUCK;

bool ans[MAX_N][MAX_N];

bool usingTrueAns;
bool trueAns[MAX_N][MAX_N];

void inputTrue()
{
    usingTrueAns = true;

    std::cin >> n;
    std::cin >> k;

    for (int i = 0; i < n; ++i)
    {
        std::string row;
        std::cin >> row;
        for (int j = 0; j < n; ++j)
        {
            trueAns[i][j] = row[j] == 'X';
        }
    }
}

void input()
{
    std::string cmd;
    std::cin >> cmd;

    if (cmd == "true")
    {
        inputTrue();
    }
    else
    {
        k = std::stoi(cmd);
    }

    MAX_KBUCK = (n - 1) / KBUCK_SIZE;
    NUM_BUCKS = MAX_KBUCK + 1;
}

bool isKnown[MAX_N][MAX_N];

void knownOutput()
{
    if (!usingTrueAns) return;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (!isKnown[i][j]) std::cout << " ";
            else if (ans[i][j]) std::cout << "X";
            else std::cout << ".";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

}

void trueOutput()
{
    knownOutput();

    int cnt = 0;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (ans[i][j] == trueAns[i][j]) ++cnt;
            if (ans[i][j]) std::cout << "X";
            else std::cout << ".";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cerr << "Correct: " << cnt << std::endl;
}

void output()
{
    if (usingTrueAns)
    {
        trueOutput();
        return;
    }

    std::cout << "Ready" << std::endl;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            std::cout << ans[i][j];
        }
        std::cout << std::endl;
    }
}

bool haveDrawn = false;

constexpr int NO_GUESS = 0;
constexpr int RIGHT_GUESS = 1;
constexpr int WRONG_GUESS = 2;

int guesses[MAX_N][MAX_N];
int totalGuessCnt[3];

std::vector<std::pair<int, int>> knownByBuck[MAX_NUM_KBUCKS][MAX_NUM_KBUCKS][2];

void setKnown(int x, int y, bool val)
{
    if (isKnown[x][y]) return;

    if (haveDrawn)
    {
        guesses[x][y] = val == ans[x][y] ? RIGHT_GUESS : WRONG_GUESS;
        ++totalGuessCnt[guesses[x][y]];
    }

    isKnown[x][y] = true;
    ans[x][y] = val;

    knownByBuck[x / KBUCK_SIZE][y / KBUCK_SIZE][ans[x][y]].emplace_back(x, y);
}

void setAllKnown(int x0, int y0, int x1, int y1, int xS, int yS)
{
    int dx = std::abs(x1 - x0);
    int sx = x0 < x1 ? 1 : -1;
    int dy = - std::abs(y1 - y0);
    int sy = y0 < y1 ? 1 : -1;
    int error = dx + dy;

    if (xS != -1 && yS != -1)
    {
        setKnown(xS, yS, true);
    }

    while (true)
    {
        if (x0 == xS && y0 == yS) break;

        setKnown(x0, y0, false);

        if (x0 == x1 && y0 == y1) break;

        int e2 = 2 * error;
        if (e2 >= dy)
        {
            if (x0 == x1) break;
            error = error + dy;
            x0 = x0 + sx;
        }
        if (e2 <= dx)
        {
            if (y0 == y1) break;
            error = error + dx;
            y0 = y0 + sy;
        }
    }
}

constexpr int INF = 1e6;

std::vector<std::pair<int, int>> getBucks(int x, int y)
{
    std::vector<std::pair<int, int>> bucks;

    int xB = x / KBUCK_SIZE;
    int yB = y / KBUCK_SIZE;

    bool down = x % KBUCK_SIZE >= KBUCK_SIZE / 2;
    bool right = y % KBUCK_SIZE >= KBUCK_SIZE / 2;

    bool up = !down;
    bool left = !right;

    if (knownByBuck[xB][yB][0].size() + knownByBuck[xB][yB][1].size() == 0)
    {
        down = right = up = left = true;
    }

    bucks.emplace_back(xB, yB);
    if (up && xB > 0) bucks.emplace_back(xB - 1, yB);
    if (down && xB < MAX_KBUCK) bucks.emplace_back(xB + 1, yB);
    if (left && yB > 0) bucks.emplace_back(xB, yB - 1);
    if (right && yB < MAX_KBUCK) bucks.emplace_back(xB, yB + 1);
    if (left && up && xB > 0 && yB > 0) bucks.emplace_back(xB - 1, yB - 1);
    if (left && down && xB < MAX_KBUCK && yB > 0) bucks.emplace_back(xB + 1, yB - 1);
    if (right && up && xB > 0 && yB < MAX_KBUCK) bucks.emplace_back(xB - 1, yB + 1);
    if (right && down && xB < MAX_KBUCK && yB < MAX_KBUCK) bucks.emplace_back(xB + 1, yB + 1);

    return bucks;
}

bool closestKnown(int x, int y)
{
    if (isKnown[x][y]) return ans[x][y];

    std::vector<std::pair<int, int>> bucks = getBucks(x, y);

    int cntByCol[2] = {0, 0};

    for (auto [xB2, yB2] : bucks)
    {
        for (int col = 0; col < 2; ++col)
        {
            cntByCol[col] += knownByBuck[xB2][yB2][col].size();
        }
    }

    if (cntByCol[0] == 0 || cntByCol[1] == 0)
    {
        return cntByCol[1] > 0;
    }

    int minDist2ByCol[2];

    minDist2ByCol[0] = KBUCK_SIZE * KBUCK_SIZE * 8;
    minDist2ByCol[1] = KBUCK_SIZE * KBUCK_SIZE * 8;

    for (auto [xB2, yB2] : bucks)
    {
        for (int col = 0; col < 2; ++col)
        {
            for (auto [x2, y2] : knownByBuck[xB2][yB2][col])
            {
                int dist2 = (x2 - x) * (x2 - x) + (y2 - y) * (y2 - y);

                if (dist2 < minDist2ByCol[col])
                {
                    minDist2ByCol[col] = dist2;
                }
            }
        }
    }

    return minDist2ByCol[1] < minDist2ByCol[0];
}

void makeTrueQuery(int x0, int y0, int x1, int y1, int& xS, int& yS)
{
    std::cout << "> " << x0 << " " << y0 << " " << x1 << " " << y1 << std::endl;

    xS = -1;
    yS = -1;

    int dx = std::abs(x1 - x0);
    int sx = x0 < x1 ? 1 : -1;
    int dy = - std::abs(y1 - y0);
    int sy = y0 < y1 ? 1 : -1;
    int error = dx + dy;

    while (true)
    {
        if (trueAns[x0][y0])
        {
            xS = x0;
            yS = y0;
            break;
        }

        if (x0 == x1 && y0 == y1) break;

        int e2 = 2 * error;
        if (e2 >= dy)
        {
            if (x0 == x1) break;
            error = error + dy;
            x0 = x0 + sx;
        }
        if (e2 <= dx)
        {
            if (y0 == y1) break;
            error = error + dx;
            y0 = y0 + sy;
        }
    }

    std::cout << "< " << xS << " " << yS << std::endl;
    std::cout << std::endl;
}

void makeQuery(int x0, int y0, int x1, int y1)
{
    int xS, yS;

    if (!usingTrueAns)
    {
        std::cout << x0 << " " << y0 << " " << x1 << " " << y1 << std::endl;
        std::cin >> xS >> yS; 
    }
    else
    {
        makeTrueQuery(x0, y0, x1, y1, xS, yS);
    }

    setAllKnown(x0, y0, x1, y1, xS, yS);
}

int currQuery;

double getQScore(int x, int y)
{
    std::vector<std::pair<int, int>> bucks = getBucks(x, y);

    int cntByCol[2] = {0, 0};

    for (auto [xB2, yB2] : bucks)
    {
        for (int col = 0; col < 2; ++col)
        {
            cntByCol[col] += knownByBuck[xB2][yB2][col].size();
        }
    }

    if (cntByCol[0] == 0 && cntByCol[1] == 0)
    {
        return INF;
    }

    if (cntByCol[1] == 0)
    {
        return 0;
    }

    // if (cntByCol[0] == 0 || cntByCol[1] == 0)
    // {
    //     return 0;
    // }

    int minDists2[2] = {INF, INF};

    for (auto [xB2, yB2] : bucks)
    {
        for (int col = 0; col < 2; ++col)
        {
            for (auto [x2, y2] : knownByBuck[xB2][yB2][col])
            {
                int dist2 = (x2 - x) * (x2 - x) + (y2 - y) * (y2 - y);
                if (dist2 < minDists2[col])
                {
                    minDists2[col] = dist2;
                }
            }
        }
    }

    if (minDists2[0] > minDists2[1] * 4) return 0;

    double min = std::min(minDists2[0], minDists2[1]);
    double max = std::max(minDists2[0], minDists2[1]);

    if (min == INF) return INF;

    double temp = 0.15 * (k - currQuery) / k;

    double whiteBoost = 1.1;

    return pow(min, temp) * sqrt(min / max) * (minDists2[0] < minDists2[1] ? whiteBoost : 1);
}

constexpr int LINE_QSCORE_LEN = 10;
constexpr int LINE_QSCORE_FREQ = 1;

double getLineQScore(int x0, int y0, int x1, int y1)
{
    double score = 0;

    int dx = std::abs(x1 - x0);
    int sx = x0 < x1 ? 1 : -1;
    int dy = - std::abs(y1 - y0);
    int sy = y0 < y1 ? 1 : -1;
    int error = dx + dy;

    int len = LINE_QSCORE_LEN;

    while (len > 0)
    {
        --len;

        if (x0 == x1 && y0 == y1) break;

        int e2 = 2 * error;
        if (e2 >= dy)
        {
            if (x0 == x1) break;
            error = error + dy;
            x0 = x0 + sx;
        }
        if (e2 <= dx)
        {
            if (y0 == y1) break;
            error = error + dx;
            y0 = y0 + sy;
        }

        if (len % LINE_QSCORE_FREQ == 0)
        {
            score += getQScore(x0, y0);
        }
    }

    return score;
}

void makeRandQuery()
{
    int x0, y0, x1, y1;

    x0 = randInt(0, n);
    y0 = randInt(0, n);

    if (currQuery > k / 10)
    {
        int maxRetries = 500;

        double bestScore = 0;

        for (int i = 0; i < maxRetries; ++i)
        {
            int x = randInt(0, n);
            int y = randInt(0, n);

            if (isKnown[x][y])
            {
                --maxRetries;
                continue;
            }

            double score = getQScore(x, y);

            if (score > bestScore)
            {
                bestScore = score;
                x0 = x;
                y0 = y;
            }
        }
    }

    if (randInt(0, 2) == 0)
    {
        x1 = x0;
        y1 = randInt(0, 2) == 0 ? 0 : n - 1;
    }
    else
    {
        x1 = randInt(0, 2) == 0 ? 0 : n - 1;
        y1 = y0;
    }

    if (currQuery > k / 10)
    {
        int maxRetries2 = 10;

        int bestScore2 = 0;
        
        for (int i = 0; i < maxRetries2; ++i)
        {
            int x;
            int y;

            if (randInt(0, 2) == 0)
            {
                x = x0;
                y = randInt(0, 2) == 0 ? 0 : n - 1;
            }
            else
            {
                x = randInt(0, 2) == 0 ? 0 : n - 1;
                y = y0;
            }

            double score = getLineQScore(x0, y0, x, y);

            if (score > bestScore2)
            {
                bestScore2 = score;
                x1 = x;
                y1 = y;
            }
        }
    }

    // double angle = randInt(0, 10000) / 10000.0 * 2 * M_PI;
    // double i1 = (511 - x0) / sin(angle);
    // double i2 = - x0 / sin(angle);
    // double i3 = (511 - y0) / cos(angle);
    // double i4 = - y0 / cos(angle);
    // if (i1 < 0) i1 = INF;
    // if (i2 < 0) i2 = INF;
    // if (i3 < 0) i3 = INF;
    // if (i4 < 0) i4 = INF;
    // double iMin = std::min(std::min(i1, i2), std::min(i3, i4));
    // if (i1 == iMin)
    // {
    //     x1 = 511;
    //     y1 = std::round(y0 + cos(angle) * i1);
    // }
    // else if (i2 == iMin)
    // {
    //     x1 = 0;
    //     y1 = std::round(y0 + cos(angle) * i2);
    // }
    // else if (i3 == iMin)
    // {
    //     y1 = 511;
    //     x1 = std::round(x0 + sin(angle) * i3);
    // }
    // else if (i4 == iMin)
    // {
    //     y1 = 0;
    //     x1 = std::round(x0 + sin(angle) * i4);
    // }

    makeQuery(x0, y0, x1, y1);

    // knownOutput();
}

constexpr int KNOWN_WEIGHT = 2;
constexpr int BLUR_RAD = 3;
constexpr double BLUR_STD = 4;

constexpr int MAX_BLUR_DIST_2 = 2 * BLUR_RAD * BLUR_RAD;

constexpr std::array<double, MAX_BLUR_DIST_2 + 1> computeBlurWeights()
{
    std::array<double, MAX_BLUR_DIST_2 + 1> weights = {};
    for (int d2 = 0; d2 <= MAX_BLUR_DIST_2; ++d2)
    {
        weights[d2] = std::exp(- d2 / (BLUR_STD * BLUR_STD));
    }
    return weights;
}

constexpr std::array<double, MAX_BLUR_DIST_2 + 1> blurWeights = computeBlurWeights();

void cleanup(int x, int y)
{
    if (isKnown[x][y]) return;

    int cnts[2] = {0, 0};
    double farCnts[2] = {0, 0};

    if (x > 0) cnts[ans[x - 1][y]] += isKnown[x - 1][y] ? KNOWN_WEIGHT : 1;
    if (x < n - 1) cnts[ans[x + 1][y]] += isKnown[x + 1][y] ? KNOWN_WEIGHT : 1;
    if (y > 0) cnts[ans[x][y - 1]] += isKnown[x][y - 1] ? KNOWN_WEIGHT : 1;
    if (y < n - 1) cnts[ans[x][y + 1]] += isKnown[x][y + 1] ? KNOWN_WEIGHT : 1;

    if (cnts[0] != cnts[1])
    {
        ans[x][y] = cnts[1] > cnts[0];
        return;
    }

    constexpr int MAX_LEN = 30;

    for (int dir = 0; dir < 4; ++dir)
    {
        int cnt = 0;

        switch (dir)
        {
        case 0:
            if (x < n - 1) cnt += ans[x + 1][y] == ans[x][y] ? 1 : -1;
            break;
        case 1:
            if (x > 0) cnt += ans[x - 1][y] == ans[x][y] ? 1 : -1;
            break;
        case 2:
            if (y < n - 1) cnt += ans[x][y + 1] == ans[x][y] ? 1 : -1;
            break;
        case 3:
            if (y > 0) cnt += ans[x][y - 1] == ans[x][y] ? 1 : -1;
            break;
        }

        if (cnt > 0) continue;

        for (int i = 0; i < MAX_LEN; ++i)
        {
            int x2 = x;
            int y2 = y;
            switch (dir)
            {
            case 0:
                x2 -= i;
                break;
            case 1:
                x2 += i;
                break;
            case 2:
                y2 -= i;
                break;
            case 3:
                y2 += i;
                break;
            }

            if (x2 < 0 || x2 >= n || y2 < 0 || y2 >= n) break;

            if (ans[x2][y2] != ans[x][y])
            {
                --cnt;
                break;
            }

            if (isKnown[x2][y2])
            {
                cnt = 100;
                break;
            }

            if (dir >= 2 && x2 > 0) cnt += ans[x2 - 1][y2] == ans[x][y] ? 1 : -1;
            if (dir >= 2 && x2 < n - 1) cnt += ans[x2 + 1][y2] == ans[x][y] ? 1 : -1;
            if (dir <= 1 && y2 > 0) cnt += ans[x2][y2 - 1] == ans[x][y] ? 1 : -1;
            if (dir <= 1 && y2 < n - 1) cnt += ans[x2][y2 + 1] == ans[x][y] ? 1 : -1;
        }

        if (cnt >= -1) continue;

        ans[x][y] = !ans[x][y];

        for (int i = 1; i < MAX_LEN; ++i)
        {
            int x2 = x;
            int y2 = y;
            switch (dir)
            {
            case 0:
                x2 -= i;
                break;
            case 1:
                x2 += i;
                break;
            case 2:
                y2 -= i;
                break;
            case 3:
                y2 += i;
                break;
            }

            if (x2 < 0 || x2 >= n || y2 < 0 || y2 >= n) break;

            if (ans[x2][y2] == ans[x][y]) break;

            ans[x2][y2] = ans[x][y];
        }
    }

    for (int i = x - BLUR_RAD; i <= x + BLUR_RAD; ++i)
    {
        for (int j = y - BLUR_RAD; j <= y + BLUR_RAD; ++j)
        {
            if (i < 0 || i >= n || j < 0 || j >= n) continue;
            if (i == x && j == y) continue;
            int dist2 = (i - x) * (i - x) + (j - y) * (j - y);
            if (dist2 > (BLUR_RAD + 0.5) * (BLUR_RAD + 0.5)) continue;
            farCnts[ans[i][j]] += blurWeights[dist2];
        }
    }

    if (farCnts[0] != farCnts[1])
    {
        ans[x][y] = farCnts[1] > farCnts[0];
        return;
    }
}

void draw()
{
    haveDrawn = true;

    std::vector<std::pair<int, int>> points;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            ans[i][j] = closestKnown(i, j);
            points.emplace_back(i, j);
        }
    }

    for (int t = 0; t < 10; ++t)
    {
        for (auto [i, j] : points)
        {
            if (t / 2 % 2 == 1) std::swap(i, j);
            cleanup(i, j);
        }
        std::reverse(points.begin(), points.end());
    }
}

void solve()
{
    for (; currQuery < k; ++currQuery)
    {
        makeRandQuery();
    }

    draw();
}

int main()
{
    std::ios::sync_with_stdio(false);
    generator.seed(1337);

    input();
    solve();
    output();

    return 0;
}
