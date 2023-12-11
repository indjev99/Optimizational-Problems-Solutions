#pragma GCC optimize("Ofast")

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
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <iterator>
#include <assert.h>

#define FORCE_INLINE __attribute__((always_inline)) inline

constexpr bool DEBUG = true;
constexpr double TIME_MULT = 1;

const double SOFT_MAX_TIME = 4.0;
const double HARD_MAX_TIME = 4.5;

double CURR_MAX_TIME;

bool DEBUG_PRINT = false;

std::ifstream inF("gerrymandering.in");
std::ofstream outF("gerrymandering.out");

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

const int MAX_N = 512;
const int MAX_D = 512;
const int MAX_P = 8;

int n, d, p;
int votes[MAX_N][MAX_N];

void input()
{
    inF >> n >> d >> p;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            inF >> votes[i][j];
            --votes[i][j];
        }
    }
}

int regions[MAX_N][MAX_N];

void output()
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            outF << regions[i][j] + 1 << " ";
        }
        outF << "\n";
    }
}

std::map<int, std::pair<int, int>> sides = {
    {3600, {12, 300}},
    {1800, {6, 300}},
    {750, {5, 150}},
    {375, {5, 75}},
    {180, {6, 30}},
    {25, {5, 5}},
    {5, {1, 5}}};

int type;
int counts[MAX_P];
double partyScores[MAX_P + 1];

void myExit(int code)
{
    std::cerr << "Code: " << code << std::endl;
    exit(code);
}

void setup()
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            ++counts[votes[i][j]];
        }
    }

    partyScores[p] = 0;
    for (int k = 0; k < p; ++k)
    {
        partyScores[k] = 1.0 * n * n / counts[k];
    }

    if (counts[0] - counts[1] > 0.2 * n * n)
        type = 1;
    else if (counts[0] > 0.4 * n * n && counts[1] > 0.4 * n * n)
        type = 2;
    else
        type = 3;

    std::cerr << "D: " << d << " Type: " << type << std::endl;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            regions[i][j] = -1;
        }
    }
}

struct Cell
{
    int idx;

    FORCE_INLINE
    Cell() : idx(0) {}

    FORCE_INLINE
    Cell(int idx) : idx(idx) {}

    FORCE_INLINE
    Cell(int x, int y) : idx(x * MAX_N + y) {}

    FORCE_INLINE
    int x()
    {
        return idx / MAX_N;
    }

    FORCE_INLINE
    int y()
    {
        return idx % MAX_N;
    }
};

FORCE_INLINE
bool willTouchRegsAfterSwap(Cell a, Cell b)
{
    if (std::abs(a.x() - b.x()) + std::abs(a.y() - b.y()) > 1) return true;

    int aReg = regions[a.x()][a.y()];
    int bReg = regions[b.x()][b.y()];

    int aAdj = 0;
    aAdj += a.x() > 0 && regions[a.x() - 1][a.y()] == bReg;
    aAdj += a.x() < n - 1 && regions[a.x() + 1][a.y()] == bReg;
    aAdj += a.y() > 0 && regions[a.x()][a.y() - 1] == bReg;
    aAdj += a.y() < n - 1 && regions[a.x()][a.y() + 1] == bReg;

    if (aAdj <= 1) return false;

    int bAdj = 0;
    bAdj += b.x() > 0 && regions[b.x() - 1][b.y()] == aReg;
    bAdj += b.x() < n - 1 && regions[b.x() + 1][b.y()] == aReg;
    bAdj += b.y() > 0 && regions[b.x()][b.y() - 1] == aReg;
    bAdj += b.y() < n - 1 && regions[b.x()][b.y() + 1] == aReg;

    if (bAdj <= 1) return false;

    return true;
}

FORCE_INLINE
bool operator==(Cell left, Cell right)
{
    return left.idx == right.idx;
}

namespace std
{
    template <>
    struct hash<Cell>
    {
        FORCE_INLINE
        std::size_t operator()(Cell c) const
        {
            return std::hash<int>()(c.idx);
        }
    };
}

const int P_THICKNESS = 3;

void makePortRegions()
{
    int regSize = n * n / d;

    int regH = sides[regSize].first;
    int regW = sides[regSize].second;
    int divsW = n / regW;

    for (int r = 0; r < d; ++r)
    {
        int rX = r / divsW;
        int rY = r % divsW;

        int begX = rX * regH;
        int begY = rY * regW;

        for (int i = 0; i < regH; ++i)
        {
            for (int j = 0; j < regW; ++j)
            {
                regions[begX + i][begY + j] = r;
            }
        }
    }

    for (int r = 0; r < d; ++r)
    {
        int rX = r / divsW;
        int rY = r % divsW;

        int begX = rX * regH;
        int begY = rY * regW;

        if (rX > 0)
        {
            for (int y = (regW % (P_THICKNESS * 2)) / 2; y + P_THICKNESS * 2 < regW; y += P_THICKNESS * 2)
            {
                int maxD = (regH - 2) / 2;
                for (int x = 0; x < maxD; ++x)
                {
                    for (int y2 = y; y2 < y + P_THICKNESS; ++y2)
                    {
                        Cell a(begX + x, begY + y2);
                        Cell b(begX - 1 - x, begY + y2 + P_THICKNESS);
                        std::swap(regions[a.x()][a.y()], regions[b.x()][b.y()]);
                    }
                }
            }
        }
    }
}

const int SNAKESIZE = 2;

std::vector<Cell> snake;

void makeSnake()
{
    std::vector<Cell> snakeBack;

    assert(n % (2 * SNAKESIZE) == 0);

    for (int i = 0; i < n; i += 2 * SNAKESIZE)
    {
        for (int i2 = 0; i2 < SNAKESIZE; ++i2)
        {
            snakeBack.emplace_back(i + i2, n - 1);
        }
        for (int i2 = 0; i2 < SNAKESIZE; ++i2)
        {
            for (int j = 0; j < n - 1; ++j)
            {
                int actJ = i2 % 2 == 0 ? j : n - j - 2;
                snake.emplace_back(i + i2, actJ);
            }
            for (int j = n - 1; j > 0; --j)
            {
                int actJ = i2 % 2 == 0 ? j : n - j;
                snakeBack.emplace_back(i + SNAKESIZE + i2, actJ);
            }
        }
        for (int i2 = 0; i2 < SNAKESIZE; ++i2)
        {
            snake.emplace_back(i + SNAKESIZE + i2, 0);
        }
    }

    for (int i = snakeBack.size() - 1; i >= 0; --i)
    {
        snake.push_back(snakeBack[i]);
    }
}

void makeSnakeRegions()
{
    int regSize = n * n / d;

    int currReg = 0;
    int currSize = 0;
    for (Cell cell : snake)
    {
        regions[cell.x()][cell.y()] = currReg;
        if (++currSize == regSize)
        {
            ++currReg;
            currSize = 0;
        }
    }
}

// cells in party P, which can be moved from region R1 to region R2
std::unordered_set<Cell> mCells[MAX_P][MAX_D][MAX_D];

FORCE_INLINE
std::vector<int> findOtherRegs(Cell c)
{
    std::vector<int> otherRegs;

    int reg = regions[c.x()][c.y()];

    if (c.x() > 0)
    {
        int otherReg = regions[c.x() - 1][c.y()];
        if (otherReg != reg)
            otherRegs.push_back(otherReg);
    }

    if (c.x() < n - 1)
    {
        int otherReg = regions[c.x() + 1][c.y()];
        if (otherReg != reg)
            otherRegs.push_back(otherReg);
    }

    if (c.y() > 0)
    {
        int otherReg = regions[c.x()][c.y() - 1];
        if (otherReg != reg)
            otherRegs.push_back(otherReg);
    }

    if (c.y() < n - 1)
    {
        int otherReg = regions[c.x()][c.y() + 1];
        if (otherReg != reg)
            otherRegs.push_back(otherReg);
    }

    return otherRegs;
}

FORCE_INLINE
bool isMCell(Cell c)
{
    int reg = regions[c.x()][c.y()];

    bool uSame = c.x() > 0 && regions[c.x() - 1][c.y()] == reg;
    bool dSame = c.x() < n - 1 && regions[c.x() + 1][c.y()] == reg;
    bool lSame = c.y() > 0 && regions[c.x()][c.y() - 1] == reg;
    bool rSame = c.y() < n - 1 && regions[c.x()][c.y() + 1] == reg;

    if (std::abs((uSame + dSame) - (lSame + rSame)) == 2)
        return false;
    if (uSame && lSame && regions[c.x() - 1][c.y() - 1] != reg)
        return false;
    if (uSame && rSame && regions[c.x() - 1][c.y() + 1] != reg)
        return false;
    if (dSame && lSame && regions[c.x() + 1][c.y() - 1] != reg)
        return false;
    if (dSame && rSame && regions[c.x() + 1][c.y() + 1] != reg)
        return false;

    bool hasOthers = (c.x() > 0 && !uSame) || (c.x() < n - 1 && !dSame) || (c.y() > 0 && !lSame) || (c.y() < n - 1 && !rSame);

    return hasOthers;
}

FORCE_INLINE
void maybeRemMCell(Cell c)
{
    if (!isMCell(c))
        return;

    int reg = regions[c.x()][c.y()];
    int party = votes[c.x()][c.y()];
    std::vector<int> otherRegs = findOtherRegs(c);

    for (int otherReg : otherRegs)
    {
        mCells[party][reg][otherReg].erase(c);
    }
}

FORCE_INLINE
void maybeAddMCell(Cell c)
{
    if (!isMCell(c))
        return;

    int reg = regions[c.x()][c.y()];
    int party = votes[c.x()][c.y()];
    std::vector<int> otherRegs = findOtherRegs(c);

    for (int otherReg : otherRegs)
    {
        mCells[party][reg][otherReg].insert(c);
    }
}

FORCE_INLINE
void maybeRemMCells3x3(Cell c)
{
    maybeRemMCell(c);
    if (c.x() > 0)
        maybeRemMCell(Cell(c.x() - 1, c.y()));
    if (c.x() < n - 1)
        maybeRemMCell(Cell(c.x() + 1, c.y()));
    if (c.y() > 0)
        maybeRemMCell(Cell(c.x(), c.y() - 1));
    if (c.y() < n - 1)
        maybeRemMCell(Cell(c.x(), c.y() + 1));
    if (c.x() > 0 && c.y() > 0)
        maybeRemMCell(Cell(c.x() - 1, c.y() - 1));
    if (c.x() > 0 && c.y() < n - 1)
        maybeRemMCell(Cell(c.x() - 1, c.y() + 1));
    if (c.x() < n - 1 && c.y() > 0)
        maybeRemMCell(Cell(c.x() + 1, c.y() - 1));
    if (c.x() < n - 1 && c.y() < n - 1)
        maybeRemMCell(Cell(c.x() + 1, c.y() + 1));
}

FORCE_INLINE
void maybeAddMCells3x3(Cell c)
{
    maybeAddMCell(c);
    if (c.x() > 0)
        maybeAddMCell(Cell(c.x() - 1, c.y()));
    if (c.x() < n - 1)
        maybeAddMCell(Cell(c.x() + 1, c.y()));
    if (c.y() > 0)
        maybeAddMCell(Cell(c.x(), c.y() - 1));
    if (c.y() < n - 1)
        maybeAddMCell(Cell(c.x(), c.y() + 1));
    if (c.x() > 0 && c.y() > 0)
        maybeAddMCell(Cell(c.x() - 1, c.y() - 1));
    if (c.x() > 0 && c.y() < n - 1)
        maybeAddMCell(Cell(c.x() - 1, c.y() + 1));
    if (c.x() < n - 1 && c.y() > 0)
        maybeAddMCell(Cell(c.x() + 1, c.y() - 1));
    if (c.x() < n - 1 && c.y() < n - 1)
        maybeAddMCell(Cell(c.x() + 1, c.y() + 1));
}

void findMCells()
{
    for (int k = 0; k < p; ++k)
    {
        for (int r = 0; r < d; ++r)
        {
            for (int r2 = 0; r2 < d; ++r2)
            {
                mCells[k][r][r2].clear();
            }
        }
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            maybeAddMCell(Cell(i, j));
        }
    }
}

int TARGET;

void printMCells()
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (isMCell(Cell(i, j)))
                std::cerr << "x";
            else
                std::cerr << " ";
        }
        std::cerr << std::endl;
    }
}

void printRegions()
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (regions[i][j] < 10)
                std::cerr << regions[i][j];
            else
                std::cerr << (char)('a' + regions[i][j] - 10);
        }
        std::cerr << std::endl;
    }
}

int regCounts[MAX_D][MAX_P];

void findRegionCounts()
{
    for (int r = 0; r < d; ++r)
    {
        for (int k = 0; k < p; ++k)
        {
            regCounts[r][k] = 0;
        }
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            ++regCounts[regions[i][j]][votes[i][j]];
        }
    }
}

FORCE_INLINE
std::vector<int> findRanking(int r)
{
    std::vector<int> ranking(p);
    std::iota(ranking.begin(), ranking.end(), 0);
    std::sort(ranking.begin(), ranking.end(), [r](int x, int y)
              {
        if (regCounts[r][x] != regCounts[r][y]) return regCounts[r][x] > regCounts[r][y];
        else return x < y; });
    return ranking;
}

FORCE_INLINE
int findWinner(int r)
{
    int maxC = 0;
    int winner = p;
    for (int k = 0; k < p; ++k)
    {
        if (regCounts[r][k] > maxC)
        {
            maxC = regCounts[r][k];
            winner = k;
        }
        else if (regCounts[r][k] == maxC)
        {
            winner = p;
        }
    }
    return winner;
}

int findWon()
{
    int won = 0;
    for (int r = 0; r < d; ++r)
    {
        if (findWinner(r) == TARGET)
            ++won;
    }
    return won;
}

double findScore()
{
    return findWon() * partyScores[TARGET];
}

bool swapsDone;

FORCE_INLINE
void swapCellRegions(Cell a, Cell b)
{
    swapsDone = true;

    int &aReg = regions[a.x()][a.y()];
    int &bReg = regions[b.x()][b.y()];

    int aVote = votes[a.x()][a.y()];
    int bVote = votes[b.x()][b.y()];

    if (DEBUG_PRINT)
    {
        std::cerr << std::endl;
        std::cerr << std::endl;
        printRegions();
        std::cerr << std::endl;
        printMCells();
        std::cerr << std::endl;

        std::cerr << "      Swapping:"
                  << "\n        " << a.x() << " " << a.y() << " in " << aReg << " with " << aVote
                  << "\n        " << b.x() << " " << b.y() << " in " << bReg << " with " << bVote << std::endl;
    }

    assert(isMCell(a));
    assert(isMCell(b));
    assert(willTouchRegsAfterSwap(a, b));
    assert(aReg != bReg);

    maybeRemMCells3x3(a);
    maybeRemMCells3x3(b);

    ++regCounts[aReg][bVote];
    --regCounts[aReg][aVote];
    ++regCounts[bReg][aVote];
    --regCounts[bReg][bVote];

    std::swap(aReg, bReg);

    maybeAddMCells3x3(a);
    maybeAddMCells3x3(b);
}

bool isAdjRegion[MAX_D][MAX_D];
std::vector<int> adjRegions[MAX_D];

void optRegion(int r)
{
    std::vector<std::pair<Cell, Cell>> swaps;

    std::vector<int> others = adjRegions[r];

    if (DEBUG_PRINT)
    {
        std::cerr << "    region: " << r << std::endl;
    }

    bool isOk[MAX_P];

    while (findWinner(r) != TARGET)
    {
        std::vector<int> ranking = findRanking(r);

        Cell give(-1);
        Cell take(-1);
        for (int other : others)
        {
            std::fill(isOk, isOk + p, true);

            int otherWinner = findWinner(other);
            if (otherWinner == TARGET)
            {
                for (int k = 0; k < p; ++k)
                {
                    isOk[k] = k == otherWinner || regCounts[other][k] + 1 < regCounts[other][otherWinner];
                }
            }

            give.idx = -1;
            for (int k : ranking)
            {
                if (k == TARGET || !isOk[k])
                    continue;
                if (!mCells[k][r][other].empty())
                {
                    give = *mCells[k][r][other].begin();
                    break;
                }
            }

            if (give.idx == -1)
                continue;

            int giveVote = votes[give.x()][give.y()];

            if (otherWinner == TARGET)
            {
                isOk[otherWinner] = true;
                for (int k = 0; k < p; ++k)
                {
                    if (k == otherWinner)
                        continue;
                    isOk[k] = true;
                    if (regCounts[other][otherWinner] - 1 <= regCounts[other][k] + (giveVote == k))
                        isOk[otherWinner] = false;
                }
            }

            take.idx = -1;
            if (isOk[TARGET])
            {
                for (Cell maybeTake : mCells[TARGET][other][r])
                {
                    if (willTouchRegsAfterSwap(maybeTake, give))
                    {
                        take = maybeTake;
                        break;
                    }
                }
            }

            if (take.idx == -1 && regCounts[r][giveVote] >= regCounts[r][TARGET])
            {
                for (int k = 0; k < p; ++k)
                {
                    if (regCounts[r][k] + 1 >= regCounts[r][TARGET] || !isOk[k])
                        continue;
                    for (Cell maybeTake : mCells[k][other][r])
                    {
                        if (willTouchRegsAfterSwap(maybeTake, give))
                        {
                            take = maybeTake;
                            break;
                        }
                    }
                    if (take.idx != -1)
                        break;
                }
            }

            if (take.idx != -1)
                break;
        }

        if (give.idx == -1 || take.idx == -1)
            break;

        swaps.emplace_back(give, take);
        swapCellRegions(give, take);
    }

    if (DEBUG_PRINT)
    {
        std::cerr << "     region: " << r << " winner: " << findWinner(r) << std::endl;
    }

    if (findWinner(r) != TARGET)
    {
        while (!swaps.empty())
        {
            swapCellRegions(swaps.back().first, swaps.back().second);
            swaps.pop_back();
        }
    }
}

void findAdjRegions()
{
    for (int r = 0; r < d; ++r)
    {
        adjRegions[r].clear();
        for (int r2 = 0; r2 < d; ++r2)
        {
            isAdjRegion[r][r2] = false;
        }
    }
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (i > 0 && regions[i - 1][j] != regions[i][j])
            {
                isAdjRegion[regions[i - 1][j]][regions[i][j]] = true;
                isAdjRegion[regions[i][j]][regions[i - 1][j]] = true;
            }
            if (j > 0 && regions[i][j - 1] != regions[i][j])
            {
                isAdjRegion[regions[i][j - 1]][regions[i][j]] = true;
                isAdjRegion[regions[i][j]][regions[i][j - 1]] = true;
            }
        }
    }
    for (int r = 0; r < d; ++r)
    {
        for (int r2 = 0; r2 < d; ++r2)
        {
            if (isAdjRegion[r][r2])
            {
                adjRegions[r].push_back(r2);
                adjRegions[r2].push_back(r);
            }
        }
    }
}

void localOpt()
{
    findAdjRegions();

    for (int r = 0; r < d && timePassed() < CURR_MAX_TIME; ++r)
    {
        optRegion(r);
    }

    if (DEBUG_PRINT)
    {
        std::cerr << "   iter won: " << findWon() << std::endl;
    }
}

void solve()
{
    setup();

    makeSnake();

    double setupTime = timePassed();

    double totalScore = 0;

    for (TARGET = 0; TARGET < p && timePassed() < SOFT_MAX_TIME; ++TARGET)
    {
        std::cerr << " target: " << TARGET << " time: " << timePassed() << std::endl;

        CURR_MAX_TIME = (TARGET + 1) * (SOFT_MAX_TIME - setupTime) / p + setupTime;

        if (partyScores[TARGET] >= 0.0 / p) makePortRegions();
        else makeSnakeRegions();
        findMCells();
        findRegionCounts();

        std::cerr << "  Initial won: " << findWon() << std::endl;

        swapsDone = true;

        while (timePassed() - setupTime < CURR_MAX_TIME && swapsDone)
        {
            swapsDone = false;
            localOpt();
        }

        std::cerr << "  Final won: " << findWon() << std::endl;

        output();

        totalScore += findScore();
    }

    std::cerr << "Total score: " << totalScore << std::endl;

    while (TARGET < p)
    {
        output();
        ++TARGET;
    }
}

int main()
{
    generator.seed(0);
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    startT = std::chrono::high_resolution_clock::now();

    input();
    solve();

    std::cerr << "Time passed: " << timePassed() << std::endl;

    return 0;
}
