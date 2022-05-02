#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <numeric>
#include <limits>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <iterator>

const bool DEBUG = false;

std::ifstream inF("cleanup.in");
std::ofstream outF("cleanup.out");

std::mt19937 generator;
int randNum(int lb, int ub)
{
    std::uniform_int_distribution<int> distribution(lb, ub - 1);
    return distribution(generator);
}

double maxTimeLimit = 4.5;
std::chrono::high_resolution_clock::time_point startT, currT;

double timePassed()
{
    using namespace std::chrono;
    currT = high_resolution_clock::now();
    return duration_cast<duration<double>>(currT - startT).count();
}

typedef unsigned long long ull;

constexpr int MAX_N = 1e5;
constexpr int MAX_S = 1e5;
constexpr int MAX_K = 50;
constexpr int MAX_D = 141419; // floor((MAX_N - 1) * sqrt(2))
constexpr int MAX_EPOCH = 317; // floor(MAX_S / (sqrt(MAX_S) - 1))
constexpr int MAX_P = 1e5;

constexpr int NUM_DIRS = 4;
constexpr int UP = 0;
constexpr int DOWN = 1;
constexpr int RIGHT = 2;
constexpr int LEFT = 3;

constexpr char dirChars[NUM_DIRS] = {'U', 'D', 'R', 'L'};

constexpr int TAKEN_PENALTY = -MAX_P;

constexpr int MAX_DEPTH = 15;
constexpr int CHUNK_SIZE = 15;

struct Solution
{
    long long score = -1;
    int dirs[MAX_S];
};

int n, k, s, xStart, yStart;
ull W, Z;
int maxD;

struct BlackHole
{
    int x, y;

    BlackHole next() const
    {
        return BlackHole{(int) ((x * W + Z) % n + 1), (int) ((y * Z + W) % n + 1)};
    }

    int dist(int _x, int _y) const
    {
        return sqrt((long long) (x - _x) * (x - _x) + (long long) (y - _y) * (y - _y));
    }
};

int p[MAX_D + 1];
BlackHole holes[MAX_EPOCH + 1][MAX_K];
int epochs[MAX_S + 1];

int pows3[MAX_DEPTH + 1];

Solution sol, bestSol;

void input()
{
    inF >> n >> k >> s >> W >> Z >> xStart >> yStart;

    maxD = (n - 1) * sqrt(2);
    for (int i = 0; i <= maxD; ++i)
    {
        inF >> p[i];
    }

    for (int i = 0; i < k; ++i)
    {
        inF >> holes[0][i].x >> holes[0][i].y;
    }
}

void output()
{
    for (int i = 0; i < s; ++i)
    {
        outF << dirChars[sol.dirs[i]];
    }
    outF << std::endl;
}

void init()
{
    int epochLen = sqrt(s);
    for (int t = 0; t <= s; ++t)
    {
        epochs[t] = t / epochLen;
    }

    for (int i = 1; i <= epochs[s]; ++i)
    {
        for (int j = 0; j < k; ++j)
        {
            holes[i][j] = holes[i - 1][j].next();
        }
    }

    pows3[0] = 1;
    for (int i = 1; i <= MAX_DEPTH; ++i)
    {
        pows3[i] = pows3[i - 1] * 3;
    }
}

struct Location
{
    int x, y, epoch;

    bool operator==(const Location& other) const
    {
        return x == other.x && y == other.y && epoch == other.epoch;
    }

    bool operator<(const Location& other) const
    {
        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
        if (epoch != other.epoch) return epoch < other.epoch;
        return false;
    }
};

namespace std
{
template <>
struct hash<Location>
{
    size_t operator()(const Location& loc) const
    {
        return loc.x ^ (loc.y << 15) ^ (loc.epoch << 22) ^ loc.epoch;
    }
};
}

struct State
{
    bool taken = false;
    int score = -1;

    int currScore() const
    {
        if (taken) return 0;
        else return score;
    }
};

struct Chunk
{
    int baseX, baseY, epoch;
    int cenX, cenY;
    int cenClosest;
    double cenSecMinD;

    State states[CHUNK_SIZE][CHUNK_SIZE];

    Chunk(int baseX, int baseY, int epoch):
        baseX(baseX),
        baseY(baseY),
        epoch(epoch)
    {
        cenX = baseX + CHUNK_SIZE / 2;
        cenY = baseY + CHUNK_SIZE / 2;

        double cenMinD = maxD;
        cenSecMinD = maxD;
        for (int i = 0; i < k; ++i)
        {
            double currD = holes[epoch][i].dist(baseX, baseY);

            if (currD < cenMinD)
            {
                cenSecMinD = cenMinD;
                cenMinD = currD;
                cenClosest = i;
            }
            else if (currD < cenSecMinD)
            {
                cenSecMinD = currD;
            }
        }

        states[cenX - baseX][cenY - baseY].score = p[(int) (cenMinD)];
    }

    State& get(int x, int y)
    {
        if (x < baseX || x >= baseX + CHUNK_SIZE || y < baseY || y >= baseY + CHUNK_SIZE) exit(4);
        State& state = states[x - baseX][y - baseY];
        if (state.score >= 0) return state;

        int d = holes[epoch][cenClosest].dist(x, y);
        int secLb = cenSecMinD - sqrt((long long) (baseX - x) * (baseX - x) + (long long) (baseY - y) * (baseY - y));
        
        if (d > secLb)
        {
            for (int i = 0; i < k; ++i)
            {
                d = std::min(d, holes[epoch][i].dist(x, y));
            }
        }

        state.score = p[d];
        return state;
    }
};

std::unordered_map<Location, Chunk> chunkMap;

Chunk& getChunk(int x, int y, int epoch)
{
    if (x <= 0 || x > n || y <= 0 || y > n) exit(2);
    Location key = {(x - 1) / CHUNK_SIZE, (y - 1) / CHUNK_SIZE, epoch};
    auto chunkIt = chunkMap.find(key);
    if (chunkIt == chunkMap.end())
    {
        chunkIt = chunkMap.insert({key, Chunk(key.x * CHUNK_SIZE + 1, key.y * CHUNK_SIZE + 1, key.epoch)}).first;
    }
    return chunkIt->second;
}

State& getState(int x, int y, int epoch)
{
    return getChunk(x, y, epoch).get(x, y);
}

int chunkCacheEpoch;
int chunkCacheX;
int chunkCacheY;

Chunk* chunkCache[3][3];

Chunk& getChunkFast(int x, int y, int epoch)
{
    if (x <= 0 || x > n || y <= 0 || y > n) exit(2);
    if (epoch != chunkCacheEpoch) return getChunk(x, y, epoch);
    int offX = (x - 1) / CHUNK_SIZE - chunkCacheX;
    int offY = (y - 1) / CHUNK_SIZE - chunkCacheY;
    if (offX < 0 || offX > 2 || offY < 0 || offY > 2) return getChunk(x, y, epoch);
    if (chunkCache[offX][offY] == nullptr) exit(5);
    return *chunkCache[offX][offY];
}

State& getStateFast(int x, int y, int epoch)
{
    return getChunkFast(x, y, epoch).get(x, y);
}

void setChunkCache(int x, int y, int epoch)
{
    int newChunkX = (x - 1) / CHUNK_SIZE - 1;
    int newChunkY = (y - 1) / CHUNK_SIZE - 1;

    if (newChunkX == chunkCacheX && newChunkY == chunkCacheY && epoch == chunkCacheEpoch) return;

    chunkCacheX = newChunkX;
    chunkCacheY = newChunkY;
    chunkCacheEpoch = epoch;

    for (int offX = 0; offX < 3; ++offX)
    {
        for (int offY = 0; offY < 3; ++offY)
        {
            int x = (chunkCacheX + offX) * CHUNK_SIZE + 1;
            int y = (chunkCacheY + offY) * CHUNK_SIZE + 1;

            if (x <= 0 || x > n || y <= 0 || y > n) chunkCache[offX][offY] = nullptr;
            else chunkCache[offX][offY] = &getChunk(x, y, epoch);
        }
    }
}

void updateSol()
{
    if (sol.score <= bestSol.score) return;

    bestSol = sol;
    std::cerr << "New score: " << bestSol.score << std::endl;
}

std::pair<int, int> findScore(int x, int y, int t, int depth)
{
    if (t > s || depth == -1) return {0, -1};

    State& state = getStateFast(x, y, epochs[t]);
    if (state.taken) return {TAKEN_PENALTY, -1};
    if (depth == 0) return {state.score, -1};

    state.taken = true;

    long long score = TAKEN_PENALTY - 1;
    int dir = -1;
    for (int i = 0; i < NUM_DIRS; ++i)
    {
        int curr = TAKEN_PENALTY - 1;
        switch (i)
        {
        case UP:
            if (x == 1) break;
            curr = findScore(x - 1, y, t + 1, depth - 1).first;
            break;
        case DOWN:
            if (x == n) break;
            curr = findScore(x + 1, y, t + 1, depth - 1).first;
            break;
        case RIGHT:
            if (y == n) break;
            curr = findScore(x , y + 1, t + 1, depth - 1).first;
            break;
        case LEFT:
            if (y == 1) break;
            curr = findScore(x, y - 1, t + 1, depth - 1).first;
            break;
        }

        if (curr > score)
        {
            score = curr;
            dir = i;
        }
    }
    score += state.score;
    
    state.taken = false;
    return {score, dir};
}

int nextDir;

int findBestDir(int x, int y, int t, int depth)
{
    long long score = TAKEN_PENALTY - 1;
    int dir = -1;
    nextDir = -1;
    for (int i = 0; i < NUM_DIRS; ++i)
    {
        std::pair<int, int> curr = {TAKEN_PENALTY - 1, -1};
        switch (i)
        {
        case UP:
            if (x == 1) break;
            curr = findScore(x - 1, y, t + 1, depth - 1);
            break;
        case DOWN:
            if (x == n) break;
            curr = findScore(x + 1, y, t + 1, depth - 1);
            break;
        case RIGHT:
            if (y == n) break;
            curr = findScore(x , y + 1, t + 1, depth - 1);
            break;
        case LEFT:
            if (y == 1) break;
            curr = findScore(x, y - 1, t + 1, depth - 1);
            break;
        }

        if (curr.first > score)
        {
            score = curr.first;
            dir = i;
            nextDir = curr.second;
        }
    }

    return dir;
}

void greedy()
{
    int depth = 17 - log(s) / log(3);
    depth = std::max(depth, 1);

    chunkMap.clear();

    int x = xStart;
    int y = yStart;

    sol.score = 0;

    int takenCnt = 0;
    ull totalNodes = 0;

    nextDir = -1;

    double initTime = timePassed();
    for (int t = 0; t < s; ++t)
    {
        double time = timePassed();

        if (time > 4.5) depth = 0;
        else if (time > 4.1) depth = 1;
        else if (time < 3.9 && t > 0)
        {
            double projectedTime = time + (time - initTime) * (s - t) * depth / totalNodes;

            if (projectedTime > 4.0)
            {
                if (depth > 1) --depth;
            }
            else
            {
                nextDir = -1;
                if (depth < MAX_DEPTH) ++depth;
            }

            if (DEBUG)
            {
                std::cerr << totalNodes << " " << projectedTime << " : " << depth << std::endl;
            }
        }

        // if (time >= 2)
        // {
        //     exit(depth);
        // }

        setChunkCache(x, y, epochs[t]);
        State& state = getStateFast(x, y, epochs[t]);
        if (state.taken) ++takenCnt;
        sol.score += state.currScore();
        state.taken = true;

        if (nextDir != -1)
        {
            sol.dirs[t] = nextDir;
            nextDir = -1;
        }
        else
        {
            sol.dirs[t] = findBestDir(x, y, t, depth);
            totalNodes += depth;
        }

        if (sol.dirs[t] == -1) exit(3);

        if (DEBUG)
        {
            std::cerr << epochs[t] << " " << x << " " << y << " : " << dirChars[sol.dirs[t]] << std::endl;
        }

        switch (sol.dirs[t])
        {
        case UP:
            --x;
            break;
        case DOWN:
            ++x;
            break;
        case RIGHT:
            ++y;
            break;
        case LEFT:
            --y;
            break;
        }

        if (x <= 0 || x > n || y <= 0 || y > n) exit(1);
    }

    State& state = getStateFast(xStart, yStart, epochs[s]);
    if (state.taken) ++takenCnt;
    sol.score += state.currScore();
    state.taken = true;

    // if (takenCnt > 0) exit((takenCnt * 100 + s / 2) / s);

    if (DEBUG)
    {
        std::cerr << epochs[s] << " " << x << " " << y << std::endl;
    }

    updateSol();
}

void solve()
{
    // while (true)
    // {
    //     bestSol.score = -1;
    //     startT = std::chrono::high_resolution_clock::now();

    //     init();
    //     greedy();

    //     std::cerr << "Total time: " << timePassed() << std::endl;
    // }

    init();
    greedy();

    std::cerr << "Total time: " << timePassed() << std::endl;
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
