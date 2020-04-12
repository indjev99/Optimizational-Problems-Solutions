#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <queue>
#include <stdlib.h>
#include <math.h>
#include <unordered_map>
using namespace std;

const int test = 5;

const string tests[] = {
    "a_solar",
    "b_dream",
    "c_soup",
    "d_maelstrom",
    "e_igloos",
    "f_glitch"
};

ifstream inFile((tests[test]+".txt").c_str());
ofstream outFile(("sol_"+tests[test]+".txt").c_str());

struct Person
{
    int c;
    int b;
    long long s1 = 0, s2 = 0;
    int x = -1, y = -1;
};

const int MAX_D = 1e5;
const int MAX_M = 2e4;
const int MAX_WH = 1e3;

int w, h;

struct Point
{
    int x, y;
};

vector<Point> devSpots;
vector<Point> manSpots;

#define WAL 0
#define DEV 1
#define MAN 2

struct Cell
{
    int tile;
    int num = -1;
};

Cell office[MAX_WH][MAX_WH];

int numDevs, numMans;
Person devs[MAX_D];
Person mans[MAX_M];

unordered_map<string, int> comps, skills;

int currComp = 0;
int enterCompany(const string& comp)
{
    if (comps.find(comp) == comps.end())
    {
        comps[comp] = currComp++;
    }
    return comps[comp];
}

int currSkill = 0;
int enterSkill(const string& skill)
{
    if (skills.find(skill) == skills.end())
    {
        skills[skill] = currSkill++;
    }
    return skills[skill];
}

void input()
{
    string row, comp, skill;
    inFile >> w >> h;
    for (int i = 0; i < h; ++i)
    {
        inFile >> row;
        for (int j = 0; j < w; ++ j)
        {
            if (row[j] == '#') office[j][i].tile = WAL;
            if (row[j] == '_')
            {
                office[j][i].tile = DEV;
                devSpots.push_back({j, i});
            }
            if (row[j] == 'M')
            {
                office[j][i].tile = MAN;
                manSpots.push_back({j, i});
            }
        }
    }
    inFile >> numDevs;
    for (int i = 0; i < numDevs; ++i)
    {
        int numSkills;
        inFile >> comp >> devs[i].b >> numSkills;
        devs[i].c = enterCompany(comp);
        for (int j = 0; j < numSkills; ++j)
        {
            inFile >> skill;
            int currSkill = enterSkill(skill);
            if (currSkill < 64) devs[i].s1 = devs[i].s1 | (1ll << currSkill);
            else devs[i].s2 = devs[i].s2 | (1ll << (currSkill - 64));
        }
    }
    inFile >> numMans;
    for (int i = 0; i < numMans; ++i)
    {
        inFile >> comp >> mans[i].b;
        mans[i].c = enterCompany(comp);
    }
}

int currAns = -1;
int ans = -1;
string outputStr;

void output()
{
    outFile << outputStr;
}

void generateOutput()
{
    ostringstream ss;
    for (int i = 0; i < numDevs; ++i)
    {
        if (devs[i].x >= 0) ss << devs[i].x << ' ' << devs[i].y;
        else ss << 'X';
        ss << '\n';
    }
    for (int i = 0; i < numMans; ++i)
    {
        if (mans[i].x >= 0) ss << mans[i].x << ' ' << mans[i].y;
        else ss << 'X';
        ss << '\n';
    }
    outputStr = ss.str();
}

inline int countSetBits(long long s)
{
    int cnt = 0;
    while (s)
    {
        s &= (s - 1);
        ++cnt;
    }
    return cnt;
}

inline int tp(const Person& p1, const Person& p2)
{
    int bp = (p1.c == p2.c) * p1.b * p2.b;
    int common = countSetBits(p1.s1 & p2.s1) + countSetBits(p1.s2 & p2.s2);
    int total = countSetBits(p1.s1 | p2.s1) + countSetBits(p1.s2 | p2.s2);
    int wp = common * (total - common);
    return bp + wp;
}

inline Person* getPerson(int x, int y, bool printNum=false)
{
    if (x >= w || y >= h || x < 0 || y < 0) return nullptr;
    const Cell& cell = office[x][y];
    if (cell.tile == WAL || cell.num < 0) return nullptr;
    if (cell.tile == DEV) return &devs[cell.num];
    else return &mans[cell.num];
}

int compVal[MAX_D];

bool customCmpDevs(int a, int b)
{
    const Person& pa = devs[a];
    const Person& pb = devs[b];
    if (pa.c != pb.c) return compVal[pa.c] < compVal[pb.c];
    else return pa.b > pb.b;
}

bool customCmpMans(int a, int b)
{
    const Person& pa = mans[a];
    const Person& pb = mans[b];
    if (pa.c != pb.c) return compVal[pa.c] < compVal[pb.c];
    else return pa.b > pb.b;
}

int currDI, currMI;
vector<int> devInds;
vector<int> manInds;

int currT = 0;
int visT[MAX_WH][MAX_WH];

void reset()
{
    currDI = 0;
    currMI = 0;
    ++currT;
    for (int i = 0; i < numDevs; ++i)
    {
        devs[i].x = -1;
        devs[i].y = -1;
    }
    for (int i = 0; i < numMans; ++i)
    {
        mans[i].x = -1;
        mans[i].y = -1;
    }
    for (int i = 0; i < devSpots.size(); ++i)
    {
        office[devSpots[i].x][devSpots[i].y].num = -1;
    }
    for (int i = 0; i < manSpots.size(); ++i)
    {
        office[manSpots[i].x][manSpots[i].y].num = -1;
    }
}

void logScore(int score)
{
    currAns = score;
    if (score > ans)
    {
        ans = score;
        cerr << "Score: " << score << endl;
        generateOutput();
        //cerr << outputStr << endl;
    }
}

void eval(bool print=false)
{
    int score = 0;
    for (int i = 0; i < w; ++i)
    {
        for (int j = 0; j < h; ++j)
        {
            Person* here = getPerson(i,j);
            if (here == nullptr) continue;
            Person* right = getPerson(i + 1, j);
            Person* down = getPerson(i, j + 1);
            if (right != nullptr) score += tp(*here, *right);
            if (down != nullptr) score += tp(*here, *down);
        }
    }
    if (print)
    {
        cerr << "Eval Score: " << score << endl;
    }
    logScore(score);
}

void solveDummy()
{
    reset();
    for (int i = 0; i < numDevs && i < devSpots.size(); ++i)
    {
        devs[i].x = devSpots[i].x;
        devs[i].y = devSpots[i].y;
        office[devSpots[i].x][devSpots[i].y].num = i;
    }
    for (int i = 0; i < numMans && i < manSpots.size(); ++i)
    {
        mans[i].x = manSpots[i].x;
        mans[i].y = manSpots[i].y;
        office[manSpots[i].x][manSpots[i].y].num = i;
    }
    eval();
}

void DFS(int x, int y)
{
    if (x < 0 || x >= w || y < 0 || y >= h || office[x][y].tile == WAL || visT[x][y] == currT) return;
    visT[x][y] = currT;
    if (office[x][y].tile == DEV && currDI < devInds.size())
    {
        devs[devInds[currDI]].x = x;
        devs[devInds[currDI]].y = y;
        office[x][y].num = devInds[currDI++];
    }
    else if (office[x][y].tile == MAN && currMI < manInds.size())
    {
        mans[manInds[currMI]].x = x;
        mans[manInds[currMI]].y = y;
        office[x][y].num = manInds[currMI++];
    }
    DFS(x + 1, y);
    DFS(x, y + 1);
    DFS(x - 1, y);
    DFS(x, y - 1);
}

void shadyBFS(int x, int y)
{
    if (x < 0 || x >= w || y < 0 || y >= h || office[x][y].tile == WAL || visT[x][y] == currT) return;
    int currComp = currDI < devInds.size() ? devs[devInds[currDI]].c : -1;
    queue<Point> q;
    q.push({x, y});
    while (!q.empty())
    {
        Point p = q.front();
        q.pop();
        x = p.x;
        y = p.y;
        if (x < 0 || x >= w || y < 0 || y >= h || office[x][y].tile == WAL || visT[x][y] == currT) continue;
        visT[x][y] = currT;
        if (office[x][y].tile == DEV && currDI < devInds.size())
        {
            devs[devInds[currDI]].x = x;
            devs[devInds[currDI]].y = y;
            office[x][y].num = devInds[currDI++];
        }
        else if (office[x][y].tile == MAN && currMI < manInds.size() && (mans[manInds[currMI]].c == currComp || currComp == -1))
        {
            mans[manInds[currMI]].x = x;
            mans[manInds[currMI]].y = y;
            office[x][y].num = manInds[currMI++];
        }
        if (currComp != -1 && (currDI == devInds.size() || devs[devInds[currDI]].c != currComp)) return;
        q.push({x + 1, y});
        q.push({x, y + 1});
        q.push({x - 1, y});
        q.push({x, y - 1});
    }
}

void solveSort(bool bfs)
{
    reset();

    random_shuffle(compVal, compVal + currComp);
    random_shuffle(devSpots.begin(), devSpots.end());
    random_shuffle(manSpots.begin(), manSpots.end());

    sort(devInds.begin(), devInds.end(), customCmpDevs);
    sort(manInds.begin(), manInds.end(), customCmpMans);

    ++currT;
    for (const Point& spot : devSpots)
    {
        if (bfs) shadyBFS(spot.x, spot.y);
        else DFS(spot.x, spot.y);
    }
    for (const Point& spot : manSpots)
    {
        if (bfs) shadyBFS(spot.x, spot.y);
        else DFS(spot.x, spot.y);
    }

    for (int i = 0; i < devSpots.size() && currDI < devInds.size(); ++i)
    {
        if (office[devSpots[i].x][devSpots[i].y].num >= 0) continue;
        devs[devInds[currDI]].x = devSpots[i].x;
        devs[devInds[currDI]].y = devSpots[i].y;
        office[devSpots[i].x][devSpots[i].y].num = devInds[currDI++];
    }
    for (int i = 0; i < manSpots.size() && currMI < manInds.size(); ++i)
    {
        if (office[manSpots[i].x][manSpots[i].y].num >= 0) continue;
        devs[manInds[currMI]].x = manSpots[i].x;
        devs[manInds[currMI]].y = manSpots[i].y;
        office[manSpots[i].x][manSpots[i].y].num = manInds[currMI++];
    }

    eval();
}

void init()
{
    for (int i = 0; i < currComp; ++i)
    {
        compVal[i] = i;
    }
    for (int i = 0; i < numDevs; ++i)
    {
        devInds.push_back(i);
    }
    for (int i = 0; i < numMans; ++i)
    {
        manInds.push_back(i);
    }
}

int valueOf(const Person& p, int x, int y, const Person* other)
{
    if (x < 0) return 0;
    int value = 0;
    Person* right = getPerson(x + 1, y, true);
    Person* down = getPerson(x, y + 1, true);
    Person* left = getPerson(x - 1, y, true);
    Person* up = getPerson(x, y - 1, true);
    if (right != nullptr && right != other) value += tp(p, *right);
    if (down != nullptr && down != other) value += tp(p, *down);
    if (left != nullptr && left != other) value += tp(p, *left);
    if (up != nullptr && up != other) value += tp(p, *up);
    return value;
}

bool trySwapDevs(int i, int j)
{
    if (devs[i].x < 0 && devs[j].x < 0) return false;
    int currValue = valueOf(devs[i], devs[i].x, devs[i].y, &devs[j]) + valueOf(devs[j], devs[j].x, devs[j].y, &devs[i]);
    int newValue = valueOf(devs[i], devs[j].x, devs[j].y, &devs[i]) + valueOf(devs[j], devs[i].x, devs[i].y, &devs[j]);
    if (newValue > currValue)
    {
        swap(devs[i].x, devs[j].x);
        swap(devs[i].y, devs[j].y);
        if (devs[i].x >= 0) office[devs[i].x][devs[i].y].num = i;
        if (devs[j].x >= 0) office[devs[j].x][devs[j].y].num = j;
        currAns += newValue - currValue;
        return true;
    }
    return false;
}

bool trySwapMans(int i, int j)
{
    if (mans[i].x < 0 && mans[j].x < 0) return false;
    int currValue = valueOf(mans[i], mans[i].x, mans[i].y, &mans[j]) + valueOf(mans[j], mans[j].x, mans[j].y, &mans[i]);
    int newValue = valueOf(mans[i], mans[j].x, mans[j].y, &mans[i]) + valueOf(mans[j], mans[i].x, mans[i].y, &mans[j]);
    if (newValue > currValue)
    {
        swap(mans[i].x, mans[j].x);
        swap(mans[i].y, mans[j].y);
        if (mans[i].x >= 0) office[mans[i].x][mans[i].y].num = i;
        if (mans[j].x >= 0) office[mans[j].x][mans[j].y].num = j;
        currAns += newValue - currValue;
        return true;
    }
    return false;
}

const int RANGE = 300;
const int NUM_ITERS = 100000;

void localOpt()
{
    for (int curr = 0; curr < NUM_ITERS; ++curr)
    {
        if (rand() % 4)
        {
            int i = rand() % (numDevs - RANGE);
            int j = i + 1 + rand() % RANGE;
            trySwapDevs(devInds[i], devInds[j]);
        }
        else
        {
            int i = rand() % (numMans - RANGE);
            int j = i + 1 + rand() % RANGE;
            trySwapMans(manInds[i], manInds[j]);
        }
    }
    eval();
}

int main()
{
    srand(55);
    input();
    init();
    for (int i = 0; i < 100; ++i)
    {
        solveSort(false);
        localOpt();
    }
    output();
    return 0;
}

