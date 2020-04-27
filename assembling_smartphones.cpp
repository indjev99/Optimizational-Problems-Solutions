#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <map>
#include <stdlib.h>
#include <math.h>
using namespace std;

//#define DEBUG
//#define DEBUG2

const int test = 5;

const string tests[] = {
    "a_example",
    "b_single_arm",
    "c_few_arms",
    "d_tight_schedule",
    "e_dense_workspace",
    "f_decentralized",
    "g_mine"
};

ifstream inFile((tests[test]+".txt").c_str());
ofstream outFile(("sol_"+tests[test]+".txt").c_str());

struct Point
{
    int x, y;
};


bool operator<(const Point& self, const Point& other)
{
    if (self.x != other.x) return self.x < other.x;
    return self.y < other.y;
}

struct Task
{
    int num;
    int s, p;
    vector<Point> aps;
    double heur;
};

#define Right 0
#define Left 1
#define Up 2
#define Down 3
#define Wait 4

const string codes = "RLUDW";

const int revIns[5] = {1, 0, 3, 2, 4};
const int insDx[5] = {1, -1, 0, 0, 0};
const int insDy[5] = {0, 0, 1, -1, 0};

struct RoboArm
{
    Point mp;
    vector<int> tskNums;
    vector<int> ins;

    int ori = 0;
    Point pos;
    vector<int> currMovs;
};

const int MAX_WH = 1e3;
const int MAX_M = 1e3;
const int MAX_T = 1e3;
const int MAX_R = 1e2;

int w, h, r, m, t, l;
Point mnts[MAX_T];
Task tsks[MAX_T];

int ansR;
RoboArm ansRobs[MAX_R];

string outputStr;

void input()
{
    inFile >> w >> h >> r >> m >> t >> l;
    for (int i = 0; i < m; ++i)
    {
        inFile >> mnts[i].x >> mnts[i].y;
    }
    for (int i = 0; i < t; ++i)
    {
        tsks[i].num = i;
        inFile >> tsks[i].s >> tsks[i].p;
        for (int j = 0; j < tsks[i].p; ++j)
        {
            int x, y;
            inFile >> x >> y;
            tsks[i].aps.push_back({x, y});
        }
    }
}

void output()
{
    outFile << outputStr;
}

void generateOutput()
{
    ostringstream ss;
    ss << ansR << "\n";
    for (int i = 0; i < ansR; ++i)
    {
        const RoboArm& curr = ansRobs[i];
        ss << curr.mp.x << " " << curr.mp.y << " " << curr.tskNums.size() << " " << curr.ins.size() << "\n";
        for (int tn : curr.tskNums)
        {
            ss << tn << " ";
        }
        ss << "\n";
        for (int ins : curr.ins)
        {
            ss << codes[ins] << " ";
        }
        ss << "\n";
    }
    outputStr = ss.str();
}

bool cmpHeur(const Task& l, const Task& r)
{
    return l.heur > r.heur;
}

void sortByHeur()
{
    for (int i = 0; i < t; ++i)
    {
        Task& curr = tsks[i];
        long long totDist = 1;
        Point prev = curr.aps[0];
        for (const Point& ap : curr.aps)
        {
            totDist += abs(ap.x - prev.x) + abs(ap.y - prev.y);
            prev = ap;
        }
        if (totDist > l) curr.heur = -1;
        else curr.heur = curr.s / totDist;
    }
    sort(tsks, tsks + t, cmpHeur);
}

int timeFreedPos[MAX_WH][MAX_WH];

void setUsedPoints()
{
    for (int i = 0; i < w; ++i)
    {
        for (int j = 0; j < h; ++j)
        {
            timeFreedPos[i][j] = 0;
        }
    }
    for (int i = 0; i < m; ++i)
    {
        timeFreedPos[mnts[i].x][mnts[i].y] = 2 * l;
    }
}

bool usedMnts[MAX_M];

int findMnt(const Point& target)
{
    int best = -1, minDist;
    for (int i = 0; i < m; ++i)
    {
        if (usedMnts[i]) continue;
        int dist = abs(target.x - mnts[i].x) + abs(target.y - mnts[i].y);
        if (i == 0 || dist < minDist)
        {
            minDist = dist;
            best = i;
        }
    }
    return best;
}

map<Point, int> timeFreedCurr;
map<Point, int> timeFreedBig;

void moveRobo(RoboArm& robo, int ins, bool recordIns=true)
{
    if (ins == Wait) return;
    Point oldPos = robo.pos;
    robo.pos.x += insDx[ins];
    robo.pos.y += insDy[ins];
    if (robo.currMovs.empty() ||
        revIns[robo.currMovs[robo.currMovs.size() - 1]] != ins)
    {
        if (recordIns) timeFreedCurr[robo.pos] = 2 * l;
        robo.currMovs.push_back(ins);
    }
    else
    {
        if (recordIns) timeFreedCurr[oldPos] = robo.ins.size();
        robo.currMovs.pop_back();
    }
    if (recordIns) robo.ins.push_back(ins);
    #ifdef DEBUG2
    if (recordIns) cerr << codes[ins] << endl;
    #endif // DEBUG2
}

void traceBackOnce(RoboArm& robo)
{
    int lastIns = robo.ins.back();
    robo.ins.pop_back();
    moveRobo(robo, revIns[lastIns], false);
    #ifdef DEBUG2
    cerr << "-" << codes[lastIns] << endl;
    #endif // DEBUG2
}

void traceBack(RoboArm& robo, int num)
{
    while (num--) traceBackOnce(robo);
}

bool isUseful(int dx, int dy, int ins)
{
    int currDx = insDx[ins];
    int currDy = insDy[ins];
    return currDx > 0 && dx > 0 || currDx < 0 && dx < 0 ||
           currDy > 0 && dy > 0 || currDy < 0 && dy < 0;
}

bool goToPoint(const Point& target, RoboArm& robo)
{
    int dx = target.x - robo.pos.x;
    int dy = target.y - robo.pos.y;
    while (!robo.currMovs.empty())
    {
        int backIns = revIns[robo.currMovs.back()];
        if (isUseful(dx, dy, backIns))
        {
            #ifdef DEBUG
            cerr << "    robot at " << robo.pos.x << " " << robo.pos.y << endl;
            #endif // DEBUG
            moveRobo(robo, backIns);
            dx = target.x - robo.pos.x;
            dy = target.y - robo.pos.y;
        }
        else break;
    }
    if (robo.currMovs.empty()) robo.ori = 0;
    else if (robo.currMovs.back() % 2 == 0) robo.ori = 1;
    else robo.ori = -1;
    bool wentRevOri = false;
    int x = robo.pos.x, y = robo.pos.y;
    bool potential;
    while (dx != 0 || dy != 0 && robo.ins.size() < l)
    {
        potential = false;
        #ifdef DEBUG
        cerr << "    robot at " << robo.pos.x << " " << robo.pos.y << endl;
        #endif // DEBUG
        if ((dx > 0 && x + 1 < w && robo.ins.size() >= timeFreedPos[x + 1][y]) &&
            (robo.ori >= 0 || dy == 0))
        {
            potential = true;
            moveRobo(robo, Right);
            ++x; --dx;
            if (robo.ori == 0) robo.ori = 1;
            else if (robo.ori < 0) wentRevOri = true;
            continue;
        }
        if ((dx < 0 && x - 1 >= 0 && robo.ins.size() >= timeFreedPos[x - 1][y]) &&
            (robo.ori <= 0 || dy == 0))
        {
            potential = true;
            moveRobo(robo, Left);
            --x; ++dx;
            if (robo.ori == 0) robo.ori = -1;
            else if (robo.ori > 0) wentRevOri = true;
            continue;
        }
        if ((dy > 0 && y + 1 < h && robo.ins.size() >= timeFreedPos[x][y + 1]) &&
            (robo.ori >= 0 || dx == 0))
        {
            potential = true;
            moveRobo(robo, Up);
            ++y; --dy;
            if (robo.ori == 0) robo.ori = 1;
            else if (robo.ori < 0) wentRevOri = true;
            continue;
        }
        if ((dy < 0 && y - 1 >= 0 && robo.ins.size() >= timeFreedPos[x][y - 1]) &&
            (robo.ori <= 0 || dx == 0))
        {
            potential = true;
            moveRobo(robo, Down);
            --y; ++dy;
            if (robo.ori == 0) robo.ori = -1;
            else if (robo.ori > 0) wentRevOri = true;
            continue;
        }
        break;
    }
    #ifdef DEBUG
    cerr << "    robot at " << robo.pos.x << " " << robo.pos.y << endl;
    #endif // DEBUG
    if (wentRevOri && dx == 0 && dy == 0)
    {
        #ifdef DEBUG
        cerr << "   restoring orientation" << endl;
        #endif // DEBUG
        int targetMov = robo.currMovs.back();
        while (robo.currMovs.back() == targetMov && robo.ins.size() <= l)
        {
            #ifdef DEBUG
            cerr << "    robot at " << robo.pos.x << " " << robo.pos.y << endl;
            #endif // DEBUG
            moveRobo(robo, revIns[targetMov]);
        }
        #ifdef DEBUG
        cerr << "    robot at " << robo.pos.x << " " << robo.pos.y << endl;
        #endif // DEBUG
    }
    #ifdef DEBUG
    if (dx == 0 && dy == 0 && robo.ins.size() <= l)
        cerr << "   reached point " << target.x << " " << target.y << endl;
    #endif // DEBUG
    return dx == 0 && dy == 0 && robo.ins.size() <= l;
}

bool solveTask(const Task& curr, RoboArm& robo)
{
    #ifdef DEBUG
    cerr << "  attempting task " << curr.num << endl;
    #endif // DEBUG
    timeFreedCurr.clear();
    const int initZ = robo.ins.size();
    bool succ = true;
    for (const Point& ap : curr.aps)
    {
        #ifdef DEBUG
        cerr << "   robot heading for " << ap.x << " " << ap.y << endl;
        #endif // DEBUG
        if (!goToPoint(ap, robo))
        {
            succ = false;
            break;
        }
    }
    #ifdef DEBUG
    if (succ) cerr << "  task completed" << endl;
    #endif // DEBUG
    if (succ)
    {
        for (const auto& pr : timeFreedCurr)
        {
            timeFreedBig[pr.first] = pr.second;
        }
    }
    if (!succ) traceBack(robo, robo.ins.size() - initZ);
    else robo.tskNums.push_back(curr.num);
    return succ;
}

bool usedTasks[MAX_T];

bool solveRob(RoboArm& robo)
{
    robo.tskNums.clear();
    robo.currMovs.clear();

    timeFreedBig.clear();
    int firstTask = 0;
    while (firstTask < t && usedTasks[firstTask]) ++firstTask;
    if (firstTask >= t) return false;
    int mntInd = findMnt(tsks[firstTask].aps[0]);
    if (mntInd == -1) return false;
    bool nothing = true;

    for (int nextTask = firstTask; nextTask < t && robo.ins.size() < l; ++nextTask)
    {
        if (usedTasks[nextTask]) continue;
        if (nothing)
        {
            mntInd = findMnt(tsks[nextTask].aps[0]);
            robo.mp = mnts[mntInd];
            robo.pos = robo.mp;
        }
        if (solveTask(tsks[nextTask], robo))
        {
            nothing = false;
            usedTasks[nextTask] = true;
        }
    }

    for (int i = robo.ins.size(); i < l && !robo.currMovs.empty(); ++i)
    {
        moveRobo(robo, revIns[robo.currMovs.back()]);
    }

    if (nothing)
    {
        usedTasks[firstTask] = true;
        robo.tskNums.clear();
        robo.currMovs.clear();
        return true;
    }

    usedMnts[mntInd] = true;
    for (const auto& pr : timeFreedBig)
    {
        timeFreedPos[pr.first.x][pr.first.y] = pr.second;
    }

    return true;
}

void solve()
{
    setUsedPoints();
    sortByHeur();

    ansR = r;
    for (int i = 0; i < ansR; ++i)
    {
        #ifdef DEBUG
        cerr << " starting robot " << i << endl;
        #endif // DEBUG
        if (!solveRob(ansRobs[i])) ansR = i;
        else if (ansRobs[i].tskNums.empty()) --i;
    }

    generateOutput();
}

int main()
{
    srand(0);
    input();
    solve();
    output();
    return 0;
}

