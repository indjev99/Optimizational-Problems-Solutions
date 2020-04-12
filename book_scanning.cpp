#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
using namespace std;

const int test = 5;

const string tests[] = {
    "a_example",
    "b_read_on",
    "c_incunabula",
    "d_tough_choices",
    "e_so_many_books",
    "f_libraries_of_the_world"
};

ifstream inFile((tests[test]+".txt").c_str());
ofstream outFile(("sol_"+tests[test]+".txt").c_str());

struct Library
{
    int libId;
    int n, t, m;
    vector<int> ids;
    double libScore;
    int canTake;
    int toTake;
    vector<int> taken;
    int takenScore = 0;
};

void swap(Library& lib1, Library& lib2)
{
    swap(lib1.libId, lib2.libId);
    swap(lib1.n, lib2.n);
    swap(lib1.t, lib2.t);
    swap(lib1.m, lib2.m);
    swap(lib1.ids, lib2.ids);
    swap(lib1.libScore, lib2.libScore);
    swap(lib1.canTake, lib2.canTake);
    swap(lib1.toTake, lib2.toTake);
    swap(lib1.taken, lib2.taken);
    swap(lib1.takenScore, lib2.takenScore);
}

const int MAX_B = 1e5;
const int MAX_L = 1e5;

int b, l, d;
int score[MAX_B];
int freq[MAX_B];
Library lib[MAX_L];
int signedUp;
int lastSignedUp;
int totalScore;
string outputStr;

void input()
{
    inFile >> b >> l >> d;
    for (int i = 0; i < b; ++i)
    {
        inFile >> score[i];
    }
    for (int i = 0; i < l; ++i)
    {
        lib[i].libId = i;
        inFile >> lib[i].n >> lib[i].t >> lib[i].m;
        for (int j = 0; j < lib[i].n; ++j)
        {
            int id;
            inFile >> id;
            ++freq[id];
            lib[i].ids.push_back(id);
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
    ss << signedUp << "\n";
    for (int i = 0; i < lastSignedUp; ++i)
    {
        if (lib[i].taken.empty()) continue;
        ss << lib[i].libId << " " << lib[i].taken.size() << "\n";
        for (int j = 0; j < lib[i].taken.size(); ++j)
        {
            if (j > 0) ss << " ";
            ss << lib[i].taken[j];
        }
        ss << "\n";
    }
    outputStr = ss.str();
}

bool isFree[MAX_B];

int getScore(int id)
{
    return isFree[id] * score[id];
}

bool cmpScore(int id1, int id2)
{
    return getScore(id1) > getScore(id2);
}

bool cmpLibScore(const Library& lib1, const Library& lib2)
{
    return lib1.libScore > lib2.libScore;
}

int bestScore = -1;

void solveIter(bool toSort)
{
    for (int i = 0; i < b; ++i)
    {
        isFree[i] = true;
    }
    if (toSort) sort(lib, lib + l, cmpLibScore);
    int currDay = 0;
    lastSignedUp = 0;
    signedUp = 0;
    totalScore = 0;
    while (lastSignedUp < l && currDay < d)
    {
        Library& curr = lib[lastSignedUp++];
        currDay += curr.t;
        if (currDay > d)
        {
            --lastSignedUp;
            break;
        }
        sort(curr.ids.begin(), curr.ids.end(), cmpScore);
        curr.taken.resize(0);
        curr.canTake = (d - currDay) * curr.m;
        if (curr.canTake > curr.ids.size()) curr.canTake = curr.ids.size();
        for (int i = 0; i < curr.canTake; ++i)
        {
            int currScore = getScore(curr.ids[i]);
            isFree[curr.ids[i]] = false;
            if (currScore > 0)
            {
                totalScore += currScore;
                curr.takenScore += currScore;
                curr.taken.push_back(curr.ids[i]);
            }
        }
        if (curr.taken.empty())
        {
            currDay -= curr.t;
        }
        else signedUp++;
    }
}

void solveIterSmart(bool toSort)
{
    for (int i = 0; i < b; ++i)
    {
        isFree[i] = true;
    }
    if (toSort) sort(lib, lib + l, cmpLibScore);
    int currDay = 0;
    lastSignedUp = 0;
    signedUp = 0;
    totalScore = 0;
    while (lastSignedUp < l && currDay < d)
    {
        Library& curr = lib[lastSignedUp++];
        currDay += curr.t;
        if (currDay > d)
        {
            --lastSignedUp;
            break;
        }
        ++signedUp;
        sort(curr.ids.begin(), curr.ids.end(), cmpScore);
        curr.taken.resize(0);
        curr.toTake = 0;
        curr.canTake = (d - currDay) * curr.m;
        if (curr.canTake > curr.ids.size()) curr.canTake = curr.ids.size();
    }
    bool took = true;
    while (took)
    {
        took = false;
        for (int i = lastSignedUp - 1; i >= 0; --i)
        {
            Library& curr = lib[i];
            if (!curr.canTake) continue;
            while (curr.toTake < curr.n && !isFree[curr.ids[curr.toTake]]) ++curr.toTake;
            if (curr.toTake < curr.n)
            {
                took = true;
                isFree[curr.ids[curr.toTake]] = false;
                totalScore += score[curr.ids[curr.toTake]];
                curr.takenScore += score[curr.ids[curr.toTake]];
                curr.taken.push_back(curr.ids[curr.toTake]);
                --curr.canTake;
            }
            else curr.canTake = 0;
        }
    }
    for (int i = 0; i < lastSignedUp; ++i)
    {
        if (lib[i].taken.empty()) --signedUp;
    }
}

double customScore(int id)
{
    return score[id] * 1.0 / freq[id];
}

void evalLib(Library& curr, int days)
{
    curr.libScore = 0;
    for (int j = 0; j < curr.ids.size(); ++j)
    {
        curr.libScore += customScore(curr.ids[j]);
    }

    if (test == 5)
    {
        curr.libScore /= pow(curr.t, 1);
        curr.libScore *= (rand() % 101) / 100.0 * 0.4 + 0.85;

    }
    else if (test == 4)
    {
        curr.libScore /= pow(curr.t, 1.75);
        curr.libScore *= pow(curr.m, 1);
        curr.libScore *= (rand() % 101) / 100.0 * 0.3 + 0.85;
    }
    else
    {
        curr.libScore /= curr.t;
    }
}

void customSort()
{
    for (int i = 0; i < b; ++i)
    {
        freq[i] = 1;
    }
    int currDay = 0;
    for (int i = 0; i < b && currDay < d + 10; ++i)
    {
        double maxS = -1;
        int maxJ = -1;
        for (int j = i; j < b; ++j)
        {
            evalLib(lib[j], d - currDay);
            if (lib[j].libScore > maxS)
            {
                maxS = lib[j].libScore;
                maxJ = j;
            }
        }
        swap(lib[i], lib[maxJ]);
        for (int id : lib[i].ids)
        {
            ++freq[id];
        }
        currDay += lib[i].t;
    }
}

const int REPS = 10;

void solve()
{
    for (int i = 0; i < b; ++i)
    {
        isFree[i] = true;
    }
    int reps = REPS;
    customSort();
    bool first = true;
    while (reps-- > 0)
    {
        solveIterSmart(!first);
        if (totalScore > bestScore)
        {
            if (!first) cerr << "not first" << endl;
            cerr << totalScore << " " << signedUp << endl;
            bestScore = totalScore;
            generateOutput();
            reps = REPS;
        }
        for (int i = 0; i < signedUp; ++i)
        {
            lib[i].libScore = lib[i].takenScore;
            lib[i].takenScore = 0;
        }
        first = false;
    }
}

void bigSolve()
{
    for (int i = 0; i < 50; ++i)
    {
        solve();
    }
}

int main()
{
    srand(-5);
    input();
    bigSolve();
    output();
    return 0;
}

