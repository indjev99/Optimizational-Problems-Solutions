#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <map>
#include <stdlib.h>
#include <math.h>
#include <unordered_map>
using namespace std;

//#define DEBUG

int test;

const string tests[] = {
    "a",
    "b",
    "c",
    "d",
    "e",
    "f"
};

ifstream inFile;
ofstream outFile;

struct Inter
{
    int num;
    std::vector<int> in;
    std::vector<int> out;

    int cars;
};

struct Street
{
    int num;
    int b, e;
    int len;

    int cars;
    int time;
};

struct Car
{
    int num;
    std::vector<int> path;
    int totLen;
};

const int MAX_I = 1e5;
const int MAX_S = 1e5;
const int MAX_V = 1e3;

int d, numI, numS, numV, f;
Inter intr[MAX_I];
Street str[MAX_S];
Car car[MAX_V];

string outputStr;

unordered_map<string, int> nameToNumMap;
vector<string> numToName;
int nameToNum(const string& name)
{
    if (nameToNumMap.find(name) == nameToNumMap.end())
    {
        nameToNumMap[name] = numToName.size();
        numToName.push_back(name);
    }
    return nameToNumMap[name];
}

void input()
{
    string name;
    int p;
    inFile >> d >> numI >> numS >> numV >> f;
    for (int i = 0; i < numS; ++i)
    {
        intr[i].num = i;
    }
    for (int i = 0; i < numS; ++i)
    {
        inFile >> str[i].b >> str[i].e >> name >> str[i].len;
        str[i].num = nameToNum(name);
        intr[str[i].b].out.push_back(str[i].num);
        intr[str[i].e].in.push_back(str[i].num);
    }
    for (int i = 0; i < numV; ++i)
    {
        car[i].num = i;
        inFile >> p;
        for (int j = 0; j < p; ++j)
        {
            inFile >> name;
            car[i].path.push_back(nameToNum(name));
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
    ss << numI << "\n";
    for (int i = 0; i < numI; ++i)
    {
        ss << intr[i].num << "\n";
        int nonZero = 0;
        for (int sNum : intr[i].in)
        {
            if (str[sNum].time <= 0) continue;
            ++nonZero;
        }
        ss << nonZero << "\n";
        random_shuffle(intr[i].in.begin(), intr[i].in.end());
        for (int sNum : intr[i].in)
        {
            if (str[sNum].time <= 0) continue;
            ss << numToName[sNum] << " " << str[sNum].time << "\n";
        }
    }
    outputStr = ss.str();
}

void findCarNums()
{
    for (int i = 0; i < numV * 0.975; ++i)
    {
        for (int j = 0; j < (int) car[i].path.size() - 1; ++j)
        {
            ++str[car[i].path[j]].cars;
        }
    }

    for (int i = 0; i < numI; ++i)
    {
        for (int sNum : intr[i].in)
        {
            intr[i].cars += str[sNum].cars;
        }
    }
}

bool cmpByTotLen(const Car& a, const Car& b)
{
    return a.totLen < b.totLen;
}

void findCarPathLens()
{
    for (int i = 0; i < numV; ++i)
    {
        car[i].totLen = 0;
        for (int j = 1; j < (int) car[i].path.size(); ++j)
        {
            car[i].totLen += str[car[i].path[j]].len;
            car[i].totLen += intr[str[car[i].path[j]].b].in.size() / 2.5;
        }
    }
    sort(car, car + numV, cmpByTotLen);
}

void solve()
{
    findCarPathLens();
    findCarNums();
    for (int i = 0; i < numI; ++i)
    {
        long long totTime = 1 * intr[i].in.size();
        for (int sNum : intr[i].in)
        {
            if (intr[i].cars == 0)
            {
                str[sNum].time = 1;
                continue;
            }
            //str[sNum].time = str[sNum].cars > 0;
            str[sNum].time = (totTime * str[sNum].cars + intr[i].cars - 1) / intr[i].cars;
            //str[sNum].time = std::min(str[sNum].time, 3);
        }
    }
    generateOutput();
}

int main()
{
    cin >> test;
    inFile.open((tests[test] + ".txt").c_str());
    outFile.open(("sol_" + tests[test] + ".txt").c_str());

    srand(-555);
    input();
    solve();
    output();
    return 0;
}
