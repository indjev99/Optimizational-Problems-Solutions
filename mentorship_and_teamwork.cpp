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
ifstream bestFile;

long long bestTotalScore = -1;
long long totalScore;

std::string outputStr;

std::mt19937 generator;
int randNum(int lb, int ub)
{
    std::uniform_int_distribution<int> distribution(lb, ub - 1);
    return distribution(generator);
}

const int MAX_CP = 1e5;
const int INF = 1e9;

unordered_map<string, int> skillToNumMap;
vector<string> numToSkill;

int skillToNum(const string& name)
{
    if (skillToNumMap.find(name) == skillToNumMap.end())
    {
        skillToNumMap[name] = numToSkill.size();
        numToSkill.push_back(name);
    }
    return skillToNumMap[name];
}

struct Contr
{
    std::string name;
    std::vector<int> levels; // skill -> level

    int freeFrom;
    std::vector<int> levelsRun; // skill -> level

    bool used;
};

struct Proj
{
    std::string name;
    int days;
    int score;
    int bestBefore;
    std::vector<std::pair<int, int>> roles; // skill level

    int startDay;
    std::vector<int> people;

    int scoreRun()
    {
        if (startDay == -1) return 0;

        int endDay = startDay + days;
        if (endDay <= bestBefore) return score;
        return std::max(score - (endDay - bestBefore), 0);
    }
};

int c, p;

Contr contrs[MAX_CP];
Proj projs[MAX_CP];

void input()
{
    inFile >> c >> p;
    int numSkills;
    std::string skill;
    int skillNum;
    int level;

    for (int i = 0; i < c; ++i)
    {
        inFile >> contrs[i].name >> numSkills;
        for (int j = 0; j < numSkills; ++j)
        {
            inFile >> skill >> level;
            skillNum = skillToNum(skill);
            if ((int) contrs[i].levels.size() <= skillNum) contrs[i].levels.resize(skillNum + 1, 0);
            contrs[i].levels[skillNum] = level;
        }
    }

    for (int i = 0; i < p; ++i)
    {
        inFile >> projs[i].name >> projs[i].days >> projs[i].score >> projs[i].bestBefore >> numSkills;
        for (int j = 0; j < numSkills; ++j)
        {
            inFile >> skill >> level;
            skillNum = skillToNum(skill);
            projs[i].roles.push_back({skillNum, level});
        }
    }

    for (int i = 0; i < c; ++i)
    {
        contrs[i].levels.resize(numToSkill.size(), 0);
    }

    std::cerr << "Number of skills: " << numToSkill.size() << std::endl;
}

int bestOrderEnd;

void loadBestOrder()
{
    std::string name;

    bestFile >> bestOrderEnd;

    for (int i = 0; i < bestOrderEnd; ++i)
    {
        bestFile >> name;
        for (int j = i; ; ++j)
        {
            if (projs[j].name == name)
            {
                std::swap(projs[i], projs[j]);
                break;
            }
        }

        for (int j = 0; j < (int) projs[i].roles.size(); ++j)
        {
            bestFile >> name;
        }
    }
}

void output()
{
    outFile << outputStr;
}

Contr bestContrs[MAX_CP];

void generateOutput()
{
    if (totalScore <= bestTotalScore) return;

    std::copy(contrs, contrs + c, bestContrs);

    bestTotalScore = totalScore;
    std::cerr << "New best score: " << bestTotalScore << std::endl;;

    ostringstream ss;

    int e = 0;
    for (int i = 0; i < p; ++i)
    {
        if (projs[i].startDay >= 0) ++e;
    }

    ss << e << "\n";
    for (int i = 0; i < p; ++i)
    {
        if (projs[i].startDay == -1) continue;

        ss << projs[i].name << "\n";
        for (int person : projs[i].people)
        {
            ss << contrs[person].name << " ";
        }
        ss << "\n";
    }

    outputStr = ss.str();
}

bool cmpProj(const Proj& left, const Proj& right)
{
    // return left.bestBefore < right.bestBefore;
    return left.days < right.days;
    // return left.score > right.score;
    // return left.roles.size() < right.roles.size();
}

const int MAX_SKILLS = 1000;
const int MAX_LEVEL = 20;

std::vector<int> peopleAtLevel[MAX_SKILLS][MAX_LEVEL];

void solveOne(bool shuffle = false)
{
    if (shuffle) std::shuffle(contrs, contrs + c, generator);

    for (int i = 0; i < MAX_SKILLS; ++i)
    {
        for (int j = 0; j < MAX_LEVEL; ++j)
        {
            peopleAtLevel[i][j].clear();
        }
    }

    totalScore = 0;

    for (int i = 0; i < c; ++i)
    {
        contrs[i].freeFrom = 0;
        contrs[i].levelsRun = contrs[i].levels;

        for (int j = 0; j < (int) numToSkill.size(); ++j)
        {
            for (int k = 0; k <= contrs[i].levels[j]; ++k)
            {
                peopleAtLevel[j][k].push_back(i);
            }
        }
    }

    for (int i = 0; i < p; ++i)
    {
        projs[i].startDay = -1;
        projs[i].people.clear();
        projs[i].people.resize(projs[i].roles.size(), -1);
    }

    for (int i = 0; i < p; ++i)
    {
        Proj& proj = projs[i];

        std::vector<int> roleOrd(proj.roles.size());
        std::iota(roleOrd.begin(), roleOrd.end(), 0);
        if (shuffle) std::shuffle(roleOrd.begin(), roleOrd.end(), generator);

        std::vector<int> maxSkillLevels(numToSkill.size(), 0);

        proj.startDay = 0;

        for (int k = 0; k < c; ++k)
        {
            contrs[k].used = false;
        }

        for (int j : roleOrd)
        {
            int skill = proj.roles[j].first;
            int level = proj.roles[j].second;

            int bestK = -1;
            int bestDay = INF;
            int bestLevel = INF;

            if (maxSkillLevels[skill] >= level) --level;

            for (int k : peopleAtLevel[skill][level])
            {
                Contr& contr = contrs[k];

                if (contr.used) continue;

                int day = std::max(contr.freeFrom, proj.startDay);
                if (day < bestDay || (day == bestDay && contr.levelsRun[skill] < bestLevel))
                // if (contr.levelsRun[skill] < bestLevel || (contr.levelsRun[skill] == bestLevel && day < bestDay))
                {
                    bestK = k;
                    bestDay = day;
                    bestLevel = level;
                }
            }

            if (bestK == -1)
            {
                proj.startDay = -1;
                break;
            }

            proj.startDay = bestDay;
            proj.people[j] = bestK;
            contrs[bestK].used = true;

            for (int t = 0; t < (int) numToSkill.size(); ++t)
            {
                maxSkillLevels[t] = std::max(contrs[bestK].levels[t], maxSkillLevels[t]);
            }
        }

        if (proj.scoreRun() <= 0) proj.startDay = -1;

        if (proj.startDay == -1) continue;

        for (int j1 = 0; j1 < (int) roleOrd.size(); ++j1)
        {
            int& p1 = proj.people[j1];
            int skill1 = proj.roles[j1].first;
            int level1 = proj.roles[j1].second;
            for (int j2 = j1 + 1; j2 < (int) roleOrd.size(); ++j2)
            {
                int& p2 = proj.people[j2];
                int skill2 = proj.roles[j2].first;
                int level2 = proj.roles[j2].second;
                if (contrs[p1].levelsRun[skill2] >= level2 - 1 && contrs[p2].levelsRun[skill1] >= level1 - 1 &&
                    (contrs[p1].levelsRun[skill2] <= level2 || contrs[p2].levelsRun[skill1] <= level1) &&
                    !(contrs[p1].levelsRun[skill1] <= level1 && contrs[p2].levelsRun[skill2] <= level2))
                {
                    std::swap(p1, p2);
                }
            }
        }

        for (int j : roleOrd)
        {
            Contr& contr = contrs[proj.people[j]];
            contr.freeFrom = proj.startDay + proj.days;

            int skill = proj.roles[j].first;
            int level = proj.roles[j].second;
            if (contr.levelsRun[skill] <= level)
            {
                ++contr.levelsRun[skill];
                peopleAtLevel[skill][contr.levelsRun[skill]].push_back(proj.people[j]);
            }
        }

        totalScore += projs[i].scoreRun();
    }

    generateOutput();
}

int SUB_TRIALS = 1;

long long testOrder()
{
    long long sum = 0;
    for (int i = 0; i < SUB_TRIALS; ++i)
    {
        solveOne();
        sum += totalScore;
    }
    return sum / SUB_TRIALS;
}

void trySwap(int i, int j)
{
    long long curr = bestTotalScore;
    std::swap(projs[i], projs[j]);
    long long other = testOrder();
    // std::cerr << i << " " << j << " : " << curr << " vs. " << other << std::endl;
    if (other < curr) std::swap(projs[i], projs[j]);
    else if (other > curr) std::cerr << " Swapping: " << i << " " << j << " : " << curr << " vs. " << other << std::endl;
}

int TRIALS;

void solve()
{
    std::shuffle(contrs, contrs + c, generator);

    // to act as local optimizer
    // loadBestOrder();

    std::sort(projs + bestOrderEnd, projs + p, cmpProj);

    solveOne();

    for (int t = 0; t < TRIALS; ++t)
    {
        int i, j;
        do
        {
            i = randNum(0, p);
            j = randNum(0, p);
        }
        while (i == j);
        trySwap(i, j);
    }

    for (int t = 0; t < TRIALS; ++t)
    {
        solveOne(true);
    }

    std::copy(bestContrs, bestContrs + c, contrs);
}

int main()
{
    cin >> test >> TRIALS;
    inFile.open((tests[test] + ".txt").c_str());
    outFile.open(("sol_" + tests[test] + ".txt").c_str());
    bestFile.open(("best_" + tests[test] + ".txt").c_str());

    generator.seed(20);

    input();
    solve();
    output();

    return 0;
}
