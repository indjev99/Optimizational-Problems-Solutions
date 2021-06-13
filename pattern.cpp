#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <numeric>
#include <limits>
#include <vector>
#include <set>
#include <algorithm>
#include <queue>
#include <string>
#include <unordered_map>
#include <iterator>
#include <assert.h>

std::ifstream inF("pattern.in");
std::ofstream outF("pattern.out");

std::mt19937 generator;
int randNum(int lb, int ub)
{
    std::uniform_int_distribution<int> distribution(lb, ub);
    return distribution(generator);
}

typedef unsigned long long ull;

int n;
std::string text;
std::string ansPattern;

void input()
{
    inF >> text;
    n = text.size();
}

void output()
{
    std::cerr << "Length: " << ansPattern.size() << std::endl;
    std::cerr << ansPattern << std::endl;
    outF << ansPattern << std::endl;
}

double timeLimit = 4.8;
std::chrono::high_resolution_clock::time_point startT, currT;
bool haveTime()
{
    using namespace std::chrono;
    currT = high_resolution_clock::now();
    double time = duration_cast<duration<double>>(currT - startT).count();
    return time < timeLimit;
}

const int MIN_LEN2 = 1;
const int MAX_LEN2 = 100;

struct RepPat
{
    std::string str;
    int reps;
};

std::string repSolve(const char* txt, int sz)
{
    if (sz <= 5)
    {
        return std::string(txt, sz);
    }

    std::string res;
    std::vector<RepPat> repSeq;

    for (int i = 0; i < sz; ++i)
    {
        int maxSave = 0;
        int bestLen = 1;
        int bestRep = 1;
        for (int len = MIN_LEN2; len <= MAX_LEN2; ++len)
        {
            int repCnt = 1;
            for (; ; ++repCnt)
            {
                bool match = true;
                for (int j = 0; j < len; ++j)
                {
                    int k = i + len * repCnt + j;
                    if (k >= sz || txt[i + j] != txt[k])
                    {
                        match = false;
                        break;
                    }
                }
                if (!match) break;
            }
            int currSave = (repCnt) * len - len - 3 - std::to_string(repCnt).size();
            if (currSave > maxSave)
            {
                maxSave = currSave;
                bestLen = len;
                bestRep = repCnt;
            }
        }
        repSeq.push_back({std::string(txt + i, bestLen), bestRep});
        i += bestLen * bestRep - 1;
    }

    for (auto& pat : repSeq)
    {
        if (pat.reps > 1)
        {
            res += '[';
            res += pat.str;
            res += ',';
            res += std::to_string(pat.reps);
            res += ']';
        }
        else res += pat.str;;
    }

    return res;
}

#define ALPH 0
#define SPECIAL 1

struct Segment
{
    std::string str;
    int type;
};

std::vector<Segment> segs;

ull BASE = 113;
ull MOD = 5112733823;

std::unordered_map<ull, int> hashCount;
std::unordered_map<ull, std::pair<int, int>> hashLasts;

std::pair<int, std::pair<int, int>> maxCounts(int len)
{
    hashCount.clear();
    hashLasts.clear();

    int maxCount = 0;
    std::pair<int, int> startPos = {-1, -1};
    ull currHash = 0;
    ull maxBase = 1;
    for (int i = 0; i < len - 1; ++i)
    {
        maxBase = maxBase * BASE % MOD;
    }

    for (int t = 0; t < segs.size(); ++t)
    {
        if (segs[t].type != ALPH) continue;
        currHash = 0;
        const std::string& txt = segs[t].str;
        for (int i = 0; i < txt.size(); ++i)
        {
            if (i >= len) currHash = (currHash - txt[i - len] * maxBase) % MOD;
            currHash = (currHash * BASE + txt[i]) % MOD;
            if (i >= len - 1 && hashLasts[currHash] <= std::make_pair(t, i - len + 1))
            {
                hashLasts[currHash] = {t, i + 1};
                ++hashCount[currHash];
                if (hashCount[currHash] > maxCount)
                {
                    maxCount = hashCount[currHash];
                    startPos = {t, i - len + 1};
                }
            }
        }
    }

    return {maxCount, startPos};
}

std::string replaceSubstring(const std::string& txt, int len, const std::string& patTxt, int startPos)
{
    std::string res;
    for (int i = 0; i < txt.size(); ++i)
    {
        bool match = true;
        for (int j = 0; j < len; ++j)
        {
            if (txt[i + j] != patTxt[startPos + j])
            {
                match = false;
                break;
            }
        }
        if (!match) res += txt[i];
        else
        {
            res += '#';
            i += len - 1;
        }
    }
    return res;
}

const int MIN_LEN = 2;
const int MAX_LEN = 20;

void replaceSolve()
{
    int maxSave = 0;
    int bestLen = -1;
    std::pair<int, int> bestStartPos;

    for (int len = MIN_LEN; len <= MAX_LEN; ++len)
    {
        auto pr = maxCounts(len);
        int currSave = pr.first * (len - 1) - 3 - len;
        if (currSave > maxSave)
        {
            maxSave = currSave;
            bestLen = len;
            bestStartPos = pr.second;
        }
    }

    std::cerr << "Best save: " << maxSave << " " << bestLen << std::endl;

    if (maxSave > 0)
    {
        ansPattern = "#=";
        for (int i = 0; i < bestLen; ++i)
        {
            ansPattern += segs[bestStartPos.first].str[bestStartPos.second + i];
        }
        ansPattern += '.';
        for (auto& seg : segs)
        {
            if (seg.type == SPECIAL) ansPattern += seg.str;
            else ansPattern += replaceSubstring(seg.str, bestLen, segs[bestStartPos.first].str, bestStartPos.second);
        }
    }
}

void solve1()
{
    std::cerr << "Solve 1" << std::endl;

    ansPattern = repSolve(text.data(), text.size());
    std::cerr << "Stage 1: " << ansPattern << std::endl;

    std::string currAlph;
    for (int i = 0; i < ansPattern.size(); ++i)
    {
        char c = ansPattern[i];
        if ('a' <= c && c <= 'z')
        {
            currAlph += c;
        }
        else
        {
            if (currAlph.size() > 0) segs.push_back({currAlph, ALPH});
            segs.push_back({std::string(1, c), SPECIAL});
            currAlph = "";
        }
    }
    if (currAlph.size() > 0) segs.push_back({currAlph, ALPH});

    replaceSolve();
    std::cerr << "Stage 2: " << ansPattern << std::endl;
}

void solve2()
{
    std::cerr << "Solve 2" << std::endl;

    ansPattern = text;
    segs.push_back({text, ALPH});
    replaceSolve();
    std::cerr << "Stage 1: " << ansPattern << std::endl;

    int start = 0;
    for (int i = 0; i < ansPattern.size(); ++i)
    {
        if (ansPattern[i] == '.')
        {
            start = i + 1;
            break;
        }
    }

    std::string pref = std::string(ansPattern.data(), start);
    ansPattern = pref + repSolve(ansPattern.data() + start, ansPattern.size() - start);
    std::cerr << "Stage 2: " << ansPattern << std::endl;
}

void solve()
{
    solve1();
    std::string currAns = ansPattern;

    segs.clear();
    solve2();

    if (currAns.size() < ansPattern.size()) ansPattern = currAns;
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
