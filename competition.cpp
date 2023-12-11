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

using ull = unsigned long long;

constexpr bool DEBUG = true;
constexpr double TIME_MULT = 1;

const double SOFT_MAX_TIME = 4.0;
const double HARD_MAX_TIME = 4.5;

bool DEBUG_PRINT = false;

std::ifstream inF("competition.in");
std::ofstream outF("competition.out");

std::mt19937 generator;

ull randNum(ull lb, ull ub)
{
    std::uniform_int_distribution<ull> distribution(lb, ub - 1);
    return distribution(generator);
}

double randReal(double lb, double ub)
{
    std::uniform_real_distribution<double> distribution(lb, ub);
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

const int MAX_NUM_BITS = 63;
const ull MAX_NUM = 1ULL << MAX_NUM_BITS;

const int MAX_L = 150;

const int ADD = 1;
const int SUB = 2;
const int MUL = 3;
const int DIV = 4;
const int MOD = 5;
const int AND = 6;
const int OR = 7;
const int XOR = 8;

const int NUM_REASONS = 4;
const int R_SETUP = 0;
const int R_SHIFT = 1;
const int R_PATTERN = 2;
const int R_UNKNOWN = 3;

const std::string REASON_NAMES[NUM_REASONS] = {
    "SETUP",
    "SHIFT",
    "PTTRN",
    "UNKNW"
};

struct Op
{
    int type;
    int j, k;
    int reason;
};

int t, currTest;

ull n;

std::vector<Op> ans;
std::vector<ull> nums;

std::vector<Op> bestAns;
std::vector<ull> bestNums;

ull totalL = 0;
ull totalScore = 0;

ull totalPerReason[NUM_REASONS];

void input()
{
    inF >> n;
}

void output()
{
    if (DEBUG_PRINT)
    {
        std::cerr << bestAns.size() << "\n";
        for (ull num : bestNums)
        {
            std::cerr << num << " ";
        }
        std::cerr << std::endl;
    }

    totalL += bestAns.size();
    totalScore += bestAns.size() * bestAns.size();

    outF << bestAns.size() << "\n";
    for (const Op& op : bestAns)
    {
        outF << op.type << " " << op.j << " " << op.k << "\n";
        totalPerReason[op.reason]++;
    }
}

void updateAns()
{
    assert((int) nums.size() >= 1);
    assert(nums.size() == ans.size() + 1);
    assert(nums.front() == 1);
    assert(nums.back() == n);

    for (int i = 0; i < (int) ans.size(); ++i)
    {
        int type = ans[i].type;
        int j = ans[i].j;
        int k = ans[i].k;
    
        assert(j >= 0 && j <= i);
        assert(k >= 0 && k <= i);

        ull res = 0;
        switch (type)
        {
        case ADD:
            res = nums[j] + nums[k];
            break;
        case SUB:
            res = nums[j] - nums[k];
            break;
        case MUL:
            res = nums[j] * nums[k];
            break;
        case DIV:
            res = nums[j] / nums[k];
            break;
        case MOD:
            res = nums[j] % nums[k];
            break;
        case AND:
            res = nums[j] & nums[k];
            break;
        case OR:
            res = nums[j] | nums[k];
            break;
        case XOR:
            res = nums[j] ^ nums[k];
            break;
        default:
            assert(false);
        }

        assert(nums[i + 1] == res);
    }

    if (ans.size() < bestAns.size() || bestNums.empty())
    {
        bestAns = ans;
        bestNums = nums;
    }
}

void bestSetup()
{
    bestAns.clear();
    bestNums.clear();
}

void setup()
{
    ans.clear();
    nums.clear();
    nums.push_back(1);
}

void addOp(int type, int j, int k, int reason = R_UNKNOWN)
{
    assert(j >= 0 && j < (int) nums.size());
    assert(k >= 0 && k < (int) nums.size());
    ull res = 0;
    switch (type)
    {
    case ADD:
        res = nums[j] + nums[k];
        break;
    case SUB:
        res = nums[j] - nums[k];
        break;
    case MUL:
        res = nums[j] * nums[k];
        break;
    case DIV:
        res = nums[j] / nums[k];
        break;
    case MOD:
        res = nums[j] % nums[k];
        break;
    case AND:
        res = nums[j] & nums[k];
        break;
    case OR:
        res = nums[j] | nums[k];
        break;
    case XOR:
        res = nums[j] ^ nums[k];
        break;
    default:
        assert(false);
    }
    ans.push_back({type, j, k, reason});
    nums.push_back(res);
}

bool checkBit(int i)
{
    return (1ULL << i) & n;
}

void solveDummy()
{
    setup();

    while (2ULL * nums.back() <= n)
    {
        addOp(ADD, nums.size() - 1, nums.size() - 1, R_SETUP);
    }

    int lastIdx = nums.size();

    for (int i = 0; i < lastIdx; ++i)
    {
        if (checkBit(i))
        {
            addOp(ADD, nums.size() - 1, i, R_PATTERN);
        }
    }

    updateAns();
}

void solveDI()
{
    setup();

    int lastPow = 0;
    for (int i = 0; i < MAX_NUM_BITS; ++i)
    {
        if (checkBit(i)) lastPow = i;
    }

    for (int i = lastPow - 1; i >= 0; --i)
    {
        addOp(ADD, nums.size() - 1, nums.size() - 1, R_SHIFT);
        if (checkBit(i)) addOp(ADD, nums.size() - 1, 0, R_PATTERN);
    }

    updateAns();
}

void solveDID()
{
    setup();

    int lastPow = 0;
    for (int i = 0; i < MAX_NUM_BITS; ++i)
    {
        if (checkBit(i)) lastPow = i;
    }

    bool flipped = false;

    if (n != MAX_NUM && lastPow >= 2 && checkBit(lastPow - 1) && checkBit(lastPow - 2))
    {
        ++lastPow;
        flipped = true;
    }

    for (int i = lastPow - 1; i >= 0; --i)
    {
        addOp(ADD, nums.size() - 1, nums.size() - 1, R_SHIFT);
        if (!flipped)
        {
            if (checkBit(i)) addOp(ADD, nums.size() - 1, 0, R_PATTERN);
            else if (i >= 2 && checkBit(i - 1) && checkBit(i - 2))
            {
                addOp(ADD, nums.size() - 1, 0, R_PATTERN);
                flipped = true;
            }
        }
        else
        {
            assert(checkBit(i));
            if (checkBit(i) && (i == 0 || !checkBit(i - 1)))
            {
                addOp(SUB, nums.size() - 1, 0, R_PATTERN);
                flipped = false;
            }
        }
    }

    updateAns();
}

void solveSID(int maxShift)
{
    assert(maxShift >= 1);

    setup();

    for (int i = 1; i <= maxShift; ++i)
    {
        addOp(ADD, nums.size() - 1, nums.size() - 1, R_SETUP);
    }

    int currShift = maxShift;

    int lastPow = 0;
    for (int i = 0; i < MAX_NUM_BITS; ++i)
    {
        if (checkBit(i)) lastPow = i;
    }

    bool flipped = false;

    if (n != MAX_NUM && lastPow >= 2 && checkBit(lastPow - 1) && checkBit(lastPow - 2))
    {
        ++lastPow;
        flipped = true;
    }

    if (currShift > lastPow) return;

    for (int i = lastPow - 1; i >= 0; --i)
    {
        if (currShift == 0)
        {
            int shift = std::min(maxShift, i + 1);
            addOp(MUL, nums.size() - 1, shift, R_SHIFT);
            currShift = shift;
        }

        --currShift;

        if (!flipped)
        {
            if (checkBit(i)) addOp(ADD, nums.size() - 1, currShift, R_PATTERN);
            else if (i >= 2 && checkBit(i - 1) && checkBit(i - 2))
            {
                addOp(ADD, nums.size() - 1, currShift, R_PATTERN);
                flipped = true;
            }
        }
        else
        {
            assert(checkBit(i));
            if (checkBit(i) && (i == 0 || !checkBit(i - 1)))
            {
                addOp(SUB, nums.size() - 1, currShift, R_PATTERN);
                flipped = false;
            }
        }
    }

    assert(currShift == 0);

    updateAns();
}

void solveSIDDP(int maxShift)
{
    assert(maxShift >= 1);

    setup();

    for (int i = 1; i <= maxShift; ++i)
    {
        addOp(ADD, nums.size() - 1, nums.size() - 1, R_SETUP);
    }

    int currShift = maxShift;

    int lastPow = 0;
    for (int i = 0; i < MAX_NUM_BITS; ++i)
    {
        if (checkBit(i)) lastPow = i;
    }

    const int DP_NOP = 0;
    const int DP_ADD = 1;
    const int DP_SUB = -1;

    std::vector<std::array<int, 2>> dp(lastPow + 1);
    std::vector<std::array<int, 2>> dpHow(lastPow + 1);

    dp[lastPow][false] = (lastPow + maxShift - 1) / maxShift;
    dpHow[lastPow][false] = DP_NOP;

    dp[lastPow][true] = (lastPow + maxShift) / maxShift;
    dpHow[lastPow][true] = DP_NOP;

    for (int i = lastPow - 1; i >= 0; --i)
    {
        if (checkBit(i))
        {
            dp[i][true] = dp[i + 1][true]; // nop
            dpHow[i][true] = DP_NOP;

            dp[i][false] = std::min(
                1 + dp[i + 1][false], // add
                1 + dp[i + 1][true] // sub
            );
            if (dp[i][false] == 1 + dp[i + 1][false]) dpHow[i][false] = DP_ADD;
            else if (dp[i][false] == 1 + dp[i + 1][true]) dpHow[i][false] = DP_SUB;
            else assert(false);
        }
        else
        {
            dp[i][false] = dp[i + 1][false]; // nop
            dpHow[i][false] = DP_NOP;

            dp[i][true] = std::min(
                1 + dp[i + 1][false], // add
                1 + dp[i + 1][true] // sub
            );
            if (dp[i][true] == 1 + dp[i + 1][false]) dpHow[i][true] = DP_ADD;
            else if (dp[i][true] == 1 + dp[i + 1][true]) dpHow[i][true] = DP_SUB;
            else assert(false);
        }
    }

    int pos = 0;
    bool flipped = false;

    std::vector<int> hows;

    while (pos < lastPow)
    {
        int currHow = dpHow[pos][flipped];
        hows.push_back(currHow);
        if (checkBit(pos))
        {
            if (!flipped && currHow == DP_SUB) flipped = true;
        }
        else
        {
            if (flipped && currHow == DP_ADD) flipped = false;
        }
        ++pos;
    }

    if (flipped)
    {
        ++lastPow;
        hows.push_back(DP_NOP);
    }

    if (currShift > lastPow) return;

    for (int i = lastPow - 1; i >= 0; --i)
    {
        if (currShift == 0)
        {
            int shift = std::min(maxShift, i + 1);
            addOp(MUL, nums.size() - 1, shift, R_SHIFT);
            currShift = shift;
        }

        --currShift;

        if (hows[i] == DP_ADD) addOp(ADD, nums.size() - 1, currShift, R_PATTERN);
        else if (hows[i] == DP_SUB) addOp(SUB, nums.size() - 1, currShift, R_PATTERN);
    }

    assert(currShift == 0);

    updateAns();
}

void removeNum(int i)
{
    assert(i > 0);
    for (int j = i; j < (int) ans.size(); ++j)
    {
        if (ans[j].j == i || ans[j].k == i)
        {
            std::cerr << "Removing: " << i << " = " << nums[i] << std::endl;
            std::cerr << "Invalid due to: " << j << " = " << nums[j + 1] << " with " << ans[j].type << " "  << ans[j].j << " = " << nums[ans[j].j] << " and " << ans[j].k << " = " << nums[ans[j].k] << std::endl;
        }

        assert(ans[j].j != i);
        assert(ans[j].k != i);

        if (ans[j].j > i) --ans[j].j;
        if (ans[j].k > i) --ans[j].k;
    }
    ans.erase(ans.begin() + i - 1);
    nums.erase(nums.begin() + i);
}

void solveSIDDP2(int maxShift)
{
    assert(maxShift >= 1);

    setup();

    std::vector<bool> usedShifts(maxShift + 1, false);

    usedShifts[0] = true;
    usedShifts[maxShift] = true;

    for (int i = 1; i <= maxShift; ++i)
    {
        addOp(ADD, nums.size() - 1, nums.size() - 1, R_SETUP);
    }

    int currShift = maxShift;

    int lastPow = 0;
    for (int i = 0; i < MAX_NUM_BITS; ++i)
    {
        if (checkBit(i)) lastPow = i;
    }

    const int DP_NOP = 0;
    const int DP_ADD = 1;
    const int DP_SUB = -1;

    std::vector<std::array<int, 2>> dp(lastPow + 1);
    std::vector<std::array<int, 2>> dpHow(lastPow + 1);

    dp[lastPow][false] = (lastPow + maxShift - 1) / maxShift;
    dpHow[lastPow][false] = DP_NOP;

    dp[lastPow][true] = (lastPow + maxShift) / maxShift;
    dpHow[lastPow][true] = DP_NOP;

    for (int i = lastPow - 1; i >= 0; --i)
    {
        if (checkBit(i))
        {
            dp[i][true] = dp[i + 1][true]; // nop
            dpHow[i][true] = DP_NOP;

            dp[i][false] = std::min(
                1 + dp[i + 1][false], // add
                1 + dp[i + 1][true] // sub
            );
            if (dp[i][false] == 1 + dp[i + 1][false]) dpHow[i][false] = DP_ADD;
            else if (dp[i][false] == 1 + dp[i + 1][true]) dpHow[i][false] = DP_SUB;
            else assert(false);
        }
        else
        {
            dp[i][false] = dp[i + 1][false]; // nop
            dpHow[i][false] = DP_NOP;

            dp[i][true] = std::min(
                1 + dp[i + 1][false], // add
                1 + dp[i + 1][true] // sub
            );
            if (dp[i][true] == 1 + dp[i + 1][false]) dpHow[i][true] = DP_ADD;
            else if (dp[i][true] == 1 + dp[i + 1][true]) dpHow[i][true] = DP_SUB;
            else assert(false);
        }
    }

    int pos = 0;
    bool flipped = false;

    std::vector<int> hows;

    while (pos < lastPow)
    {
        int currHow = dpHow[pos][flipped];
        hows.push_back(currHow);
        if (checkBit(pos))
        {
            if (!flipped && currHow == DP_SUB) flipped = true;
        }
        else
        {
            if (flipped && currHow == DP_ADD) flipped = false;
        }
        ++pos;
    }

    if (flipped)
    {
        ++lastPow;
        hows.push_back(DP_NOP);
    }

    if (currShift > lastPow) return;

    for (int i = lastPow - 1; i >= 0; --i)
    {
        if (currShift == 0)
        {
            int shift = std::min(maxShift, i + 1);
            addOp(MUL, nums.size() - 1, shift, R_SHIFT);
            currShift = shift;
            usedShifts[shift] = true;

            // std::cerr << std::endl;
            // for (int i = 0; i < maxShift - shift; ++i) std::cerr << ".";
        }

        --currShift;

        // if (hows[i] == DP_NOP) std::cerr << ".";
        // else if (hows[i] == DP_ADD) std::cerr << "+";
        // else if (hows[i] == DP_SUB) std::cerr << "-";

        if (hows[i] != DP_NOP) usedShifts[currShift] = true;

        if (hows[i] == DP_ADD) addOp(ADD, nums.size() - 1, currShift, R_PATTERN);
        else if (hows[i] == DP_SUB) addOp(SUB, nums.size() - 1, currShift, R_PATTERN);
    }

    // std::cerr << std::endl;
    // std::cerr << std::endl;

    assert(currShift == 0);

    int numRemoved = 0;

    for (int i = maxShift; i >= 3; --i)
    {
        if (usedShifts[i]) continue;
        bool succ = true;
        for (int j = i + 1; j <= maxShift - numRemoved; ++j)
        {
            if (ans[j - 1].j != i && ans[j - 1].k != i) continue;
            succ = false;
            for (int k1 = 1; k1 < j && !succ; k1++)
            {
                if (k1 == i) continue;
                for (int k2 = 1; k2 < j && !succ; k2++)
                {
                    if (k2 == i) continue;
                    if (nums[k1] * nums[k2] == nums[j])
                    {
                        ans[j - 1] = {MUL, k1, k2, R_SETUP};
                        succ = true;
                    }
                }
            }
            if (!succ)
            {
                // std::cerr << "Blocked: " << nums[i] << " by " << nums[j] << std::endl;
                break;
            }
        }
        if (!succ) continue;
        // std::cerr << "Removing: " << nums[i] << std::endl;
        removeNum(i);
        ++numRemoved;
        // break;
    }

    updateAns();
}

void solve()
{
    bestSetup();
    // solveDI();
    // solveDID();
    // for (int maxShift = 1; maxShift <= 8; ++maxShift)
    // {
    //     solveSID(maxShift);  
    // }
    // for (int maxShift = 1; maxShift <= 8; ++maxShift)
    // {
    //     solveSIDDP(maxShift);  
    // }
    for (int maxShift = 1; maxShift <= 30; ++maxShift)
    {
        solveSIDDP2(maxShift);  
    }
}

int main()
{
    generator.seed(0);
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    startT = std::chrono::high_resolution_clock::now();

    inF >> t;

    // t = 1;
    // DEBUG_PRINT = true;

    for (currTest = 1; currTest <= t; ++currTest)
    {
        input();
        solve();
        output();
    }

    long double ml = totalL;
    ml /= t;

    long double smsl = totalScore;
    smsl /= t;
    smsl = sqrtl(smsl);

    std::cerr << "ML: " << ml << std::endl;
    std::cerr << "SMSL: " << smsl << std::endl;
    std::cerr << std::endl;

    for (int i = 0; i < NUM_REASONS; ++i)
    {
        long double mpr = totalPerReason[i];
        mpr /= t;
        std::cerr << REASON_NAMES[i] << ": " << mpr << std::endl;
    }

    std::cerr << std::endl;
    std::cerr << "Total score: " << totalScore << std::endl;
    std::cerr << "Time passed: " << timePassed() << std::endl;

    return 0;
}
