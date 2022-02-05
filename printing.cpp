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
#include <iterator>
#include <assert.h>

const bool DEBUG = false;

std::ifstream inF("printing.in");
std::ofstream outF("printing.out");

std::mt19937 generator;
int randNum(int lb, int ub)
{
    std::uniform_int_distribution<int> distribution(lb, ub - 1);
    return distribution(generator);
}

double maxTimeLimit = 4.5;
std::chrono::high_resolution_clock::time_point startT, currT;

double timeLeft(double timeLimit)
{
    using namespace std::chrono;
    currT = high_resolution_clock::now();
    double time = duration_cast<duration<double>>(currT - startT).count();
    return timeLimit - time;
}

typedef unsigned long long ull;

const int NUM_LETS = 27;
const int MAX_N = 1e4;
const ull INF = ULLONG_MAX;

int n;
int rotLeftCost, rotRightCost, basePrintCost, changeCost, stackCost;
int numSeqCosts[NUM_LETS + 1];

std::string textStr;

int printCosts[MAX_N + 1];
int text[MAX_N];

int numUnique;
int unique[NUM_LETS];

int numNonUsed;
int nonUsed[NUM_LETS];

int optNumSeqCosts[MAX_N + 1];
int optNumSeqs[MAX_N + 1];

void input()
{
    inF >> n;
    inF >> textStr;
    inF >> rotLeftCost >> rotRightCost >> basePrintCost >> changeCost >> stackCost;
    for (int i = 1; i <= NUM_LETS; ++i)
    {
        inF >> numSeqCosts[i];
    }
}

void init()
{
    printCosts[0] = 0;
    for (int i = 1; i <= n; ++i)
    {
        printCosts[i] = printCosts[i - 1] + basePrintCost / i;
    }

    numUnique = 0;

    for (int i = 0; i < n; ++i)
    {
        if (textStr[i] == '_') text[i] = NUM_LETS - 1;
        else text[i] = textStr[i] - 'a';

        if (std::find(unique, unique + numUnique, text[i]) == unique + numUnique)
        {
            unique[numUnique++] = text[i];
        }
    }

    numNonUsed = 0;

    for (int i = 0; i < NUM_LETS; ++i)
    {
        if (std::find(unique, unique + numUnique, i) == unique + numUnique)
        {
            nonUsed[numNonUsed++] = i;
        }
    }

    for (int i = 1; i <= numUnique; ++i)
    {
        optNumSeqCosts[i] = numSeqCosts[i];
        optNumSeqs[i] = i;

        for (int j = 1; j <= numNonUsed; ++j)
        {
            if (numSeqCosts[i + j] < optNumSeqCosts[i])
            {
                optNumSeqCosts[i] = numSeqCosts[i + j];
                optNumSeqs[i] = i + j;
            }
        }
    }

}

struct Sequence
{
    int size;
    int lets[NUM_LETS];

    int inPos[NUM_LETS];
};

#define T_CHANGEROT 0
#define T_ONESTACKU 1
#define T_ONESTACKS 2
#define T_MAXPRINTC 3
#define T_MAXPRINTR 4

const std::string TYPE_NAMES[] = {"T_CHANGEROT", "T_ONESTACKU", "T_ONESTACKS", "T_MAXPRINTC", "T_MAXPRINTR"};

struct Solution
{
    int type;

    int numSeqs;
    Sequence seqs[NUM_LETS];

    int inSeq[NUM_LETS];

    ull cost = INF;

    bool operator<(const Solution& other) const
    {
        return cost < other.cost;
    }
};

#define ROT_LEFT 0
#define ROT_RIGHT 1
#define PRINT 2
#define CHANGE 3
#define PUSH 4
#define POP 5

void fixSolution(Solution& sol)
{
    int left = 0;
    int right = NUM_LETS - 1;
    while (left <= right)
    {
        if (sol.seqs[left].size > 0) ++left;
        else if (sol.seqs[right].size == 0) --right;
        else
        {
            std::swap(sol.seqs[left], sol.seqs[right]);
            ++left;
            --right;
        }
    }

    sol.numSeqs = left;
    
    int startIn = -1;
    for (int i = 0; i < sol.numSeqs && startIn == -1; ++i)
    {
        Sequence& seq = sol.seqs[i];
        for (int j = 0; j < seq.size && startIn == -1; ++j)
        {
            if (seq.lets[j] == text[0])
            {
                startIn = i;
            }
        }
    }

    std::swap(sol.seqs[0], sol.seqs[startIn]);

    for (int i = 0; i < sol.numSeqs; ++i)
    {
        Sequence& seq = sol.seqs[i];
        for (int j = 0; j < seq.size; ++j)
        {
            seq.inPos[seq.lets[j]] = j;
            sol.inSeq[seq.lets[j]] = i;
        }
    }
}

ull dpMPC[MAX_N + 1];
int dpMPCPrev[MAX_N + 1];
int dpMPCNext[MAX_N + 1];
int dpMPCNumLetsPrev[MAX_N + 1];
int dpMPCNumLets[MAX_N + 1];

void findOptMPC()
{
    dpMPC[0] = optNumSeqCosts[numUnique];
    for (int i = 1; i <= n; ++i)
    {
        dpMPC[i] = INF;
    }

    for (int i = 0; i < n; ++i)
    {
        int numLets = 0;
        bool repeat = false;
        for (int j = i; j < n; ++j)
        {
            if (!repeat && std::find(text + i, text + j, text[j]) == text + j) ++numLets;
            else if (text[i + (j - i) % numLets] == text[j]) repeat = true;
            else break;

            ull currCost = dpMPC[i];
            if (i > 0) currCost += changeCost;
            currCost += (numLets - 1) * (changeCost + 2 * stackCost);
            currCost += printCosts[j - i + 1];
            if (j < n - 1) currCost += (numLets - 1) * (changeCost + 2 * stackCost);

            if (currCost < dpMPC[j + 1])
            {
                dpMPC[j + 1] = currCost;
                dpMPCPrev[j + 1] = i;
                dpMPCNumLetsPrev[j + 1] = numLets;
            }
        }
    }

    int curr = n;
    while (curr > 0)
    {
        int prev = dpMPCPrev[curr];
        dpMPCNext[prev] = curr;
        dpMPCNumLets[prev] = dpMPCNumLetsPrev[curr];
        curr = prev;
    }

    std::cerr << "MPC cost: " << dpMPC[n] << std::endl;
}

double expRotCost[NUM_LETS + 1];

double dpMPR[MAX_N + 1];
int dpMPRPrev[MAX_N + 1];
int dpMPRNext[MAX_N + 1];

void findOptMPR()
{
    for (int i = 1; i <= numUnique; ++i)
    {
        expRotCost[i] = 0;
        for (int j = 0; j < i; ++j)
        {
            int leftCost = j * rotLeftCost;
            int rightCost = (i - j) * rotRightCost;
            int pushpopCost = 2 * j * stackCost + rotLeftCost + rotRightCost;
            expRotCost[i] += std::min(std::min(leftCost, rightCost), pushpopCost);
        }
        expRotCost[i] /= i;
    }

    dpMPR[0] = optNumSeqCosts[1];
    for (int i = 1; i <= n; ++i)
    {
        dpMPR[i] = INF;
    }

    for (int i = 0; i < n; ++i)
    {
        double accumCost = dpMPR[i];
        for (int j = i; j < n; ++j)
        {
            if (std::find(text + i, text + j, text[j]) != text + j) break;

            accumCost += expRotCost[numUnique - (j - i)];
            accumCost += 2 * stackCost;

            double currCost = accumCost + printCosts[j - i + 1];

            if (currCost < dpMPR[j + 1])
            {
                dpMPR[j + 1] = currCost;
                dpMPRPrev[j + 1] = i;
            }
        }
    }

    int curr = n;
    while (curr > 0)
    {
        int prev = dpMPRPrev[curr];
        dpMPRNext[prev] = curr;
        curr = prev;
    }

    std::cerr << "MPR cost: " << (ull) round(dpMPR[n]) << std::endl;
}

struct Instruction
{
    int code;
    int arg;
};

template<bool genInstrs>
std::vector<Instruction> findInstructionsMPC(Solution& sol)
{
    sol.cost = dpMPC[n];

    int inSeq[NUM_LETS];
    std::copy(sol.inSeq, sol.inSeq + NUM_LETS, inSeq);
    std::vector<int> freeSeqs;

    std::vector<Instruction> instrs;

    if constexpr (!genInstrs) return instrs;

    int curr = 0;

    while (true)
    {
        int next = dpMPCNext[curr];
        int numLets = dpMPCNumLets[curr];

        if (curr > 0) instrs.push_back({CHANGE, inSeq[text[curr]]});

        for (int i = 0; i < numLets - 1; ++i)
        {
            instrs.push_back({PUSH});
            freeSeqs.push_back(inSeq[text[curr + i]]);
            instrs.push_back({CHANGE, inSeq[text[curr + i + 1]]});
        }

        for (int i = 0; i < numLets - 1; ++i)
        {
            instrs.push_back({POP});
        }

        instrs.push_back({PRINT, next - curr});

        if (next == n) break;

        for (int i = 0; i < numLets - 1; ++i)
        {
            instrs.push_back({PUSH});
        }

        inSeq[text[next - 1]] = inSeq[text[curr + numLets - 1]];

        for (int i = 0; i < numLets - 1; ++i)
        {
            instrs.push_back({CHANGE, freeSeqs.back()});
            instrs.push_back({POP});
            inSeq[text[next - i - 2]] = freeSeqs.back();
            freeSeqs.pop_back();
        }

        curr = next;
    }

    return instrs;
}

template<bool genInstrs>
std::vector<Instruction> findInstructionsMPR(Solution& sol)
{
    int stackSize = 0;
    int stack[NUM_LETS];
    std::copy(sol.seqs[0].lets, sol.seqs[0].lets + numUnique, stack);

    int lastPush = 0;

    std::vector<Instruction> instrs;

    sol.cost = optNumSeqCosts[sol.numSeqs];

    int curr = 0;

    while (curr < n)
    {
        int next = dpMPRNext[curr];

        for (int i = curr; i < next; ++i)
        {
            if (stack[stackSize] != text[i])
            {
                int pos = -1;
                for (int j = stackSize; j < numUnique; ++j)
                {
                    if (stack[j] == text[i])
                    {
                        pos = j;
                        break;
                    }
                }

                int leftMoves = pos - stackSize;
                int rightMoves = numUnique - pos;

                int leftCost = leftMoves * rotLeftCost;
                int rightCost = rightMoves * rotRightCost;

                int pushpopCost = leftMoves * stackCost * 2 + rotLeftCost + rotRightCost;

                if (pushpopCost < leftCost && pushpopCost < rightCost)
                {
                    sol.cost += pushpopCost;
                    lastPush = 0;
                    if constexpr (genInstrs)
                    {
                        for (int j = 0; j < leftMoves; ++j)
                        {
                            instrs.push_back({PUSH});
                        }
                        instrs.push_back({ROT_LEFT});
                        for (int j = 0; j < leftMoves; ++j)
                        {
                            instrs.push_back({POP});
                        }
                        instrs.push_back({ROT_RIGHT});
                    }
                    std::rotate(stack + pos, stack + pos + 1, stack + numUnique);
                    std::rotate(stack + stackSize, stack + numUnique - 1, stack + numUnique);
                }
                else
                {
                    int cost, moves, op;
                    if (leftCost < rightCost)
                    {
                        cost = leftCost;
                        moves = leftMoves;
                        op = ROT_LEFT;
                    }
                    else
                    {
                        cost = rightCost;
                        moves = rightMoves;
                        op = ROT_RIGHT;
                    }

                    sol.cost += cost;
                    lastPush = 0;
                    if constexpr (genInstrs)
                    {
                        for (int j = 0; j < moves; ++j)
                        {
                            instrs.push_back({op});
                        }
                    }
                    std::rotate(stack + stackSize, stack + pos, stack + numUnique);
                }
            }

            sol.cost += stackCost;
            ++lastPush;
            if constexpr (genInstrs)
            {
                instrs.push_back({PUSH});
            }
            ++stackSize;
        }

        for (int i = curr; i < next; ++i)
        {
            if (lastPush > 0)
            {
                sol.cost -= stackCost;
                --lastPush;
                instrs.pop_back();
            }
            else
            {
                sol.cost += stackCost;
                if constexpr (genInstrs)
                {
                    instrs.push_back({POP});
                }
            }
            --stackSize;
        }

        sol.cost += printCosts[next - curr];
        if constexpr (genInstrs)
        {
            instrs.push_back({PRINT, next - curr});
        }

        std::rotate(stack + stackSize, stack + stackSize + (next - curr), stack + numUnique);
        curr = next;
    }

    return instrs;
}

template<bool genInstrs>
std::vector<Instruction> findInstructions(Solution& sol)
{
    fixSolution(sol);

    if (sol.type == T_MAXPRINTC)
    {
        return findInstructionsMPC<genInstrs>(sol);
    }
    else if (sol.type == T_MAXPRINTR)
    {
        return findInstructionsMPR<genInstrs>(sol);
    }

    int currSeq = 0;
    int currOffsets[NUM_LETS];
    std::fill(currOffsets, currOffsets + NUM_LETS, 0);

    int stackSize = 0;
    int stack[NUM_LETS];
    int stackOffset = 0;
    std::copy(sol.seqs[0].lets, sol.seqs[0].lets + numUnique, stack);

    int lastPrint = 0;

    std::vector<Instruction> instrs;

    sol.cost = optNumSeqCosts[sol.numSeqs];

    for (int i = 0; i < n; ++i)
    {
        int c = text[i];

        if (sol.type == T_CHANGEROT)
        {
            int nextSeq = sol.inSeq[c];
            if (currSeq != nextSeq)
            {
                sol.cost += changeCost;
                lastPrint = 0;
                if constexpr (genInstrs)
                {
                    instrs.push_back({CHANGE, nextSeq});
                }
                currSeq = nextSeq;
            }

            int sz = sol.seqs[currSeq].size;
            int pos = sol.seqs[currSeq].inPos[c];
            int offset = currOffsets[currSeq];
            if (offset != pos)
            {
                int leftMoves = (pos - offset + sz) % sz;
                int rightMoves = sz - leftMoves;

                int leftCost = leftMoves * rotLeftCost;
                int rightCost = rightMoves * rotRightCost;

                int cost, op, moves;

                if (leftCost < rightCost)
                {
                    cost = leftCost;
                    op = ROT_LEFT;
                    moves = leftMoves;
                }
                else
                {
                    cost = rightCost;
                    op = ROT_RIGHT;
                    moves = rightMoves;
                }

                sol.cost += cost;
                lastPrint = 0;
                if constexpr (genInstrs)
                {
                    for (int j = 0; j < moves; ++j)
                    {
                        instrs.push_back({op});
                    }
                }
                currOffsets[currSeq] = pos;
            }
        }
        else if (stack[stackSize] != c)
        {
            if (sol.type == T_ONESTACKS && stackOffset != 0)
            {
                int rightMoves = stackOffset % (numUnique - stackSize);
                int leftMoves = numUnique - stackSize - rightMoves;

                int leftCost = leftMoves * rotLeftCost;
                int rightCost = rightMoves * rotRightCost;

                if (leftCost < rightCost)
                {
                    sol.cost += rotLeftCost;
                    lastPrint = 0;
                    if constexpr (genInstrs)
                    {
                        instrs.push_back({ROT_LEFT});
                    }
                    std::rotate(stack + stackSize, stack + stackSize + 1, stack + numUnique);
                    stackOffset = (stackOffset + 1) % (numUnique - stackSize);
                }
                else
                {
                    sol.cost += rotRightCost;
                    lastPrint = 0;
                    if constexpr (genInstrs)
                    {
                        instrs.push_back({ROT_RIGHT});
                    }
                    std::rotate(stack + stackSize, stack + numUnique - 1, stack + numUnique);
                    stackOffset = (stackOffset - 1 + numUnique - stackSize) % (numUnique - stackSize);
                }

                --i;
                continue;
            }

            int pos = -1;

            for (int j = 0; j < numUnique; ++j)
            {
                if (stack[j] == c)
                {
                    pos = j;
                    break;
                }
            }

            int stackMoveCost = std::abs(stackSize - pos) * stackCost;
        
            int rightMoves = numUnique - pos;
            int rightCost = rightMoves * rotRightCost;

            if (stackSize < pos && rightCost < stackCost)
            {
                sol.cost += rightCost;
                lastPrint = 0;
                if constexpr (genInstrs)
                {
                    for (int j = 0; j < rightMoves; ++j)
                    {
                        instrs.push_back({ROT_RIGHT});
                    }
                }
                std::rotate(stack + stackSize, stack + pos, stack + numUnique);
                stackOffset = (stackOffset - rightMoves + numUnique - stackSize) % (numUnique - stackSize);
            }
            else
            {
                sol.cost += stackMoveCost;
                lastPrint = 0;
                if constexpr (genInstrs)
                {
                    for (int i = pos; i < stackSize; ++i)
                    {
                        instrs.push_back({POP});
                    }
                    for (int i = stackSize; i < pos; ++i)
                    {
                        instrs.push_back({PUSH});
                    }
                }
                stackSize = pos;
            }
        }

        if (lastPrint == 0)
        {
            sol.cost += printCosts[1];
            lastPrint = 1;
            if constexpr (genInstrs)
            {
                instrs.push_back({PRINT, 1});
            }
        }
        else
        {
            sol.cost += printCosts[lastPrint + 1] - printCosts[lastPrint];
            ++lastPrint;
            if constexpr (genInstrs)
            {
                ++instrs.back().arg;
            }
        }

        if (sol.type == T_CHANGEROT)
        {
            int sz = sol.seqs[currSeq].size;
            currOffsets[currSeq] = (currOffsets[currSeq] + 1) % sz;
        }
        else
        {
            for (int i = stackSize; i < numUnique - 1; ++i)
            {
                stack[i] = stack[i + 1];
            }
            stack[numUnique - 1] = c;
            stackOffset = (stackOffset + 1) % (numUnique - stackSize);
        }
    }

    return instrs;
}

void mutateSol(Solution& sol, int depth = 0)
{
    if (depth == 5) return;

    int ssIdx = randNum(0, sol.type != T_CHANGEROT ? 1 : sol.numSeqs);
    int tsIdx = randNum(0, sol.type != T_CHANGEROT ? 1 : std::min(sol.numSeqs + 1, NUM_LETS));

    Sequence& sourceSeq = sol.seqs[ssIdx];
    Sequence& targetSeq = sol.seqs[tsIdx];

    if (sourceSeq.size == 1 && targetSeq.size == 0)
    {
        mutateSol(sol, depth + 1);
        return;
    }

    int sourcePos = randNum(0, sourceSeq.size);
    int c = sourceSeq.lets[sourcePos];

    --sourceSeq.size;
    for (int i = sourcePos; i < sourceSeq.size; ++i)
    {
        sourceSeq.lets[i] = sourceSeq.lets[i + 1];
    }

    int targetPos = randNum(0, targetSeq.size + 1);

    ++targetSeq.size;
    for (int i = targetSeq.size; i > targetPos; --i)
    {
        targetSeq.lets[i] = targetSeq.lets[i - 1];
    }
    targetSeq.lets[targetPos] = c;

    if (ssIdx == tsIdx && sourcePos == targetPos)
    {
        mutateSol(sol, depth + 1);
        return;
    }
}

Solution bestSol;

void output()
{
    std::vector<Instruction> instrs = findInstructions<true>(bestSol);

    outF << optNumSeqs[bestSol.numSeqs] << "\n";
    for (int i = 0; i < bestSol.numSeqs; ++i)
    {
        Sequence& seq = bestSol.seqs[i];
        for (int j = 0; j < seq.size; ++j)
        {
            int c = seq.lets[j];
            if (c == NUM_LETS - 1) outF << "_";
            else outF << (char) (c + 'a');
        }
        outF << "\n";
    }

    for (int i = bestSol.numSeqs; i < optNumSeqs[bestSol.numSeqs]; ++i)
    {
        outF << (char) (nonUsed[i - bestSol.numSeqs] + 'a') << "\n";
    }

    outF << instrs.size() << "\n";
    for (Instruction instr : instrs)
    {
        if (instr.code == ROT_LEFT) outF << "cl";
        else if (instr.code == ROT_RIGHT) outF << "cr";
        else if (instr.code == PRINT) outF << "w " << instr.arg;
        else if (instr.code == CHANGE) outF << "cd " << instr.arg + 1;
        else if (instr.code == PUSH) outF << "push";
        else if (instr.code == POP) outF << "pop";
        else outF << "UNKNOWN";
        outF << "\n";
    }

    std::cerr << "Best cost: " << bestSol.cost << std::endl;

    if constexpr (DEBUG)
    {
        ull cost = numSeqCosts[optNumSeqs[bestSol.numSeqs]];
        for (Instruction instr : instrs)
        {
            if (instr.code == ROT_LEFT) cost += rotLeftCost;
            else if (instr.code == ROT_RIGHT) cost += rotRightCost;
            else if (instr.code == PRINT) cost += printCosts[instr.arg];
            else if (instr.code == CHANGE) cost += changeCost;
            else if (instr.code == PUSH) cost += stackCost;
            else if (instr.code == POP) cost += stackCost;
        }

        if (cost != bestSol.cost)
        {
            std::cerr << "Non-matching cost: " << cost << std::endl;
            while (true) {}
        }
    }
}

Solution sol, sol2;

void genSingletonsSol()
{
    for (int i = 0; i < numUnique; ++i)
    {
        sol.seqs[i].size = 1;
        sol.seqs[i].lets[0] = unique[i];
    }

    for (int i = numUnique; i < NUM_LETS; ++i)
    {
        sol.seqs[i].size = 0;
    }
}

void genOneSeqSol()
{
    sol.seqs[0].size = numUnique;
    for (int i = 0; i < numUnique; ++i)
    {
        sol.seqs[0].lets[i] = unique[i];
    }

    for (int i = 1; i < NUM_LETS; ++i)
    {
        sol.seqs[i].size = 0;
    }
}

void genRandSol()
{
    for (int i = 0; i < NUM_LETS; ++i)
    {
        sol.seqs[i].size = 0;
    }

    if (sol.type == T_CHANGEROT)
    {
        for (int i = 0; i < numUnique; ++i)
        {
            int inSeq = randNum(0, numUnique);
            sol.seqs[inSeq].lets[sol.seqs[inSeq].size++] = unique[i];
        }
    }
    else if (sol.type == T_ONESTACKU || sol.type == T_ONESTACKS || sol.type == T_MAXPRINTR)
    {
        sol.seqs[0].size = numUnique;
        for (int i = 0; i < numUnique; ++i)
        {
            sol.seqs[0].lets[i] = unique[i];
        }
        std::shuffle(sol.seqs[0].lets, sol.seqs[0].lets + numUnique, generator);
    }
    else
    {
        std::cerr << "Cannot randomize type " << TYPE_NAMES[sol.type] << std::endl;
    }
}

void updateSol()
{
    findInstructions<false>(sol);
    if (sol.cost < bestSol.cost)
    {
        bestSol = sol;
        if constexpr (DEBUG)
        {
            std::cerr << "Cost: " << sol.cost << " Type: " << TYPE_NAMES[sol.type] << std::endl;
        }
    }
}

const int RES = 1e6;

void solve()
{
    init();

    findOptMPC();
    findOptMPR();

    int itersDone = 0;

    sol.type = T_ONESTACKS;
    genOneSeqSol();
    updateSol();

    sol.type = T_CHANGEROT;
    genSingletonsSol();
    updateSol();

    genOneSeqSol();
    updateSol();

    while (bestSol.type == T_CHANGEROT && timeLeft(maxTimeLimit / 4) > 0)
    {
        ++itersDone;
    
        genRandSol();
        updateSol();
    }

    std::cerr << "Rand iters: " << itersDone << std::endl;

    itersDone = 0;

    bool middlePassed = false;

    sol = bestSol;
    while (timeLeft(maxTimeLimit * 3 / 4) > 0)
    {
        if (timeLeft(maxTimeLimit * 3 / 8) <= 0 && !middlePassed)
        {
            middlePassed = true;

            std::cerr << "Annealing 1 iters: " << itersDone << std::endl;

            itersDone = 0;

            sol.type = T_ONESTACKU;
            genOneSeqSol();
            updateSol();

            for (int i = 0; i < 10; ++i)
            {
                ++itersDone;
    
                genRandSol();
                updateSol();
            }

            sol.type = T_MAXPRINTR;
            genOneSeqSol();
            updateSol();

            for (int i = 0; i < 10; ++i)
            {
                ++itersDone;
    
                genRandSol();
                updateSol();
            }

            sol = bestSol;
        }

        ++itersDone;
    
        sol2 = sol;
        mutateSol(sol2);
        findInstructions<false>(sol2);
    
        if (sol2.cost < sol.cost)
        {
            sol = sol2;
        }
        else
        {
            double TNorm = timeLeft(maxTimeLimit * 3 / 4) * 2 / maxTimeLimit;
            double T = std::max(TNorm * bestSol.cost / 200, 1.0);
            double P = exp( - (double) (sol2.cost - sol.cost) / T);
            if (randNum(0, RES) < P * RES) sol = sol2;
        }

        updateSol();
    }

    std::cerr << "Annealing 2 iters: " << itersDone << std::endl;

    itersDone = 0;

    while (timeLeft(maxTimeLimit) > 0)
    {
        ++itersDone;
    
        sol = bestSol;
        mutateSol(sol);
        updateSol();
    }

    std::cerr << "Local opt iters: " << itersDone << std::endl;

    sol.type = T_MAXPRINTC;
    genSingletonsSol();
    updateSol();
}

int main()
{
    generator.seed(0);

    startT = std::chrono::high_resolution_clock::now();

    input();
    solve();
    output();

    return 0;
}
