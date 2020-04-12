#include<iostream>
#include<fstream>
#include<vector>
#include<queue>
#include<stack>
#include<algorithm>
#include<math.h>
#include<unordered_map>
#include<conio.h>
#include<chrono>

//#define LOG_FREQS

std::ifstream inF("imrec.in");
std::ofstream outF("imrec.out");

const int MAX_WH = 512;
const int NUM_COLS = 16;
const int LEVELS = 256;
const int NUM_CHANS = 3;
const double SWAP_PROB = 0.2;
const int DIFF_OFFSETS = 3;

/// JEPG STARTS

const int BL_SZ = 8;

const int NUM_CHAN_TYPES = 2;
const int CHAN_TYPES[NUM_CHANS] = {0, 1, 1};

const int Q_BITS = 12;
const double Q_SCALE = 256.0 / (1 << Q_BITS);
const double Q_MULTS[NUM_CHAN_TYPES] = {1, 0.9};

double uncenter(double x)
{
    return x + (LEVELS - 1) * 0.5;
}

int roundClamp(double x)
{
    int a = round(x);
    if (a < 0) return 0;
    if (a >= LEVELS) return LEVELS - 1;
    return a;
}

double center(double a)
{
    return a - (LEVELS - 1) * 0.5;
}

const double RANGE = 4.0 / 3;

struct Color
{
    int r, g, b;
    double xyz[NUM_CHANS];

    void RGB2XYZ()
    {
        const double luma = (r + g + b) / 3.0;
        xyz[0] = center(luma);
        xyz[1] = (r - luma) / RANGE;
        xyz[2] = (g - luma) / RANGE;
    }

    void XYZ2RGB()
    {
        const double luma = uncenter(xyz[0]);
        r = roundClamp(luma + xyz[1] * RANGE);
        g = roundClamp(luma + xyz[2] * RANGE);
        b = roundClamp(luma - (xyz[1] + xyz[2]) * RANGE);
    }
};

#define OKAY 0
#define GOTDC 1
#define NOTHING 2

struct Block
{
    int chanType;
    double vals[BL_SZ][BL_SZ];
    double dct[BL_SZ][BL_SZ];

    void doDCT();
    void genRep(double q, int& lastQDC, std::vector<bool>& data);
    int repSize(double q, int& lastQDC);
    int parseRep(double q, int& lastQDC, const std::vector<bool>& data, int& ptr);
    void doIDCT();
};

/// JEPG ENDS

struct BlockRow
{
    int q;
    int cap;
    std::vector<bool> data;
};

bool run;
int h, w;
int bh, bw;

Color img[MAX_WH][MAX_WH];
int encImg[MAX_WH][MAX_WH];
bool encBits[MAX_WH][MAX_WH * 3];
int pal[NUM_COLS * NUM_CHANS];
BlockRow brows[MAX_WH];
Block blocks[MAX_WH][MAX_WH][NUM_CHANS];

int offset;
int bitBelongs[MAX_WH][MAX_WH * 3];
int bitEnds[MAX_WH], minBitEnds;
int origBPR, origHalfH, origW;
int xOffset, xOffset2, yOffset[DIFF_OFFSETS], browOffset;
int tag, tagX, tagY, wParity;

bool calculated[MAX_WH][MAX_WH][NUM_CHANS];
bool transposed;

void warning(const std::string& s)
{
    std::cerr << "\n WARNING: " << s << "\n" << std::endl;
}

void input()
{
    inF >> run;
    if (!run)
    {
        inF >> w >> h;
        for (int i = 0; i < h; ++i)
        {
            for (int j = 0; j < w; ++j)
            {
                inF >> img[i][j].r >> img[i][j].g >> img[i][j].b;
            }
        }
    }
    else
    {
        for (int i = 0; i < NUM_COLS * NUM_CHANS; ++i)
        {
            inF >> pal[i];
        }
        inF >> w >> h;
        for (int i = 0; i < h; ++i)
        {
            for (int j = 0; j < w; ++j)
            {
                inF >> encImg[i][j];
            }
        }
    }
}

void output()
{
    if (run)
    {
        for (int i = 0; i < h; ++i)
        {
            for (int j = 0; j < w; ++j)
            {
                if (img[i][j].r < 0 || img[i][j].r >= LEVELS)
                {
                    img[i][j].r = LEVELS / 2;
                    warning("Red outside range!");
                }
                if (img[i][j].g < 0 || img[i][j].g >= LEVELS)
                {
                    img[i][j].g = LEVELS / 2;
                    warning("Green outside range!");
                }
                if (img[i][j].b < 0 || img[i][j].b >= LEVELS)
                {
                    img[i][j].b = LEVELS / 2;
                    warning("Blue outside range!");
                }
                outF << img[i][j].r << ' ' << img[i][j].g << ' ' << img[i][j].b << ' ';
            }
            outF << '\n';
        }
    }
    else
    {
        for (int i = 0; i < NUM_COLS * NUM_CHANS; ++i)
        {
            if (pal[i] < 0 || pal[i] >= LEVELS)
            {
                pal[i] = LEVELS / 2;
                warning("Palette color outside range!");
            }
            outF << pal[i] << '\n';
        }
        for (int i = 0; i < h; ++i)
        {
            for (int j = 0; j < w; ++j)
            {
                if (encImg[i][j] < 0 || encImg[i][j] >= NUM_COLS)
                {
                    encImg[i][j] = 0;
                    warning("Encoded pixel outside range!");
                }
                outF << encImg[i][j] << ' ';
            }
            outF << '\n';
        }
    }
}

void error(const std::string& s)
{
    std::cerr << "\n ERROR: " << s << "\n" << std::endl;
    output();
    exit(0);
}

Color img2[MAX_WH][MAX_WH];
int encImg2[MAX_WH][MAX_WH];

void transpose()
{
    for (int i = 0; i < h; ++i)
    {
        for (int j = 0; j < w; ++j)
        {
            img2[i][j] = img[i][j];
            encImg2[i][j] = encImg[i][j];
        }
    }
    for (int i = 0; i < h; ++i)
    {
        for (int j = 0; j < w; ++j)
        {
            img[j][i] = img2[i][j];
            encImg[j][i] = encImg2[i][j];
        }
    }
    std::swap(h, w);
}

void setBHBW(int h, int w)
{
    bh = (h + BL_SZ - 1) / BL_SZ;
    bw = (w + BL_SZ - 1) / BL_SZ;
}

const int TBPP = 4;
const int RPP_OPTS = 2;
const int RPP[RPP_OPTS] = {2, 1};
const int BASE_BPP = TBPP - RPP[0];
const int BPC = 8;

int typeGroup(int x, int y)
{
    const int res = (offset + y - x) % DIFF_OFFSETS;
    return res + (res < 0) * DIFF_OFFSETS;
}

int typeGroupBits(int tpg)
{
    return RPP[tpg & 1];
}

const int NUM_START_BITS = 32;
const int NUM_OGHH_BITS = 8;
const int NUM_TAG_BITS = 2;
const int NUM_W_BITS = 1;
const int NUM_INIT_BITS = NUM_START_BITS + NUM_OGHH_BITS + NUM_TAG_BITS + NUM_W_BITS;
const bool START_BITS[NUM_START_BITS] = {1,0,0,0,0,1,0,1, 1,1,0,1,0,0,0,0, 1,0,1,1,1,1,0,0, 0,1,1,1,0,1,0,1};
const int NUM_SEC_START_BITS = 3;
const bool SEC_START_BITS[NUM_SEC_START_BITS] = {1,1,0};
const int SEC_LB = -2;
const int SEC_UB = 2;

int getBPR(int w)
{
    return (w + 1) / 2 * BASE_BPP + w / 6;
}

void encode()
{
    for (int i = 0; i < bh; ++i)
    {
        while (brows[i].data.size() < brows[i].cap)
        {
            brows[i].data.push_back(rand() % 2);
        }
    }

    int initPos;
    int pos;
    const int halfH = (h + 1) / 2;
    const int halfW = (w + 1) / 2;
    const int BPR = getBPR(w);
    if (BPR < NUM_START_BITS) error("Start bits can't fit!");
    for (int i = 0; i < h; ++i)
    {
        if (i % BL_SZ == 0)
        {
            if (i > 0 && pos < brows[(i - 1) / BL_SZ].data.size())
            {
                warning("Row data didn't fit!");
            }
            pos = 0;
            initPos = -1;
        }
        std::vector<bool> rowBits;
        int tagLoc = -1;
        int k = 0;

        if (i % halfH == 0) initPos = 0;
        else if (i % halfH <= SEC_UB || (i % halfH >= SEC_LB + halfH && i < halfH))
        {
            for (; k < NUM_SEC_START_BITS; ++k)
            {
                rowBits.push_back(SEC_START_BITS[k]);
            }
        }
        if (initPos >= 0 && initPos < NUM_INIT_BITS)
        {
            for (; initPos < NUM_INIT_BITS && k < BPR; ++k, ++initPos)
            {
                if (initPos < NUM_START_BITS)
                {
                    rowBits.push_back(START_BITS[initPos]);
                }
                else if (initPos < NUM_START_BITS + NUM_OGHH_BITS)
                {
                    rowBits.push_back(halfH >> initPos - NUM_START_BITS & 1);
                }
                else if (initPos < NUM_INIT_BITS - NUM_W_BITS)
                {
                    if (tagLoc < 0) tagLoc = rowBits.size();
                    rowBits.push_back(0);
                }
                else
                {
                    rowBits.push_back(w % 2);
                }
            }
        }

        int currBrow = i / BL_SZ;
        for (; k < BPR; ++k)
        {
            if (pos < brows[currBrow].data.size())
            {
                rowBits.push_back(brows[currBrow].data[pos++]);
            }
            else rowBits.push_back(rand() % 2);
        }

        for (int i = 0; i < BPR * 2; ++i)
        {
            rowBits.push_back(rowBits[i % BPR]);
        }
        if (tagLoc >= 0)
        {
            rowBits[tagLoc] = i >= halfH;
            rowBits[BPR + tagLoc] = i >= halfH;
            rowBits[BPR * 2 + tagLoc] = i >= halfH;

            rowBits[tagLoc + 1] = 0;
            rowBits[BPR + tagLoc + 1] = 1;
            rowBits[BPR * 2 + tagLoc + 1] = 0;
        }
        int pos2 = 0;
        for (int j = 0; j < w; ++j)
        {
            const int BPP = TBPP - typeGroupBits(typeGroup(i, j));
            int val = 0;
            for (int k = 0; k < BPP; ++k, ++pos2)
            {
                if (pos2 < rowBits.size()) val = val * 2 + rowBits[pos2];
                else error("Not enough row bits generated!");
            }
            encImg[i][j] = val;
        }
    }
    if (pos < brows[(h - 1) / BL_SZ].data.size())
    {
        warning("Row data didn't fit!");
    }
    for (int i = 0; i < h; ++i)
    {
        for (int j = 0; j < w; ++j)
        {
            const int tpg = typeGroup(i, j);
            encImg[i][j] = encImg[i][j] << typeGroupBits(tpg) | tpg;
        }
    }
}

bool suitable(int x, int y, int x2, int y2)
{
    return typeGroup(x, y) == (encImg[x2][y2] & (1 << typeGroupBits(encImg[x2][y2]) - 1));
}

void fixOneDFS(int x, int y)
{
    if (x < 0 || x >= h || y < 0 || y >= w || suitable(x, y, x, y)) return;
    bool suitableAdj[4];
    suitableAdj[0] = x > 0 && suitable(x, y, x - 1, y) && suitable(x - 1, y, x, y);
    suitableAdj[1] = x < h - 1 && suitable(x, y, x + 1, y) && suitable(x + 1, y, x, y);
    suitableAdj[2] = y > 0 && suitable(x, y, x, y - 1) && suitable(x, y - 1, x, y);
    suitableAdj[3] = y < w - 1 && suitable(x, y, x, y + 1) && suitable(x, y + 1, x, y);
    int totalSuitable = suitableAdj[0] + suitableAdj[1] + suitableAdj[2] + suitableAdj[3];
    if (totalSuitable == 0)
    {
        error("No suitable adjacent cell!");
        return;
    }
    if (totalSuitable > 1)
    {
        return;
    }
    if (suitableAdj[0])
    {
        std::swap(encImg[x - 1][y], encImg[x][y]);
        fixOneDFS(x + 1, y);
        fixOneDFS(x, y - 1);
        fixOneDFS(x, y + 1);
        fixOneDFS(x - 1, y - 1);
        fixOneDFS(x - 1, y + 1);
        fixOneDFS(x - 2, y);
        return;
    }
    if (suitableAdj[1])
    {
        std::swap(encImg[x + 1][y], encImg[x][y]);
        fixOneDFS(x - 1, y);
        fixOneDFS(x, y - 1);
        fixOneDFS(x, y + 1);
        fixOneDFS(x + 1, y - 1);
        fixOneDFS(x + 1, y + 1);
        fixOneDFS(x + 2, y);
        return;
    }
    if (suitableAdj[2])
    {
        std::swap(encImg[x][y - 1], encImg[x][y]);
        fixOneDFS(x, y + 1);
        fixOneDFS(x - 1, y);
        fixOneDFS(x + 1, y);
        fixOneDFS(x - 1, y - 1);
        fixOneDFS(x + 1, y - 1);
        fixOneDFS(x, y - 2);
        return;
    }
    if (suitableAdj[3])
    {
        std::swap(encImg[x][y + 1], encImg[x][y]);
        fixOneDFS(x, y - 1);
        fixOneDFS(x - 1, y);
        fixOneDFS(x + 1, y);
        fixOneDFS(x - 1, y + 1);
        fixOneDFS(x + 1, y + 1);
        fixOneDFS(x, y + 2);
        return;
    }
    warning("Reached end of fixOne!");
}

void findRedOffset()
{
    int offsetCnts[DIFF_OFFSETS];
    for (int i = 0; i < DIFF_OFFSETS; ++i)
    {
        offsetCnts[i] = 0;
    }
    for (int i = 0; i < h; ++i)
    {
        for (int j = 0; j < w; ++j)
        {
            for (offset = 0; offset < DIFF_OFFSETS; ++offset)
            {
                offsetCnts[offset] += suitable(i, j, i, j);
            }
        }
    }
    int maxCnt = 0;
    for (int i = 0; i < DIFF_OFFSETS; ++i)
    {
        if (offsetCnts[i] > maxCnt)
        {
            maxCnt = offsetCnts[i];
            offset = i;
        }
    }
}

void fixSwaps()
{
    findRedOffset();
    for (int i = 0; i < h; ++i)
    {
        for (int j = 0; j < w; ++j)
        {
            fixOneDFS(i, j);
        }
    }
    for (int i = 0; i < h; ++i)
    {
        for (int j = 0; j < w; ++j)
        {
            if (!suitable(i, j, i, j))
            {
                error("Couldn't fix swaps!");
                break;
            }
        }
    }
}

void removeDiffsAndUnrollRows()
{
    for (int i = 0; i < h; ++i)
    {
        int j2 = 0;
        for (int j = 0; j < w; ++j)
        {
            const int tpgBits = typeGroupBits(encImg[i][j]);
            const int BPP = TBPP - tpgBits;
            const int num = encImg[i][j] >> tpgBits;
            for (int k = 0; k < BPP; ++k)
            {
                bitBelongs[i][j2] = j;
                encBits[i][j2++] = (num >> (BPP - 1 - k)) & 1;
            }
        }
        if (i == 0 || j2 < minBitEnds) minBitEnds = j2;
        bitEnds[i] = j2;
    }
}

const int TOTAL_TAG_BITS = NUM_TAG_BITS * 2;

bool equalCols(int y1, int y2, int& diffs)
{
    for (int x = 0; x < h && diffs <= TOTAL_TAG_BITS; ++x)
    {
        if (encBits[x][y1] != encBits[x][y2]) ++diffs;
    }
    return diffs <= TOTAL_TAG_BITS;
}

void findWH()
{
    for (int i = 1; i <= minBitEnds; ++i)
    {
        bool match = true;
        int diffs = 0;
        for (int j = 0; i + j < minBitEnds && match; ++j)
        {
            match = equalCols(j, i + j, diffs);
        }
        if (match)
        {
            origBPR = i;
            break;
        }
    }
}

bool testOffset(int x, int y, const int NUM_START_BITS, const bool START_BITS[])
{
    for (int i = 0; i < NUM_START_BITS; ++i)
    {
        const int y2 = (y + i) % origBPR;
        const bool bit = encBits[x][y2];
        if (bit != START_BITS[i]) return false;
    }
    return true;
}

void findOffsets()
{
    bool done = false;
    for (int i = 0; i < h && !done; ++i)
    {
        for (int j = 0; j < origBPR && !done; ++j)
        {
            if (testOffset(i, j, NUM_START_BITS, START_BITS))
            {
                xOffset = i;
                yOffset[i % DIFF_OFFSETS] = j;
                done = true;
            }
        }
    }
    if (!done) error("Offsets not found!");

    for (int i = SEC_LB; i <= SEC_UB; ++i)
    {
        if (i == 0) continue;
        const int x = xOffset + i;
        if (x < 0 || x >= h) continue;
        done = false;
        for (int j = -1; j <= 1 && !done; ++j)
        {
            const int y = (yOffset[xOffset % DIFF_OFFSETS] + j + origBPR) % origBPR;
            if (testOffset(x, y, NUM_SEC_START_BITS, SEC_START_BITS))
            {
                yOffset[x % DIFF_OFFSETS] = y;
                done = true;
            }
        }
        if (!done) error("Secondary offset indicator not found!");
    }
}

void findInitBits()
{
    tag = 0;
    tagX = -1;
    tagY = -1;
    origHalfH = 0;
    xOffset2 = -1000;
    for (int i = 0; i < NUM_OGHH_BITS + NUM_TAG_BITS + NUM_W_BITS; ++i)
    {
        const int j = i + NUM_START_BITS;
        const int x = (xOffset + j / origBPR);
        const int y = (yOffset[x % DIFF_OFFSETS] + j) % origBPR;
        const bool bit = encBits[x][y];
        if (i < NUM_OGHH_BITS) origHalfH |= bit << i;
        else if (i < NUM_OGHH_BITS + NUM_TAG_BITS)
        {
            tag = tag * 2 + bit;
            if (tagX < 0) tagX = x;
            if (tagY < 0) tagY = y;
        }
        else wParity = bit;
    }
    if (tag / 2 == 0)
    {
        if (xOffset == 0) xOffset2 = origHalfH;
        else warning("Tag doesn't make sense.");
    }
    else xOffset2 = xOffset - origHalfH;
    browOffset = std::min(xOffset, xOffset2) % BL_SZ;
    if (browOffset < 0) browOffset += BL_SZ;
}

void generateData()
{
    for (int i = 0; i < h; ++i)
    {
        for (int j = 0; j < origBPR; ++j)
        {
            const int y = (yOffset[i % DIFF_OFFSETS] + j) % origBPR;
            const int brow = (i + BL_SZ - browOffset) / BL_SZ;
            const bool bit = encBits[i][y];

            if ((i == xOffset || i == xOffset2) && j < NUM_INIT_BITS) continue;
            if (((i >= xOffset + SEC_LB && i <= xOffset + SEC_UB) ||
                 (i >= xOffset2 + SEC_LB && i <= xOffset2 + SEC_UB)) &&
                  j < NUM_SEC_START_BITS) continue;
            brows[brow].data.push_back(bit);
        }
    }
}

void decodeW()
{
    origW = origBPR * 6 / 7;
    switch (origBPR % 7)
    {
    case 0:
        if (wParity) warning("Invalid w parity!");
        break;
    case 2:
    case 4:
        origW += 1 - wParity;
        break;
    case 6:
        if (!wParity) warning("Invalid w parity!");
        break;
    default:
        warning("Invalid BPR caught!");
    }
}

void decode()
{
    removeDiffsAndUnrollRows();
    findWH();
    findOffsets();
    findInitBits();
    generateData();
    decodeW();
}

void alignImage()
{
    int origTagY = -1;
    const int tagBitPos = NUM_START_BITS + NUM_OGHH_BITS;
    const int targetX = (tag / 2 ? origHalfH : 0) + tagBitPos / origBPR;
    int firstTagBit = tagBitPos % origBPR;
    if (tag % 2) firstTagBit += origBPR;
    else if (tagY > firstTagBit) firstTagBit += origBPR * 2;
    int nextBit = 0;
    offset = 0;
    for (int i = 0; ; ++i)
    {
        nextBit += TBPP - typeGroupBits(typeGroup(targetX, i));
        if (nextBit > firstTagBit)
        {
            origTagY = i;
            break;
        }
    }
    const int yOff2 = origTagY - bitBelongs[tagX][tagY];
    if (yOff2 < 0) error("Incorrect offsets!");
    const int xOff2 = BL_SZ - browOffset;
    for (int i = 0; i < h; ++i)
    {
        for (int j = 0; j < w; ++j)
        {
            img[i][j] = img[i + xOff2][j + yOff2];
        }
    }
}

void setCapacities()
{
    for (int i = 0; i < bh; ++i)
    {
        int currH = BL_SZ;
        if (i == bh - 1) currH -= bh * BL_SZ - h;
        brows[i].cap = currH * getBPR(w);
    }
    const int FIRST_HIT = std::min(getBPR(w), NUM_INIT_BITS);
    for (int i = SEC_LB; i <= SEC_UB; ++i)
    {
        int x1 = i;
        int x2 = (h + 1) / 2 + i;
        if (x1 >= 0 && x1 < h)
        {
            if (i != 0) brows[x1 / BL_SZ].cap -= NUM_SEC_START_BITS;
            else
            {
                brows[x1 / BL_SZ].cap -= FIRST_HIT;
                brows[(x1 + 1) / BL_SZ].cap -= NUM_INIT_BITS - FIRST_HIT;
            }
        }
        if (x2 >= 0 && x2 < h)
        {
            if (i != 0) brows[x2 / BL_SZ].cap -= NUM_SEC_START_BITS;
            else
            {
                brows[x2 / BL_SZ].cap -= FIRST_HIT;
                brows[(x2 + 1) / BL_SZ].cap -= NUM_INIT_BITS - FIRST_HIT;
            }
        }
    }
}

/// JEPG STARTS

const int NUM_ACS = 10 * 16 + 2;
const int NUM_DCS = 12;

const int AC_FREQS[NUM_CHAN_TYPES][NUM_ACS] = {
    {40416 , 15592 , 7829 , 4995 , 3351 , 2719 , 1996 , 1633 , 1209 , 1000 , 912 , 658 , 501 , 443 , 270 , 167 , 14994 , 2977 , 962 , 572 , 342 , 264 , 220 , 169 , 114 , 114 , 75 , 69 , 59 , 29 , 16 , 11 , 4762 , 455 , 80 , 44 , 25 , 25 , 18 , 21 , 13 , 11 , 10 , 8 , 15 , 5 , 9 , 1 , 1432 , 74 , 4 , 16 , 1 , 2 , 1 , 19 , 1 , 2 , 4 , 34 , 2 , 2 , 1 , 1 , 455 , 13 , 19 , 20 , 1 , 2 , 4 , 26 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 72 , 2 , 1 , 6 , 1 , 1 , 1 , 6 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 9 , 2 , 4 , 3 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 7 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 12181 , 802},
    {9939 , 3560 , 939 , 962 , 506 , 386 , 286 , 275 , 227 , 119 , 80 , 157 , 95 , 53 , 29 , 12 , 3327 , 604 , 103 , 132 , 57 , 22 , 15 , 46 , 18 , 4 , 8 , 38 , 5 , 1 , 1 , 1 , 1009 , 70 , 31 , 37 , 1 , 2 , 8 , 24 , 6 , 1 , 1 , 8 , 6 , 1 , 1 , 1 , 164 , 2 , 6 , 10 , 1 , 1 , 2 , 9 , 2 , 1 , 1 , 2 , 1 , 1 , 1 , 1 , 14 , 1 , 1 , 8 , 1 , 1 , 1 , 5 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 10 , 1 , 1 , 4 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 5 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 24540 , 99}
};

const int DC_FREQS[NUM_CHAN_TYPES][NUM_DCS] = {
    {1738 , 2286 , 2664 , 2398 , 1650 , 933 , 431 , 121 , 35 , 25 , 1 , 1},
    {10824 , 7018 , 3666 , 1916 , 856 , 203 , 49 , 17 , 2 , 1 , 1 , 1}
};

#ifdef LOG_FREQS

int acFreqs[NUM_CHAN_TYPES][NUM_ACS];
int dcFreqs[NUM_CHAN_TYPES][NUM_DCS];

void loadPrevFreqs()
{
    std::ifstream freqFile("freqs.txt");
    char comma;
    for (int t = 0 ; t < NUM_CHAN_TYPES; ++t)
    {
        for (int i = 0; i < NUM_ACS; ++i)
        {
            if (freqFile)
            {
                if (i) freqFile >> comma;
                freqFile >> acFreqs[t][i];
            }
            else acFreqs[t][i] = 1;
        }
        for (int i = 0; i < NUM_DCS; ++i)
        {
            if (freqFile)
            {
                if (i) freqFile >> comma;
                freqFile >> dcFreqs[t][i];
            }
            else dcFreqs[t][i] = 1;
        }
    }
}

void saveCurrFreqs()
{
    std::ofstream freqFile("freqs.txt");
    for (int t = 0 ; t < NUM_CHAN_TYPES; ++t)
    {
        for (int i = 0; i < NUM_ACS; ++i)
        {
            if (i) freqFile << " , ";
            freqFile << acFreqs[t][i];
        }
        freqFile << '\n';
        for (int i = 0; i < NUM_DCS; ++i)
        {
            if (i) freqFile << " , ";
            freqFile << dcFreqs[t][i];
        }
        freqFile << '\n';
    }
}

#endif // LOG_FREQS

struct HuffmanNode
{
    int code;
    int left, right;
};

int ACHuffRoot[NUM_CHAN_TYPES];
HuffmanNode ACHuffNodes[NUM_CHAN_TYPES][NUM_ACS * 2];
std::vector<bool> ACHuffReps[NUM_CHAN_TYPES][NUM_ACS];

int DCHuffRoot[NUM_CHAN_TYPES];
HuffmanNode DCHuffNodes[NUM_CHAN_TYPES][NUM_DCS * 2];
std::vector<bool> DCHuffReps[NUM_CHAN_TYPES][NUM_DCS];

void findHuffReps(int curr, std::vector<bool>& currRep, const HuffmanNode huffNodes[], std::vector<bool> huffReps[])
{
    if (huffNodes[curr].code >= 0)
    {
        huffReps[curr] = currRep;
        return;
    }
    currRep.push_back(false);
    findHuffReps(huffNodes[curr].left, currRep, huffNodes, huffReps);
    currRep[currRep.size() - 1] = true;
    findHuffReps(huffNodes[curr].right, currRep, huffNodes, huffReps);
    currRep.pop_back();
}

void precomputeHuffman(const int NUM, const int FREQS[], int& huffRoot, HuffmanNode huffNodes[], std::vector<bool> huffReps[])
{
    std::priority_queue<std::pair<int, int>> pq;
    for (huffRoot = 0; huffRoot < NUM; ++huffRoot)
    {
        huffNodes[huffRoot] = {huffRoot, -1, -1};
        pq.push({-FREQS[huffRoot], huffRoot});
    }
    while (pq.size() >= 2)
    {
        std::pair<int, int> left = pq.top();
        pq.pop();
        std::pair<int, int> right = pq.top();
        pq.pop();
        huffNodes[huffRoot] = {-1, left.second, right.second};
        pq.push({left.first + right.first, huffRoot++});
    }
    --huffRoot;
    std::vector<bool> currRep;
    findHuffReps(huffRoot, currRep, huffNodes, huffReps);
}

const double T[BL_SZ][BL_SZ] = {
    { 0.35355339,  0.35355339,  0.35355339,  0.35355339,  0.35355339,  0.35355339,  0.35355339,  0.35355339},
    { 0.49039264,  0.41573481,  0.27778512,  0.09754516, -0.09754516, -0.27778512, -0.41573481, -0.49039264},
    { 0.46193977,  0.19134172, -0.19134172, -0.46193977, -0.46193977, -0.19134172,  0.19134172,  0.46193977},
    { 0.41573481, -0.09754516, -0.49039264, -0.27778512,  0.27778512,  0.49039264,  0.09754516, -0.41573481},
    { 0.35355339, -0.35355339, -0.35355339,  0.35355339,  0.35355339, -0.35355339, -0.35355339,  0.35355339},
    { 0.27778512, -0.49039264,  0.09754516,  0.41573481, -0.41573481, -0.09754516,  0.49039264, -0.27778512},
    { 0.19134172, -0.46193977,  0.46193977, -0.19134172, -0.19134172,  0.46193977, -0.46193977,  0.19134172},
    { 0.09754516, -0.27778512,  0.41573481, -0.49039264,  0.49039264, -0.41573481,  0.27778512, -0.09754516},
};

void Block::doDCT()
{
    for (int i = 0; i < BL_SZ; ++i)
    {
        for (int j = 0; j < BL_SZ; ++j)
        {
            double curr = 0;
            for (int t = 0; t < BL_SZ; ++t)
            {
                for (int k = 0; k < BL_SZ; ++k)
                {
                    curr += T[i][t] * T[j][k] * vals[t][k];
                }
            }
            dct[i][j] = curr;
        }
    }
}

void Block::doIDCT()
{
    for (int i = 0; i < BL_SZ; ++i)
    {
        for (int j = 0; j < BL_SZ; ++j)
        {
            double curr = 0;
            for (int t = 0; t < BL_SZ; ++t)
            {
                for (int k = 0; k < BL_SZ; ++k)
                {
                    curr += T[t][i] * T[k][j] * dct[t][k];
                }
            }
            vals[i][j] = curr;
        }
    }
}

void appendVector(std::vector<bool>& vec1, const std::vector<bool>& vec2)
{
    for (bool bit : vec2)
    {
        vec1.push_back(bit);
    }
}

void repFixedBits(int val, int bits, std::vector<bool>& data)
{
    while (bits--)
    {
        data.push_back(val & 1);
        val >>= 1;
    }
}

int repVarBits(int val, std::vector<bool>& data)
{
    if (val == 0) return 0;
    int bits = 1;
    while (val > 1)
    {
        data.push_back(val & 1);
        val >>= 1;
        ++bits;
    }
    return bits;
}

void repDC(int val, std::vector<bool>& data, int chanType)
{
    const bool sign = val < 0;
    val = abs(val);
    std::vector<bool> valRep;
    const int bits = repVarBits(val, valRep);
    #ifdef LOG_FREQS
        ++dcFreqs[chanType][bits];
    #endif
    if (bits < 0 || bits >= NUM_DCS) error("Invalid DC code.");
    appendVector(data, DCHuffReps[chanType][bits]);
    appendVector(data, valRep);
    if (val) data.push_back(sign);
}

void repAC(int val, int zeros, std::vector<bool>& data, int chanType)
{
    if (val)
    {
        const bool sign = val < 0;
        val = abs(val);
        std::vector<bool> valRep;
        const int bits = repVarBits(val, valRep);
        const int code = bits - 1 << 4 | zeros;
        if (zeros < 0 || zeros >= 16 || code < 0 || code >= NUM_ACS - 2) error("Invalid AC code.");
        #ifdef LOG_FREQS
            ++acFreqs[chanType][code];
        #endif
        appendVector(data, ACHuffReps[chanType][code]);
        appendVector(data, valRep);
        data.push_back(sign);
    }
    else if (zeros == 0)
    {
        #ifdef LOG_FREQS
            ++acFreqs[chanType][NUM_ACS - 2];
        #endif
        appendVector(data, ACHuffReps[chanType][NUM_ACS - 2]);
    }
    else if (zeros == 15)
    {
        #ifdef LOG_FREQS
            ++acFreqs[chanType][NUM_ACS - 1];
        #endif
        appendVector(data, ACHuffReps[chanType][NUM_ACS - 1]);
    }
    else error("Invalid AC code.");
}

const int ZigZagInv[BL_SZ * BL_SZ][2] = {
    {0, 0}, {0, 1}, {1, 0}, {2, 0}, {1, 1}, {0, 2}, {0, 3}, {1, 2},
    {2, 1}, {3, 0}, {4, 0}, {3, 1}, {2, 2}, {1, 3}, {0, 4}, {0, 5},
    {1, 4}, {2, 3}, {3, 2}, {4, 1}, {5, 0}, {6, 0}, {5, 1}, {4, 2},
    {3, 3}, {2, 4}, {1, 5}, {0, 6}, {0, 7}, {1, 6}, {2, 5}, {3, 4},
    {4, 3}, {5, 2}, {6, 1}, {7, 0}, {7, 1}, {6, 2}, {5, 3}, {4, 4},
    {3, 5}, {2, 6}, {1, 7}, {2, 7}, {3, 6}, {4, 5}, {5, 4}, {6, 3},
    {7, 2}, {7, 3}, {6, 4}, {5, 5}, {4, 6}, {3, 7}, {4, 7}, {5, 6},
    {6, 5}, {7, 4}, {7, 5}, {6, 6}, {5, 7}, {6, 7}, {7, 6}, {7, 7}
};

void Block::genRep(double q, int& lastQDC, std::vector<bool>& data)
{
    double qInv = 1.0 / q;
    const int QDC = round(dct[0][0] * qInv);
    repDC(QDC - lastQDC, data, chanType);
    lastQDC = QDC;
    int zeros = 0;
    for (int k = 1; k < BL_SZ * BL_SZ; ++k)
    {
        const int curr = round(dct[ZigZagInv[k][0]][ZigZagInv[k][1]] * qInv);
        if (curr)
        {
            while (zeros > 15)
            {
                zeros -= 16;
                repAC(0, 15, data, chanType);
            }
            repAC(curr, zeros, data, chanType);
            zeros = 0;
        }
        else ++zeros;
    }
    if (zeros) repAC(0, 0, data, chanType);
}

int countBits(int val)
{
    int bits = 0;
    while (val)
    {
        val >>= 1;
        ++bits;
    }
    return bits;
}

const int DC_ARGS = 1 << 12;
int DCSizeCache[NUM_CHAN_TYPES][DC_ARGS];
int repDCSize(int val, int chanType)
{
    const int arg = abs(val) << 1 | (val < 0);
    if (DCSizeCache[chanType][arg]) return DCSizeCache[chanType][arg];
    const int bits = countBits(abs(val));
    return DCSizeCache[chanType][arg] = DCHuffReps[chanType][bits].size() + bits;
}

const int AC_ARGS = 1 << 15;
int ACSizeCache[NUM_CHAN_TYPES][AC_ARGS];
int repACSize(int val, int zeros, int chanType)
{
    if (val)
    {
        const int arg = abs(val) << 5 | (val < 0) << 4 | zeros;
        if (ACSizeCache[chanType][arg]) return ACSizeCache[chanType][arg];
        val = abs(val);
        const int bits = countBits(val);
        return ACSizeCache[chanType][arg] = ACHuffReps[chanType][bits - 1 << 4 | zeros].size() + bits;
    }
    else if (zeros == 0) return ACHuffReps[chanType][NUM_ACS - 2].size();
    else return ACHuffReps[chanType][NUM_ACS - 1].size();
}

int Block::repSize(double q, int& lastQDC)
{
    const double qInv = 1.0 / q;
    const int QDC = round(dct[0][0] * qInv);
    int totalSize = repDCSize(QDC - lastQDC, chanType);
    lastQDC = QDC;
    int zeros = 0;
    for (int k = 1; k < BL_SZ * BL_SZ; ++k)
    {
        const int curr = round(dct[ZigZagInv[k][0]][ZigZagInv[k][1]] * qInv);
        if (curr)
        {
            while (zeros > 15)
            {
                zeros -= 16;
                totalSize += repACSize(0, 15, chanType);
            }
            totalSize += repACSize(curr, zeros, chanType);
            zeros = 0;
        }
        else ++zeros;
    }
    if (zeros) totalSize += repACSize(0, 0, chanType);
    return totalSize;
}

bool parseBit(const std::vector<bool>& data, int& ptr)
{
    if (ptr == data.size()) throw 0;
    return data[ptr++];
}

int parseHuffmanCode(int curr, const HuffmanNode huffNodes[], const std::vector<bool>& data, int& ptr)
{
    const HuffmanNode& currNode = huffNodes[curr];
    if (currNode.code >= 0) return currNode.code;
    const bool bit = parseBit(data, ptr);
    if (bit) return parseHuffmanCode(currNode.right, huffNodes, data, ptr);
    else return parseHuffmanCode(currNode.left, huffNodes, data, ptr);
}

int parseBits(int bits, const std::vector<bool>& data, int& ptr, bool skipMSOne)
{
    if (bits == 0) return 0;
    bits -= skipMSOne;
    if (ptr + bits - 1 == data.size()) throw 0;
    int val = skipMSOne ? 1 << bits : 0;
    for (int i = 0; i < bits; ++i)
    {
        val |= parseBit(data, ptr) << i;
    }
    return val;
}

int parseDC(const std::vector<bool>& data, int& ptr, int chanType)
{
    const int bits = parseHuffmanCode(DCHuffRoot[chanType], DCHuffNodes[chanType], data, ptr);
    const int val = parseBits(bits, data, ptr, true);
    const bool sign = val && parseBit(data, ptr);
    if (sign) return -val;
    else return val;
}

std::pair<int, int> parseAC(const std::vector<bool>& data, int& ptr, int chanType)
{
    const int code = parseHuffmanCode(ACHuffRoot[chanType], ACHuffNodes[chanType], data, ptr);
    const int zeros = code < NUM_ACS - 2 ? code & 0b1111 : (code + 2 - NUM_ACS) * 15;
    const int val = code < NUM_ACS - 2 ? parseBits((code >> 4) + 1, data, ptr, true) : 0;
    const bool sign = code < NUM_ACS - 2 && parseBit(data, ptr);
    if (sign) return {-val, zeros};
    else return {val, zeros};
}

int Block::parseRep(double q, int& lastQDC, const std::vector<bool>& data, int& ptr)
{
    try
    {
        lastQDC += parseDC(data, ptr, chanType);
    }
    catch (int)
    {
        return NOTHING;
    }
    dct[0][0] = lastQDC * q;
    int pos = 1;
    while (pos < BL_SZ * BL_SZ)
    {
        std::pair<int, int> valueZeros;
        try
        {
            valueZeros = parseAC(data, ptr, chanType);
        }
        catch (int)
        {
            return GOTDC;
        }
        if (valueZeros.second == 0 && valueZeros.first == 0) break;
        pos += valueZeros.second;
        dct[ZigZagInv[pos][0]][ZigZagInv[pos][1]] = valueZeros.first * q;
        ++pos;
    }
    return OKAY;
}

void extendImage()
{
    int lastX = bh * BL_SZ;
    int lastY = bw * BL_SZ;
    for (int i = 0; i < h; ++i)
    {
        for (int j = w; j < lastY; ++j)
        {
            img[i][j] = img[i][w - 1];
        }
    }
    for (int j = 0; j < lastY; ++j)
    {
        for (int i = h; i < lastX; ++i)
        {
            img[i][j] = img[h - 1][j];
        }
    }
}

void setBlockChans()
{
    for (int i = 0; i < bh; ++i)
    {
        for (int j = 0; j < bw; ++j)
        {
            for (int k = 0; k < NUM_CHANS; ++k)
            {
                blocks[i][j][k].chanType = CHAN_TYPES[k];
            }
        }
    }
}

void RGB2XYZ()
{
    for (int i = 0; i < bh * BL_SZ; ++i)
    {
        for (int j = 0; j < bw * BL_SZ; ++j)
        {
            img[i][j].RGB2XYZ();
        }
    }
}

void XYZ2RGB()
{
    for (int i = 0; i < bh * BL_SZ; ++i)
    {
        for (int j = 0; j < bw * BL_SZ; ++j)
        {
            img[i][j].XYZ2RGB();
        }
    }
}

void setBlockVals()
{
    for (int i = 0; i < bh; ++i)
    {
        for (int j = 0; j < bw; ++j)
        {
            for (int i2 = 0; i2 < BL_SZ; ++i2)
            {
                for (int j2 = 0; j2 < BL_SZ; ++j2)
                {
                    Color& curr = img[i * BL_SZ + i2][j * BL_SZ + j2];
                    for (int k = 0; k < NUM_CHANS; ++k)
                    {
                        blocks[i][j][k].vals[i2][j2] = curr.xyz[k];
                    }
                }
            }
        }
    }
}

void getBlockVals()
{
    for (int i = 0; i < bh; ++i)
    {
        for (int j = 0; j < bw; ++j)
        {
            for (int i2 = 0; i2 < BL_SZ; ++i2)
            {
                for (int j2 = 0; j2 < BL_SZ; ++j2)
                {
                    Color& curr = img[i * BL_SZ + i2][j * BL_SZ + j2];
                    for (int k = 0; k < NUM_CHANS; ++k)
                    {
                        curr.xyz[k] = blocks[i][j][k].vals[i2][j2];
                    }
                }
            }
        }
    }
}

void doDCT()
{
    for (int i = 0; i < bh; ++i)
    {
        for (int j = 0; j < bw; ++j)
        {
            for (int k = 0; k < NUM_CHANS; ++k)
            {
                blocks[i][j][k].doDCT();
            }
        }
    }
}

void doIDCT()
{
    for (int i = 0; i < bh; ++i)
    {
        for (int j = 0; j < bw; ++j)
        {
            for (int k = 0; k < NUM_CHANS; ++k)
            {
                blocks[i][j][k].doIDCT();
            }
        }
    }
}

std::pair<int, int> getJK(int pos)
{
    return {(pos + bw / 3) % bw, pos / bw};
}

int repSize(int brow, double q)
{
    int totalSize = Q_BITS;
    int lastQDCs[NUM_CHANS] = {0, 0, 0};
    q *= Q_SCALE;
    for (int pos = 0; pos < bw * NUM_CHANS; ++pos)
    {
        const std::pair<int, int> jk = getJK(pos);
        Block& curr = blocks[brow][jk.first][jk.second];
        totalSize += curr.repSize(q * Q_MULTS[curr.chanType], lastQDCs[jk.second]);
    }
    return totalSize;
}

void genRep(int brow)
{
    int lastQDCs[NUM_CHANS] = {0, 0, 0};
    double q = brows[brow].q * Q_SCALE;
    repFixedBits(brows[brow].q, Q_BITS, brows[brow].data);
    for (int pos = 0; pos < bw * NUM_CHANS; ++pos)
    {
        const std::pair<int, int> jk = getJK(pos);
        Block& curr = blocks[brow][jk.first][jk.second];
        curr.genRep(q * Q_MULTS[curr.chanType], lastQDCs[jk.second], brows[brow].data);
    }
}

void parseRep(int brow)
{
    int ptr = 0;
    int lastQDCs[NUM_CHANS] = {0, 0, 0};
    try
    {
        brows[brow].q = parseBits(Q_BITS, brows[brow].data, ptr, false);
    }
    catch (int)
    {
        return;
    }
    double q = brows[brow].q * Q_SCALE;
    int status = OKAY;
    for (int pos = 0; pos < bw * NUM_CHANS; ++pos)
    {
        const std::pair<int, int> jk = getJK(pos);
        Block& curr = blocks[brow][jk.first][jk.second];
        status = curr.parseRep(q * Q_MULTS[curr.chanType], lastQDCs[jk.second], brows[brow].data, ptr);
        if (status != NOTHING) calculated[brow][jk.first][jk.second] = true;
        if (status != OKAY) return;
    }
}

void findQ(int brow)
{
    int left = ceil(1 / Q_SCALE) - 1;
    int right = 255 / Q_SCALE;
    while (right - left > 1)
    {
        const int mid = (left + right + 1) / 2;
        if (repSize(brow, mid) <= brows[brow].cap) right = mid;
        else left = mid;
    }
    brows[brow].q = right;
}

void fillMissing()
{
    for (int i = 0; i < bh; ++i)
    {
        for (int j = 0; j < bw; ++j)
        {
            for (int k = 0; k < NUM_CHANS; ++k)
            {
                if (calculated[i][j][k]) continue;
                const int iOther = i ? i * BL_SZ - 1 : BL_SZ;
                for (int j2 = 0; j2 < BL_SZ; ++j2)
                {
                    const int y = j * BL_SZ + j2;
                    const double val = img[iOther][y].xyz[k];
                    for (int i2 = 0; i2 < BL_SZ; ++i2)
                    {
                        const int x = i * BL_SZ + i2;
                        img[x][y].xyz[k] = val;
                    }
                }
            }
        }
    }
}

int diff(const Color& a, const Color& b)
{
    int d = (a.r - b.r) * (a.r - b.r);
    d += (a.g - b.g) * (a.g - b.g);
    d += (a.b - b.b) * (a.b - b.b);
    return d;
}

const double BLUR = 60 / 1e3;
const double MAX_WEIGHT = 0.25;

void postEffect()
{
    for (int i = 0; i < bh * BL_SZ; ++i)
    {
        for (int j = 0; j < bw * BL_SZ; ++j)
        {
            int cnt = 0;
            double totalDiff = 0;
            for (int i2 = -1; i2 <= 1; ++i2)
            {
                for (int j2 = -1; j2 <= 1; ++j2)
                {
                    if (i2 == 0 && j2 == 0) continue;
                    if (i + i2 < 0 || i + i2 >= bh * BL_SZ || j + j2 < 0 || j + j2 > bw * BL_SZ) continue;
                    const Color& other = img[i + i2][j + j2];
                    totalDiff += diff(img[i][j], other);
                    ++cnt;
                }
            }
            const double noiseAppr = brows[i / BL_SZ].q * Q_SCALE * std::sqrt(totalDiff / cnt);
            double r, g, b, totalW = 1;
            r = img[i][j].r;
            g = img[i][j].g;
            b = img[i][j].b;
            for (int i2 = -1; i2 <= 1; ++i2)
            {
                for (int j2 = -1; j2 <= 1; ++j2)
                {
                    if (i2 == 0 && j2 == 0) continue;
                    if (i + i2 < 0 || i + i2 >= bh * BL_SZ || j + j2 < 0 || j + j2 > bw * BL_SZ) continue;
                    const Color& other = img[i + i2][j + j2];
                    const double colDist = 1 + diff(img[i][j], other);
                    const double physDist = std::sqrt(i2 * i2 + j2 * j2);
                    const double w = std::min(BLUR * noiseAppr / (physDist * colDist), MAX_WEIGHT);
                    r += other.r * w;
                    g += other.g * w;
                    b += other.b * w;
                    totalW += w;
                }
            }
            img2[i][j].r = roundClamp(r / totalW);
            img2[i][j].g = roundClamp(g / totalW);
            img2[i][j].b = roundClamp(b / totalW);
        }
    }
    for (int i = 0; i < bh * BL_SZ; ++i)
    {
        for (int j = 0; j < bw * BL_SZ; ++j)
        {
            img[i][j] = img2[i][j];
        }
    }
}

void compress()
{
    extendImage();
    RGB2XYZ();
    setBlockChans();
    setBlockVals();
    doDCT();
    for (int i = 0; i < bh; ++i)
    {
        findQ(i);
        genRep(i);
    }
    for (int i = 1; i < NUM_COLS * NUM_CHANS; i += NUM_CHANS)
    {
        pal[i] = rand() % (LEVELS - 30);
    }
}

void decompress()
{
    int ptr = 0;
    setBlockChans();
    for (int i = 1; i < bh; ++i)
    {
        parseRep(i);
    }
    doIDCT();
    getBlockVals();
    fillMissing();
    XYZ2RGB();
    postEffect();
}

/// JPEG ENDS

bool shouldTranspose()
{
    if (getBPR(w) < NUM_INIT_BITS && getBPR(h) < NUM_INIT_BITS) return h > w;
    if (getBPR(w) < NUM_INIT_BITS) return true;
    if (getBPR(h) < NUM_INIT_BITS) return false;
    return w > h;
}

void solveFirst()
{
    #ifdef LOG_FREQS
        loadPrevFreqs();
    #endif

    transposed = shouldTranspose();
    if (transposed) transpose();
    pal[0] = transposed;
    setBHBW(h, w);
    setCapacities();
    compress();
    encode();
    if (transposed) transpose();

    #ifdef LOG_FREQS
        saveCurrFreqs();
    #endif
}

void solveSecond()
{
    transposed = pal[0];
    if (transposed) transpose();
    fixSwaps();
    decode();
    bh = (h - browOffset) / BL_SZ + 2;
    bw = (origW + BL_SZ - 1) / BL_SZ;
    decompress();
    alignImage();
    if (transposed) transpose();
}

void solve()
{
    precomputeHuffman(NUM_ACS, AC_FREQS[0], ACHuffRoot[0], ACHuffNodes[0], ACHuffReps[0]);
    precomputeHuffman(NUM_DCS, DC_FREQS[0], DCHuffRoot[0], DCHuffNodes[0], DCHuffReps[0]);
    precomputeHuffman(NUM_ACS, AC_FREQS[1], ACHuffRoot[1], ACHuffNodes[1], ACHuffReps[1]);
    precomputeHuffman(NUM_DCS, DC_FREQS[1], DCHuffRoot[1], DCHuffNodes[1], DCHuffReps[1]);

    if (!run) solveFirst();
    else solveSecond();
}

std::chrono::high_resolution_clock::time_point startT, currT;

int main()
{
    startT = std::chrono::high_resolution_clock::now();
    currT = std::chrono::high_resolution_clock::now();
    srand(420);

    input();
    solve();
    output();

    currT = std::chrono::high_resolution_clock::now();

    std::cerr << "\n TOTAL TIME: " << std::chrono::duration_cast<std::chrono::duration<double>>(currT - startT).count() << "\n" << std::endl;

    return 0;
}
