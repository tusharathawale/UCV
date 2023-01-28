
#include <iostream>
#include <vector>
#include <bitset>

void traverse3d(int depth, int id, const int numPoints,
                std::vector<int> &path)
{

    if (depth == numPoints)
    {
        // print path
        std::string bitsetstring;
        for (int i = 0; i < path.size(); i++)
        {
            std::cout << path[i] << ",";
            bitsetstring = bitsetstring + std::to_string(path[i]);
        }
        std::cout << std::endl;
        // transfer path to the bitset
        std::bitset<8> bitset(bitsetstring);

        std::cout << "bitset " << bitset.to_ulong() << std::endl;
        return;
    }
    // two branches for current node

    path.push_back(0);

    // put left
    traverse3d(depth + 1, 1 + (id << 1), numPoints, path);

    path.pop_back();

    path.push_back(1);

    // put right
    traverse3d(depth + 1, id << 1, numPoints, path);

    path.pop_back();

    return;
}

void intToBisStr()
{
    for (uint i = 0; i < 256; i++)
    {
        std::bitset<8> b(i);
        // get each option and transfer to int
        for (uint j = 0; j < 8; j++)
        {
            if (b.test(j))
            {
                std::cout << 1 << ",";
            }
            else
            {
                std::cout << 0 << ",";
            }
        }
        std::cout << std::endl;
    }
}

void intToBisStr2()
{
    for (uint i = 0; i < 256; i++)
    {

        // get each option and transfer to int
        for (uint j = 0; j < 8; j++)
        {
            if (i & (1 << j))
                std::cout << 1 << ",";
            else
                std::cout << 0 << ",";
        }
        std::cout << std::endl;
    }
}

void BitToInt()
{
    // construct the case table
    // test all 0
    // output case is 0
    // default is all zero
    uint caseVal = 0;

    // test all 1
    // output case is 1111111 in binary
    // which is 255 in int
    for (uint i = 0; i < 8; i++)
    {
        // set ith digit as 1
        caseVal = (1 << i) | caseVal;
    }

    std::cout << "caseVal all 1 " << caseVal << std::endl;

    // half 1
    caseVal = 0;
    for (uint i = 0; i < 8; i++)
    {
        if (i <= 3)
        {
            // set ith digit as 1
            caseVal = (1 << i) | caseVal;
        }
    }
    std::cout << "caseVal half 1 " << caseVal << std::endl;
}

int main()
{
    // std::vector<int> path;
    // traverse3d(0, 0, 8, path);
    // intToBisStr();
    intToBisStr2();
    BitToInt();
}