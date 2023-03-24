#include <iostream>
#include <bitset>

#include <math.h>
#include <random>
#include <limits>
#include <exception>
#include <map>
#include <bit>
#include <ctime>

double rastriginFunction(const double x){
	return 10 + x * x - 10 * cos(2 * M_PI * x);
}

uint64_t number_to_gray(uint64_t n) {
    return n ^ (n >> 1);
}

uint64_t gray_to_number(uint64_t g) {
    uint64_t n = 0;
    for (; g!=0; g >>= 1) { n ^= g; }
    return n;
}

double randomPoint(const double left, const double right){
	double len = right - left;
	return left + len * (std::rand() / (double)(RAND_MAX));
}

void mutate(std::vector<uint64_t> & g, int k) //k - сколько битов меняем
{
    for(int i = 0; i < g.size(); ++i)
    {
        for(int j = 0; j < k; ++j)
        {
            size_t r = rand() % 64;
            g[i] ^= (1 << r);
        }
    }
}

uint64_t cross(uint64_t left, uint64_t right)
{
    int k = 1 + rand() % 63;
    std::bitset<64> leftBits = static_cast<std::bitset<64>>(left);
    std::bitset<64> rightBits = static_cast<std::bitset<64>>(right);
    std::bitset<64> answBits;
    for(int i = 0; i < 64; ++i)
    {
        if(i < k) answBits[i] = leftBits[i];
        else answBits[i] = rightBits[i];
    }
    return answBits.to_ullong();
}

void makeCrossG(std::vector<uint64_t> & g, int len)
{
    int half = g.size() / 2;
    std::vector<uint64_t> g1, g2;
    for(int i = 0; i < g.size(); ++i)
    {
        if(i < half) g1.push_back(g[i]);
        else g2.push_back(g[i]);
    }
    std::vector<uint64_t> ans;
    for(int i = 0; i < len; ++i)
    {
        int a = rand() % g1.size();
        int b = rand() % g2.size();
        uint64_t xCross = cross(g1[a], g2[b]);
        ans.push_back(xCross);
    }
    std::swap(g, ans);
}

std::vector<double> crossAndMutate(std::vector<double> xs)
{
    std::vector<uint64_t> g(xs.size());
    for(int i = 0; i < xs.size(); ++i)
    {
        g[i] = number_to_gray(std::bit_cast<uint64_t>(xs[i]));
    }
    mutate(g, 5);
    makeCrossG(g, xs.size());

    std::vector<double> h(xs.size());
    for(int i = 0; i < xs.size(); ++i)
    {
        h[i] = static_cast<double>(gray_to_number(g[i]));
    }
    return h;
    
}

double genetic(double (*func)(double), 
                const double left, const double right,
                const int epochs,
                const int randomSampleSize, const int selectionSampleSize)
{
    if(randomSampleSize < selectionSampleSize) throw std::exception();

    std::map<double, double> points; // первая коорлината y, вторая x

    for(int i = 0; i < randomSampleSize; ++i) {
        double x = randomPoint(left, right);
        points[func(x)] = x;
    }

    std::map<double, double> bestPoints;
    std::map<double, double>::const_iterator it = points.begin();
    for(int i = 0; i < selectionSampleSize; ++i)
    {
        bestPoints[it->first] = it->second;
        ++it;
    }

    double min = bestPoints.begin()->second;

    for(int i = 0; i < epochs; ++i)
    {
        std::vector<double> xs;
        for(std::map<double, double>::const_iterator itBest = bestPoints.begin(); itBest != bestPoints.end(); ++itBest) xs.push_back(itBest->second);
        xs = crossAndMutate(xs);

        std::map<double, double>::const_iterator it = points.begin();
        for(int i = 0; i < selectionSampleSize; ++i)
        {
            bestPoints[it->first] = it->second;
            ++it;
        }
        min = bestPoints.begin()->second;
    }
    return min;

}
int main()
{
    srand ( time(NULL) );
    std::cout << genetic(rastriginFunction, -10, 10, 1000, 200, 100) << std::endl;
    return 0;
}