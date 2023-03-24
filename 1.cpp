#include <iostream>
#include <math.h>
#include <random>
#include <limits>
#include <algorithm>

double rastriginFunction(const double x){
	return 10 + x * x - 10 * cos(2 * M_PI * x);
}

double randomPoint(const double left, const double right){
	double len = right - left;
	return left + len * (std::rand() / (double)(RAND_MAX));
}


double randomSearch(double (*func)(double), 
                    const double left, const double right,
                    const double minDV, const int maxSteps)
{
	double minV = func((left + right) / 2 + randomPoint(-1, 1));
	int step = 0;
	do
    {  
		double point = randomPoint(left, right);
		double value = func(point);
		if (value < minV) minV = value;		
		++step;
	} while (abs(minV) > minDV && step <  maxSteps);
    std::cout << "Iteraions = " << step << " answer = ";
	return minV;
}



double lipo(double (*func)(double), 
            const double left, const double right,
            double minDV, double minDx, int maxSteps)
{
	double xMax = randomPoint(left, right);
	double yMax = -func(xMax);
    std::vector<float> points;
    points.push_back(xMax);
    double L = 100;

    int finishI;

    for (int i = 0; i < maxSteps; ++i)
	{
        double x = randomPoint(left, right);
        double minV = -func(points.front()) + L * abs(x - points.front());
        for (size_t j = 1; j < points.size(); j++) minV = std::min(minV, -func(points[j]) + L * abs(x - points[j]));

        if (minV >= yMax)
        {
            points.push_back(x);
            double y = -1 * func(x);
			if (y > yMax)
			{
				if (abs(xMax - x) < minDx || abs(yMax - y) < minDV) {
					xMax = x;
                    finishI = i;
					break;
				}
				xMax = x;
				yMax = y;
			}   
        }
    }
    std::cout << "iterations = " << finishI << " answer = ";
    return xMax;
}

double negative_Rastricin(double x){
	return -1 *  rastriginFunction(x);
}

int main()
{
	std::cout << "Random search: " << randomSearch(rastriginFunction, -10, 10, 1e-4, 1e4) << std::endl;
	std::cout << "Lipo : " << lipo(rastriginFunction, -10, 10, 1e-4, 1e-2, 1e4) << std::endl;;
	return 0;
}