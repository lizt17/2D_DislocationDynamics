#include "MobilityLaw_W.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;
const double mu_SI = 161e9;     // [Pa], shear modulus
const double b_SI = 2.74e-10;   // [m], burgers vector
const double unitSIF = mu_SI * sqrt(b_SI);

int main()
{
    mobilityLaw_W mobilityLaw;
    ofstream outputFile("output_vel.txt");
    
    double Kapp = 1e6 / unitSIF; // Example applied stress in Pa
    double temperature = 200; // Example temperature in Kelvin

    for(int i = 0; i < 1000; i++)
    {
        double r = (i+1) * 2.0;
        double rss = Kapp / sqrt(2*M_PI*r);
        auto result = mobilityLaw.velocity(rss, temperature);
        outputFile << r << " " << result.second << " " << result.first << endl;
    }
    outputFile.close();
    return 0;
}