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
    ofstream outputFile("output_vel_800K.txt");
    
    double Kapp = 1.6e6 / unitSIF; // Example applied stress in Pa
    double temperature = 800; // Example temperature in Kelvin
    double tau_friction = 500e6 / mu_SI; // Friction stress in Pa

    for(int i = 0; i < 1000; i++)
    {
        double r = (i+1) * 2.0;
        double rss = Kapp / sqrt(2*M_PI*r);
        rss = rss <= tau_friction ? 0.0 : rss - tau_friction;
        auto result = mobilityLaw.velocity(rss, temperature);
        auto dGkp = mobilityLaw.dG_kinkpair(rss, temperature);
        outputFile << r << " " << result.second << " " << result.first;
        outputFile << " " << dGkp << endl;
    }
    outputFile.close();
    return 0;
}
