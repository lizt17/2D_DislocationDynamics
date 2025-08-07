// dislocation_dynamics.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <map>
#include <memory>
#include <iomanip>
#include <filesystem>
#include <cstdlib>
#include "MobilityLaw_W.hpp"
#include "DislocationProperty.hpp"

/*------------------------------------------------------------------------------*/
// Global constants loaded from "matPara_W.mat"
const double mu_SI = 161e9;     // [Pa], shear modulus
const double b_SI = 2.74e-10;   // [m], burgers vector
const double rho_SI = 19250.0;  // [kg/m^3]
const double cs_SI = sqrt(mu_SI/rho_SI); // [m/s], speed of shear wave
const double tauP_SI = 2030e6;  // [Pa]
const double KG_SI = 2.1e6; // [Pa], clevage fracture toughness
const double Ke_SI = 1.8e6;
const double unitSIF = mu_SI * sqrt(b_SI);
const double unitSIFrate = mu_SI * cs_SI / sqrt(b_SI);
const double unitTime = b_SI / cs_SI;
/*------------------------------------------------------------------------------*/
// 
const double mu = 1.0;
const double b = 1.0;

// Parameters
const double tau_nuc = 500e6 / mu_SI;
const double tau_friction_SI = 1000.0e6; // [Pa], friction stress
const double tau_friction = tau_friction_SI / mu_SI; 
const double r_source = 200.0;
const double T = 400;
const double crack_tip = 0.0;
const double KappDot = 0e6 / unitSIFrate;
const double Kapp0 = 1.6e6 / unitSIF;
// const double Kapp0 = 1e6 / unitSIF;
const long int Nsteps = 2e7;
const int outputNum = 200;
const int outputInterval = Nsteps / outputNum;
const double shearWaveFraction = 1e-5;
const int maxNumDis = 100;
const double Ke = Ke_SI / unitSIF;
const double rs = std::pow(mu,2)*std::pow(b,2) / (8*M_PI*std::pow(Ke,2));
// const double r_source = rs;

double dxMax = 10.0;

// Stress functions
double tau_interaction(double ri, double rj) 
{
    return mu * b / (2 * M_PI) / (ri - rj) * sqrt(rj/ri);
}
double tau_image(double r) 
{
    return mu * b / (4 * M_PI) / (r - crack_tip);
}
double tau_applied(double Kapp, double r) 
{
    return Kapp / sqrt(2 * M_PI * r);
}

int main() {
    std::cout << "Dislocation dynamics" << std::endl;
    std::cout << "SIF [mu*b] = " << unitSIF << std::endl;
    std::cout << "SIF rate [mu*b/cs] = " << unitSIFrate << std::endl;
    std::cout << "Time unit [b/cs] = " << unitTime << std::endl;
    mobilityLaw_W mobilityLaw;
    
    int Nd = 0, nextDisID = 0;
    std::vector<std::shared_ptr<Dislocation>> disArr;
    double K0nuc = tau_nuc * sqrt(2.0 * M_PI * r_source) + mu * b / 2.0 / sqrt(2 * M_PI * r_source);
    std::cout << "nucleation Kapp0 [MPa m^0.5]= " << K0nuc * unitSIF / 1e6<< std::endl;
    std::cout << "distance to source [b] = " << r_source << std::endl;

    std::string outputDir = "outputVars/";
    std::ofstream outfile(outputDir+"outputVars_"+std::to_string(int(T))+"K.csv", std::ios::trunc);

    int outFileIndex = 0;
    // outfile << "time" << std::endl << "Number of dislocations" << std::endl << "RSS on source" << std::endl
    // << "Kapp" << std::endl << "Ktip";
    // std::ofstream fileRSSsource(outputDir+"rss_source.csv");
    // std::ofstream fileKtip(outputDir+"Ktip.csv");
    // std::ofstream fileBackStress(outputDir+"back_stress.csv");

    double time = 0.0;
    double Kapp;
    for (int kInc = 0; kInc < Nsteps; ++kInc) 
    {
        Kapp = Kapp0 + KappDot * time;
        // std::cout << "Step " << kInc << ", time = " << time << std::endl;
        // Back stress
        double back_stress = 0.0;
        double KD = 0.0;
        for (int nd = 0; nd < Nd; ++nd) 
        {
            auto P = disArr[nd]->getPosition();
            back_stress += tau_interaction(r_source, P);
            KD += mu * b / sqrt(2 * M_PI * P);
        }
        // back_stress = std::max(back_stress, -0.01);

        double rss_source = tau_applied(Kapp, r_source) + back_stress - tau_image(r_source);

        if ( rss_source >= tau_nuc + tau_friction && Nd < maxNumDis) 
        {
            // Nucleate a new dislocation
            Dislocation nucleatedDis = {nextDisID, r_source+b, 0.0};
            nextDisID++;
            disArr.push_back(std::make_shared<Dislocation>(nucleatedDis));
            std::cout << "Nucleating dis " << nucleatedDis.getId() << ", rss_source = " << rss_source << ", Ktip [MPa m^0.5] = "
            << Kapp-KD << ", increment = " << kInc << std::endl;
            ++Nd;
        }
        double dt = 10.0;
        double dx = dxMax;
        if (Nd > 0) 
        {
            std::vector<double> currP(Nd), tau_int(Nd, 0.0);
            std::vector<int> disID;
            double KD = 0.0;
            for (int i = 0; i < Nd; ++i) 
            {
                currP[i] = disArr[i]->getPosition();
                double dxDis = (i+1<Nd)? (currP[i] - disArr[i+1]->getPosition())/2.0 : dxMax;
                dx = std::min(dx, dxDis);
                if(dx < 0.0) std::exit(0);
                disID.push_back(disArr[i]->getId());
                KD += mu * b / sqrt(2 * M_PI * currP[i]);
                for (int j = 0; j < Nd; ++j) 
                {
                    if (i != j) tau_int[i] += tau_interaction(currP[i], disArr[j]->getPosition());
                }
            }
            std::vector<double> tau_app(Nd), tau_im(Nd), rss(Nd), currV(Nd, 0.0),
                athermal(Nd, false);
            double Vmax = 1e-3;
            for (int i = 0; i < Nd; ++i) 
            {
                // tau_app[i] = tau_applied(Kapp-KD, currP[i]);
                tau_app[i] = tau_applied(Kapp, currP[i]);
                tau_im[i] = tau_image(currP[i]);
                rss[i] = tau_app[i] + tau_int[i] - tau_im[i];
                rss[i] = rss[i]<=tau_friction? 0.0 : rss[i] - tau_friction;
                // auto vel = mobilityLaw.velocity_Q(rss[i], T);
                // auto vel = mobilityLaw.velocity_Inf(rss[i], T);
                auto vel = mobilityLaw.velocity(rss[i], T);
                currV[i] = vel.second;
                athermal[i] = vel.first;
                Vmax = std::max(Vmax, std::abs(currV[i]));
            }
            
            dt = Vmax>shearWaveFraction? dx/Vmax : dx/shearWaveFraction;

            std::vector<double> newP(Nd);
            std::vector<int> toRemove;
            for (int i = 0; i < Nd; ++i)
            {
                newP[i] = currP[i] + currV[i] * dt;
                // if (newP[i] <= crack_tip)
                if (newP[i] <= r_source + b)
                {
                    newP[i] = r_source + b;
                    if(currV[i] < 0.0)
                    {
                        toRemove.push_back(i);
                    }
                }
                // if(currV[i] < 0.0)
                // {
                //     toRemove.push_back(i);
                //     newP[i] = r_source; // Reset position to source if moving towards crack tip
                // }
            }
            bool isDescending = std::is_sorted(newP.begin(), newP.end(), std::greater<double>());
            if(!isDescending)
            {
                std::cout << "Dislocations not sorted in descending order, sorting..." << std::endl;
                for (int i = 0; i < Nd; ++i)
                {
                    std::cout << "Dislocation " << i << ": " << newP[i] << "velocity = "
                    << currV[i] <<std::endl;
                } 
                break;
            }

            for (int i = 0; i < Nd; ++i) 
            {
                disArr[i]->setPosition(newP[i]);
                disArr[i]->setVelocity(currV[i]);
                disArr[i]->setAthermal(athermal[i]);
                disArr[i]->setRss(rss[i]);
            }

            for (int i = toRemove.size() - 1; i >= 0; --i) 
            {
                disArr.erase(disArr.begin() + toRemove[i]);
                --Nd;
            }
            if(toRemove.size() > 0)
            {
                std::cout << "Annihilated " << toRemove.size() << " dislocations." << std::endl;
            }

        } 
        else 
        {
            dt = dxMax/shearWaveFraction;
        }

        // Store positions for plotting
        if (kInc % outputInterval == 0) 
        {
            outfile << time << ", " << Nd << ", " << rss_source << ", " << Kapp * unitSIF / 1e6
            << ", " << (Kapp-KD) * unitSIF / 1e6 << ", " << back_stress << ", " << dx << ", ";
            if(disArr.size() > 0) 
            {
                outfile << disArr[0]->getPosition();
                outfile << ", " << disArr[0]->getVelocity();
                outfile << ", " << disArr[0]->getRss();
            }
            else
            {
                outfile << 0.0;
                outfile << ", " << 0.0;
                outfile << ", " << 0.0;
            }
            outfile << std::endl;
            std::ofstream fileDis("output/dislocation_" + std::to_string(outFileIndex)+".csv", std::ios::trunc);
            for (int i = 0; i < Nd; ++i)
            {
                fileDis << disArr[i]->getId() << ", " << disArr[i]->isAthermal() << ", ";
                fileDis << std::scientific << std::setprecision(9)<< disArr[i]->getPosition() << ", "
                        << disArr[i]->getVelocity() << std::endl;
            }

            // if(outFileIndex == 20)
            // {
            //     std::cout << "dislocation "<< Nd-1 << ", P = " << disArr[Nd-1]->getPosition() << std::endl;
            //     break;
            // }

            outFileIndex++;
            fileDis.close();
        }

        time += dt;
        // if ((Kapp-KD) * unitSIF >= KG_SI)
        // {
        //     std::cout << "Ktip reached fracture toughness, stopping simulation." << std::endl;
        //     break;
        // }
    }

    std::cout << "Number of dislocations = " << Nd << std::endl;
    std::cout << "current SIF [MPa m^0.5] = " << Kapp * unitSIF / 1e6 << std::endl;
    std::cout << "Annihilated dislocations = " << nextDisID - Nd << std::endl;
    outfile.close();

    return 0;
}
