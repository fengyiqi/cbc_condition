#ifndef USER_DEFINITION_H
#define USER_DEFINITION_H

#include <stdexcept>

constexpr static int X1_inner = 64;
constexpr static int X2_inner = 64;
constexpr static int X3_inner = 64;
constexpr static int CutoffWN = (X1_inner / 2) + 1;
constexpr static double PI = 3.14159265358979323846;
constexpr static int Block_Inner_Cell_Size = 64 * 64 * 64;


double AmplificationCoefficient(const int setup)
{
    switch (setup)
    {
    case 42:
        return 2.02;
    case 98:
        return 2.06;
    case 171:
        return 2.2;
    case 1:
        return 2.4;
    default:
        throw std::invalid_argument("invalid setup");
        break;
    }
}

double DummyEnergyAtWavenumber1(const int setup)
{
    switch (setup)
    {
    case 42:
        return 0.0062;
    case 98:
        return 0.0043;
    case 171:
        return 0.0068;
    case 1:
        return 2.4;
    default:
        throw std::invalid_argument("invalid setup");
        break;
    }
}

#endif