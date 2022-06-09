#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include "user_definition.h"

void InitWavenumbers(int (&k1)[X1_inner + 1], int (&k2)[X2_inner], int (&k3)[X3_inner])
{
    for (int i = 0; i <= X1_inner / 2; i++)
        k1[i] = i;

    for (int i = 0; i <= X2_inner / 2; i++)
        k2[i] = i;
    for (int i = X2_inner / 2 + 1; i < X2_inner; i++)
        k2[i] = (i - X2_inner);

    for (int i = 0; i <= X3_inner / 2; i++)
        k3[i] = i;
    for (int i = X3_inner / 2 + 1; i < X3_inner; i++)
        k3[i] = (i - X3_inner);
}

double randomfloat()
{
    double scale = RAND_MAX + 1.;
    double base = rand() / scale;
    double fine = rand() / scale;
    return base + fine / scale;
}

double powB(double base, double ex)
{

    double ret;
    if ((base == 1.0) || (ex == 1.0))
    {
        ret = base;
    }
    else if (base == 0.0)
        ret = 0.0;
    else if (ex < 0)
    {
        ret = 1.0 / powB(base, -ex);
    }
    else
    {
        ret = std::exp(ex * std::log(std::abs(base)));
    }
    if (base >= 0.0)
        return (ret);
    else
        return -ret;
}



#endif