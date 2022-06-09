#ifndef ENERGY_SPECTRUM_H
#define ENERGY_SPECTRUM_H

#include "user_definition.h"
#include <cmath>
#include <stdexcept>
#include "fftw3.h"
#include "utils.h"
#include <fstream>
#include <iostream>

constexpr double q2 = 3.0 / 2.0;

constexpr double ll = 55.0;
constexpr double const1 = 2.0 * PI / ll;
constexpr double const2 = const1 / (27.1893 * 27.1893);
constexpr double Amp = 1.0;

// Compute energy for CBC initial condition
double ComputeEnergyAtWavenumber(const double coef[7], const int wavenumber)
{
    double sum = 0.0;
    if (wavenumber >= 1.0 / 1000.0)
    {
        for (int a = 0; a < 7; a++)
            sum += coef[a] * std::pow(std::log(const1 * wavenumber), a);

        return const2 * std::exp(sum) / Amp;
    }
    else
        return 0.0;
}


// Implemented in BONSAI code where the polynomial interpolation was used
double EnergyAtWavenumberPolynomial(const int setup, const int wavenumber)
{

    double EofK = 0.0;
    switch (setup)
    {
    case 42: // Comte-Bellot and Corrosin spectrum M = 42
    {
        const double pol_cof[7] = {5.6102, -1.1236, -0.30961, 0.33172, -0.10959, -0.022320, 0.0066575};
        EofK = ComputeEnergyAtWavenumber(pol_cof, wavenumber);
        break;
    }
    case 98: // Comte-Bellot and Corrosin spectrum M = 98
    {
        const double pol_cof[7] = {0.43649e+01, -0.11793e+01, 0.67320e-01, 0.88554e-01, -0.13372e+00, 0.28385e-01, -0.99708e-02};
        EofK = ComputeEnergyAtWavenumber(pol_cof, wavenumber);
        break;
    }
    case 171: // Comte-Bellot and Corrosin spectrum M = 171
    {
        const double pol_cof[7] = {0.36567e+01, -0.11641e+01, -0.51571e-02, 0.38064e-02, -0.10484e+00, 0.12676e-01, -0.54562e-02};
        EofK = ComputeEnergyAtWavenumber(pol_cof, wavenumber);
        break;
    }
    case 1: // Hu's scale separation paper
    {
        EofK = std::pow(wavenumber, 4.0) * std::exp(- wavenumber * wavenumber / 4.0);
        break;
    }
    default:
        throw std::invalid_argument("invalid setup");
        break;
    }

    return EofK;
}

// Implemented with GP
double EnergyAtWavenumberGP(const int setup, const int wavenumber)
{

    switch (setup)
    {
    case 42: // Comte-Bellot and Corrosin spectrum M = 42
    {
        const double Eofk[56] = {
            0.0000, 0.0048, 0.0305, 0.0629, 0.0743, 0.0715, 0.0638, 0.0557, 0.0486, 0.0426,
            0.0377, 0.0337, 0.0304, 0.0276, 0.0253, 0.0233, 0.0215, 0.0200, 0.0187, 0.0175, 
            0.0164, 0.0154, 0.0145, 0.0137, 0.0130, 0.0123, 0.0117, 0.0111, 0.0106, 0.0101, 
            0.0096, 0.0092, 0.0088, 0.0084, 0.0080, 0.0077, 0.0074, 0.0071, 0.0068, 0.0065, 
            0.0062, 0.0060, 0.0058, 0.0055, 0.0053, 0.0051, 0.0049, 0.0047, 0.0046, 0.0044, 
            0.0042, 0.0041, 0.0039, 0.0038, 0.0037, 0.0035
            };
        return Eofk[wavenumber];
    }
    default:
        throw std::invalid_argument("invalid setup");
        break;
    }
}


void InitializeVelocity(const int setup)
{
    srand(1);
    const int maxWN = std::rint(std::sqrt(CutoffWN * CutoffWN + (X2_inner / 2) * (X2_inner / 2) + (X3_inner / 2) * (X3_inner / 2)));
    
    int k1[X1_inner + 1], k2[X2_inner], k3[X3_inner];
    InitWavenumbers(k1, k2, k3);

    double TargetEk[maxWN];
    for (unsigned int i = 0; i < maxWN; i++)
        TargetEk[i] = EnergyAtWavenumberGP(setup, i);

    double TargetEkin_total = 0.0;
    double Ekin_tot_test = 0.0;
    double Ekin_tot_testBACK = 0.0;

    for (int i = 0; i < CutoffWN; i++)
        TargetEkin_total += TargetEk[i];

    double spectrum_Ek_test[maxWN], spectrum_Ek_testBACK[maxWN];
    const int cells_cmplx = X2_inner * X3_inner * (X1_inner / 2 + 1);

    double *u_vel = new double[Block_Inner_Cell_Size];
    double *v_vel = new double[Block_Inner_Cell_Size];
    double *w_vel = new double[Block_Inner_Cell_Size];

    // spectral space velocity vectors
    fftw_complex *u_1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * cells_cmplx);
    fftw_complex *u_2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * cells_cmplx);
    fftw_complex *u_3 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * cells_cmplx);
    // row major order:
    fftw_plan plan_backward_1 = fftw_plan_dft_c2r_3d(X3_inner, X2_inner, X1_inner, u_1, u_vel, FFTW_ESTIMATE);
    fftw_plan plan_backward_2 = fftw_plan_dft_c2r_3d(X3_inner, X2_inner, X1_inner, u_2, v_vel, FFTW_ESTIMATE);
    fftw_plan plan_backward_3 = fftw_plan_dft_c2r_3d(X3_inner, X2_inner, X1_inner, u_3, w_vel, FFTW_ESTIMATE);

    fftw_complex *u_cmplx = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * cells_cmplx);
    fftw_complex *v_cmplx = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * cells_cmplx);
    fftw_complex *w_cmplx = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * cells_cmplx);

    fftw_plan plan_forward_1 = fftw_plan_dft_r2c_3d(X3_inner, X2_inner, X1_inner, u_vel, u_cmplx, FFTW_ESTIMATE);
    fftw_plan plan_forward_2 = fftw_plan_dft_r2c_3d(X3_inner, X2_inner, X1_inner, v_vel, v_cmplx, FFTW_ESTIMATE);
    fftw_plan plan_forward_3 = fftw_plan_dft_r2c_3d(X3_inner, X2_inner, X1_inner, w_vel, w_cmplx, FFTW_ESTIMATE);

    double Ekin_real_tot = 0.0;
    int count = 0;
    // main loop
    do
    {
        for (int k = 0; k < X3_inner; k++)
            for (int j = 0; j < X2_inner; j++)
                for (int i = 0; i < CutoffWN; i++)
                {
                    const int index = i + CutoffWN * (j + X2_inner * k);

                    const double wavenumber_r2 = k1[i] * k1[i] + k2[j] * k2[j] + k3[k] * k3[k];
                    const double wavenumber_r = std::sqrt(wavenumber_r2);
                    const int wavenumber_i = std::rint(wavenumber_r);

                    double f = 0.0;
                    if (wavenumber_r != 0.0)
                    {
                        f = std::sqrt(TargetEk[wavenumber_i] / (4.0 * PI) / wavenumber_r2); // Rogallo et al.
                        if (i != 0 && i != maxWN - 1)
                        {
                            double Amplification = std::sqrt(AmplificationCoefficient(setup)); // Possible for spectrum 55 and 64^3 grid.
                            f = f * Amplification;
                        }
                    }
                    const double theta1 = 2.0 * PI * randomfloat();
                    const double theta2 = 2.0 * PI * randomfloat();
                    const double phi = 2.0 * PI * randomfloat();

                    const double alpha_r = f * std::cos(theta1) * std::cos(phi);
                    const double alpha_i = f * std::sin(theta1) * std::cos(phi);
                    const double beta_r = f * std::cos(theta2) * std::cos(phi);
                    const double beta_i = f * std::sin(theta2) * std::cos(phi);

                    const double k12 = std::sqrt(k1[i] * k1[i] + k2[j] * k2[j]);
                    if (k12 == 0.0)
                    {
                        u_1[index][0] = alpha_r;
                        u_1[index][1] = alpha_i;
                        u_2[index][0] = beta_r;
                        u_2[index][1] = beta_i;
                        u_3[index][0] = 0.0;
                        u_3[index][1] = 0.0;
                    }
                    else
                    {
                        u_1[index][0] = (alpha_r * wavenumber_r * k2[j] + beta_r * k1[i] * k3[k]) / (wavenumber_r * k12);
                        u_1[index][1] = (alpha_i * wavenumber_r * k2[j] + beta_i * k1[i] * k3[k]) / (wavenumber_r * k12);

                        u_2[index][0] = (beta_r * k2[j] * k3[k] - alpha_r * wavenumber_r * k1[i]) / (wavenumber_r * k12);
                        u_2[index][1] = (beta_i * k2[j] * k3[k] - alpha_i * wavenumber_r * k1[i]) / (wavenumber_r * k12);

                        u_3[index][0] = (beta_r * k12) / wavenumber_r;
                        u_3[index][1] = (beta_i * k12) / wavenumber_r;
                    }
                }
        for (int i = 0; i < maxWN; i++)
            spectrum_Ek_test[i] = 0.0;

        for (int k = 0; k < X3_inner; k++)
            for (int j = 0; j < X2_inner; j++)
                for (int i = 0; i < CutoffWN; i++)
                {
                    const int index = i + CutoffWN * (j + X2_inner * k);
                    const double wavenumber_r2 = k1[i] * k1[i] + k2[j] * k2[j] + k3[k] * k3[k];
                    const double wavenumber_r = std::sqrt(wavenumber_r2);
                    const int wavenumber_i = std::rint(wavenumber_r);

                    double Phi_k = u_1[index][0] * u_1[index][0] + u_1[index][1] * u_1[index][1] + u_2[index][0] * u_2[index][0] + u_2[index][1] * u_2[index][1] + u_3[index][0] * u_3[index][0] + u_3[index][1] * u_3[index][1];
                    if (i != 0 && i != maxWN - 1)
                        Phi_k = Phi_k * 2.0;

                    spectrum_Ek_test[wavenumber_i] += 0.5 * Phi_k;
                }

        Ekin_tot_test = 0.0;
        for (int i = 0; i < CutoffWN; i++)
            Ekin_tot_test += spectrum_Ek_test[i];

        // TRANSFORM komplex velocities to physical space -- START
        for (int k = 0; k < X3_inner; k++)
            for (int j = 0; j < X2_inner; j++)
                for (int i = 0; i < X1_inner; i++)
                {
                    const int index = i + X1_inner * (j + X2_inner * k);
                    u_vel[index] = 0.0;
                    v_vel[index] = 0.0;
                    w_vel[index] = 0.0;
                }
        // Fourier transform
        fftw_execute(plan_backward_1);
        fftw_execute(plan_backward_2);
        fftw_execute(plan_backward_3);

        fftw_execute(plan_forward_1);
        fftw_execute(plan_forward_2);
        fftw_execute(plan_forward_3);

        for (int k = 0; k < X3_inner; k++)
            for (int j = 0; j < X2_inner; j++)
                for (int i = 0; i < CutoffWN; i++)
                {
                    int index = i + CutoffWN * (j + X2_inner * k);
                    for (int n = 0; n < 2; n++)
                    {
                        u_cmplx[index][n] /= (Block_Inner_Cell_Size);
                        v_cmplx[index][n] /= (Block_Inner_Cell_Size);
                        w_cmplx[index][n] /= (Block_Inner_Cell_Size);
                    }
                }
        for (int i = 0; i < maxWN; i++)
            spectrum_Ek_testBACK[i] = 0.0;

        for (int k = 0; k < X3_inner; k++)
            for (int j = 0; j < X2_inner; j++)
                for (int i = 0; i < CutoffWN; i++)
                {
                    int index = i + CutoffWN * (j + X2_inner * k);
                    double wavenumber_r2 = k1[i] * k1[i] + k2[j] * k2[j] + k3[k] * k3[k];
                    double wavenumber_r = std::sqrt(wavenumber_r2);
                    int wavenumber_i = std::rint(wavenumber_r);

                    double Phi_k = u_cmplx[index][0] * u_cmplx[index][0] + u_cmplx[index][1] * u_cmplx[index][1] + v_cmplx[index][0] * v_cmplx[index][0] + v_cmplx[index][1] * v_cmplx[index][1] + w_cmplx[index][0] * w_cmplx[index][0] + w_cmplx[index][1] * w_cmplx[index][1];

                    if (i != 0 && i != maxWN - 1)
                        Phi_k = Phi_k * 2.0;

                    spectrum_Ek_testBACK[wavenumber_i] += 0.5 * Phi_k;
                }

        Ekin_tot_testBACK = 0.0;
        for (int i = 0; i < CutoffWN; i++)
            Ekin_tot_testBACK += spectrum_Ek_testBACK[i];

        Ekin_real_tot = 0.0;
        for (int k = 0; k < X3_inner; k++)
            for (int j = 0; j < X2_inner; j++)
                for (int i = 0; i < X1_inner; i++)
                {
                    const int index = i + (X1_inner) * (j + X2_inner * k);
                    Ekin_real_tot += 0.5 * (u_vel[index] * u_vel[index] + v_vel[index] * v_vel[index] + w_vel[index] * w_vel[index]);
                }

        Ekin_real_tot /= (X1_inner * X2_inner * X3_inner); // Normalization

        // generating output file and appending to it

        std::cout << "------------------------------------------" << std::endl;
        std::cout << "CBC Mode = " << setup << std::endl;
        std::cout << "Iterations passed: " << count << std::endl;
        std::cout << "TargetEkin_total: " << TargetEkin_total << std::endl;

        std::cout << "Ekin_tot_test: " << Ekin_tot_test << std::endl;
        std::cout << "Ekin_real_tot: " << Ekin_real_tot << std::endl;
        std::cout << "Ekin_tot_testBACK: " << Ekin_tot_testBACK << std::endl;
        std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        count++;
    } while (std::fabs(Ekin_real_tot - TargetEkin_total) > 0.00001 && count <= 10000);

    // for file name
    std::string setup_suffix = std::to_string(setup);

    std::string spectrum_file_name = "energy_spec_at_ini_flow_ini_" + setup_suffix + ".dat";
    std::ofstream out_spec(spectrum_file_name);
    // defining header for tecplot(plot software)
    out_spec << "variables = k, E(k)<sub>imposed</sub>, E(k)<sub>is</sub>, E(k)<sub>real</sub>"
             << "\n";
    out_spec << "ZONE T='" << 0.0 << "', I= " << CutoffWN << ", DATAPACKING=POINT"
             << "\n";
    out_spec << "SOLUTIONTIME = " << 0.0 << std::endl;

    for (int i = 0; i < CutoffWN; i++)
        out_spec << i << " " << TargetEk[i] << " " << spectrum_Ek_test[i] << " " << spectrum_Ek_testBACK[i] << " ";
    out_spec << "\n";
    out_spec.close();

    std::string file_name = "velocity_field_" + setup_suffix + ".csv";
    std::ofstream csv_file(file_name);
    csv_file << "u,v,w" << std::endl;
    csv_file.setf(std::ios::fixed);
    csv_file.precision(16);
    for (unsigned int i = 0; i < Block_Inner_Cell_Size; i++)
        csv_file << u_vel[i] << "," << v_vel[i] << "," << w_vel[i] << std::endl;
    csv_file.close();

    delete[] u_vel;
    delete[] v_vel;
    delete[] w_vel;

    fftw_destroy_plan(plan_backward_1);
    fftw_destroy_plan(plan_backward_2);
    fftw_destroy_plan(plan_backward_3);
    fftw_free(u_1);
    fftw_free(u_2);
    fftw_free(u_3);
    fftw_destroy_plan(plan_forward_1);
    fftw_destroy_plan(plan_forward_2);
    fftw_destroy_plan(plan_forward_3);
    fftw_free(u_cmplx);
    fftw_free(v_cmplx);
    fftw_free(w_cmplx);
}

#endif