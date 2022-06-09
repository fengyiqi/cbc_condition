# Non-dimensionalization of The CBC Condition

## Generate Target CBC Spectra

Comete-Bellot and Corrsin performed a [turbulent experiment](https://acoustique.ec-lyon.fr/publi/comte-bellot_jfm71.pdf) at 1971. The recorded data have been to validate numerical methods in many papers.

The 5.08cm grid setup is used in this file. The collected kinetic energy at different wavelength are listed as (also see table 3 of their paper)

| $k, cm^{-1}$ | $\frac{tU_0}{M}=42$ | $\frac{tU_0}{M}=98$ | $\frac{tU_0}{M}=171$ |
|--------------|---------------------|---------------------|----------------------|
| 0.15         | -                   | -                   | 49.7                 |
| 0.20         | 129                 | 106                 | 92                   |
| 0.25         | 230                 | 196                 | 120                  |
| 0.30         | 322                 | 195                 | 125                  |
| 0.40         | 435                 | 202                 | 98                   |
| 0.50         | 457                 | 168                 | 81.5                 |
| 0.70         | 380                 | 127                 | 60.2                 |
| 1.00         | 270                 | 79.2                | 39.4                 |
| 1.50         | 168                 | 47.8                | 24.1                 |
| 2.00         | 120                 | 34.6                | 16.5                 |
| 2.50         | 89                  | 28.6                | 12.5                 |
| 3.00         | 70.3                | 23.1                | 9.12                 |
| 4.00         | 47                  | 14.3                | 5.62                 |
| 6.00         | 24.7                | 5.95                | 1.69                 |
| 8.00         | 12.6                | 2.23                | 0.52                 |
| 10.00        | 7.42                | 0.9                 | 0.161                |
| 12.50        | 3.96                | 0.363               | 0.052                |
| 15.00        | 2.33                | 0.162               | 0.0141               |
| 17.50        | 1.34                | 0.066               | -                    |
| 20.00        | 0.8                 | 0.033               | -                    | 

where $U_0=10m/s$ and $M=5.08cm$. Hence, we could compute the end time as 0.21336, 0.49784 and 0.86868, respectively. The Taylor microscale Reynolds numbers $Re_{\lambda} = u^`\lambda/\nu$ are 71.6, 65.1, 60.7, respectively. The grid Reynolds number of the experiment is $Re=34000$, then we have the dynamic viscosity about 0.000015.

In the simulation, the velocity field should be initialized based on the energy spectra at $\frac{tU_0}{M}=42$, and then is evolved to $\frac{tU_0}{M}=98$ and $\frac{tU_0}{M}=171$. If we set the start time at $\frac{tU_0}{M}=42$ as 0.0, the end time for next two state should be 0.28448 and 0.65532, respectively.

To simulate this isotropic flow in a $(2\pi)^3$ domain, the length, velocity and time should be scaled by the reference length $L_{ref}=10.8M/2\pi$, the reference velocity $U_{ref}=\sqrt{3/2}*22.2cm/s$ and the reference time $t_{ref}=L_{ref}/U_{ref}$. 

See [exp_data.ipynb](./exp_data.ipynb) for more details.

## Run Velocity Field Generation Code

1. Copy the generated target spectra into `double EnergyAtWavenumberGP(const int setup, const int wavenumber)` in `src/energy_spectrum.h`
2. make `build` folder and enter
3. `cmake ..`
4. `make`
5. `./Spectrum_init`
   
Finally you can find a `.csv` file in which u, v and w are stored. 