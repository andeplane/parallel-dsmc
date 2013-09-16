#include <unitconverter.h>

UnitConverter::UnitConverter()
{
    m0 = 1.66053886e-27;  // SI
    L0 = 1e-6;            // SI
    E0 = 1.65088e-21;     // SI
    E0ev = 1.0318e-2;     // eV
    kb = 1.3806503e-23;   // SI

    t0 = L0*sqrt(m0/E0);
    F0 = E0/L0;
    T0 = E0/kb;
    P0 = m0/(t0*t0*L0);
    v0 = L0/t0;
    a0 = v0/t0;
    visc0 = P0*t0;
    diff0 = L0*L0/t0;
    perm0 = L0*L0;
    number_density0 = 1.0/(L0*L0*L0);
}

double UnitConverter::pressure_to_SI(double P) { return P0*P; }
double UnitConverter::pressure_from_SI(double P) { return P/P0; }

double UnitConverter::temperature_to_SI(double T) { return T0*T; }
double UnitConverter::temperature_from_SI(double T) { return T/T0; }

double UnitConverter::mass_to_SI(double m) { return m0*m; }
double UnitConverter::mass_from_SI(double m) { return m/m0; }

double UnitConverter::length_to_SI(double L) { return L0*L; }
double UnitConverter::length_from_SI(double L) { return L/L0; }

double UnitConverter::force_to_SI(double F) { return F0*F; }
double UnitConverter::force_from_SI(double F) { return F/F0; }

double UnitConverter::energy_to_SI(double E) { return E0*E; }
double UnitConverter::energy_from_SI(double E) { return E/E0; }

double UnitConverter::energy_to_eV(double E) { return E0ev*E; }
double UnitConverter::energy_from_eV(double E) { return E/E0ev; }

double UnitConverter::time_to_SI(double t) { return t0*t; }
double UnitConverter::time_from_SI(double t) { return t/t0; }

double UnitConverter::velocity_to_SI(double v) { return v*v0; }
double UnitConverter::velocity_from_SI(double v) { return v/v0; }

double UnitConverter::acceleration_to_SI(double a) { return a*a0; }
double UnitConverter::acceleration_from_SI(double a) { return a/a0; }

double UnitConverter::viscosity_to_SI(double v) { return v*visc0; }
double UnitConverter::viscosity_from_SI(double v) { return v/visc0; }

double UnitConverter::permeability_to_SI(double k) { return k*perm0; }
double UnitConverter::permeability_from_SI(double k) { return k/perm0; }

double UnitConverter::diffusion_to_SI(double d) { return d*diff0; }
double UnitConverter::diffusion_from_SI(double d) { return d/diff0; }

double UnitConverter::number_density_to_SI(double d) { return d*number_density0; }
double UnitConverter::number_density_from_SI(double d) { return d/number_density0; }


