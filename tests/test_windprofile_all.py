# test_windprofile_all.py

## Create a fast and robust radial profile of the tropical cyclone rotating wind from inputs Vmax, R34kt, latitude, and Vtrans.


import os
import sys

# Add parent directory to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import math

############################################################
# NHC/Best Track Operational Inputs
VmaxNHC_kt = 97.19  #100 [kt]; NHC storm intensity (point-max wind speed)
Vtrans_kt = 0    #20 [kt]; storm translation speed, usually estimated from adjacent track points; used to estimate azimuthal-mean Vmax (Vmaxmean_ms = VmaxNHC_ms - 0.55*Vtrans_ms)
lat = 20  #20 [degN]; default 20N; storm-center latitude
R34ktNHCquadmax_nautmi = 127.04 #(135 + 150 + 145 + 150) / 4 #average NHC R34kt radius (here 4 quadrants)
                                                        # these are officially the MAXIMUM radii of this wind speed in each quadrant;
                                                        # value is reduced by factor 0.85 within the code to estimate the mean radius (see Chavas Knaff Klotzbach 2025 for more info)
Penv_mb = 1008      #[mb]; environmental pressure, to create full pressure profile
## Default values: VmaxNHC_kt=100 kt, R34ktNHCquadmax_nautmi= 145.0 naut mi, lat = 20 --> unadjusted Rmax=38.1 km (sanity check)

# Rmax
Rmax_km = 30     #30 [km]; input value
#Rmax_km = None  #None: estimate Rmax from R34kt -- ref Chavas and Knaff 2022 WAF)
############################################################


################################################################
## Calculate wind and pressure profiles and associated data
"""
Full modeling pipeline:
- If no Rmax input: estimate Rmax from R34kt -- ref Chavas and Knaff 2022 WAF
- Estimate R0 from R34kt: ref Tao et al. (2025, GRL); approximate version of outer model of refs Emanuel 2004 / Chavas et al. 2015 JAS / Chavas and Lin 2016 JAS
- Generate wind profile: Analytic complete wind profile: ref Tao et al. (2025, GRL)
    1) eye: r<Rmax (linear model);
    2) inner-core: Rmax to R34kt (linear-M model; Tao+ 2023 GRL);
    3) intermediate radii: R34kt to transition radius (modified Rankine model; Tao+ 2023 GRL, Klotzbach+ 2022 JGRA); and
    4) large radii: transition radius to outer radius (Ekman suction model; Emanuel 2004; Chavas+ 2015/2016 JAS).
- Estimate Pmin: ref Chavas Knaff Klotzbach 2025 WAF
- Generate pressure profile that matches Pmin: ref Chavas Knaff Klotzbach 2025 WAF
"""
from tcwindprofile.windprofile_all import run_full_wind_model


run_kwargs = dict(
    VmaxNHC_kt=VmaxNHC_kt,
    Vtrans_kt=Vtrans_kt,
    R34kt_quad_max_nautmi=R34ktNHCquadmax_nautmi,
    lat=lat,
    Penv_mb=Penv_mb,
    plot=True,
)

# Only pass Rmax_km if user specified it
if Rmax_km is not None:
    run_kwargs["Rmax_km"] = Rmax_km
    mode_tag = "givenRmax"
else:
    mode_tag = "estRmax"

tc_wind_and_pressure_profile = run_full_wind_model(**run_kwargs)

if Rmax_km is None:
    print(f"[Estimated Rmax] Rmax = {tc_wind_and_pressure_profile['Rmax_km']:.1f} km")
else:
    print(f"[Given Rmax]     Rmax = {tc_wind_and_pressure_profile['Rmax_km']:.1f} km")

print(f"R0 = {tc_wind_and_pressure_profile['R0_km']:.1f} km")
print(f"Pmin = {tc_wind_and_pressure_profile['Pmin_mb']:.1f} hPa")
################################################################

################################################################
################################################################
## Plot that data
from tcwindprofile.plot_windprofile import plot_wind_and_pressure

def plot_profile(tc_profile, save_path):
    # unpack
    rr_km         = tc_profile['rr_km']
    vv_ms         = tc_profile['vv_ms']
    pp_mb         = tc_profile['pp_mb']
    Vmaxmean_ms   = tc_profile['Vmaxmean_ms']
    Rmax_km       = tc_profile['Rmax_km']
    R34ktmean_km  = tc_profile['R34ktmean_km']
    V34kt_ms      = tc_profile['V34kt_ms']
    R0_km         = tc_profile['R0_km']
    Pmin_mb       = tc_profile['Pmin_mb']

    # then:
    Renv_km = R0_km
    plot_wind_and_pressure(
        rr_km, vv_ms,
        Rmax_km, Vmaxmean_ms,
        R34ktmean_km, V34kt_ms,
        R0_km, lat,
        rr_km, pp_mb,
        Renv_km, Penv_mb, Pmin_mb,
        save_path=save_path
    )

plot_profile(tc_wind_and_pressure_profile, save_path=f'tc_wind_pressure_profiles_{mode_tag}')
################################################################
