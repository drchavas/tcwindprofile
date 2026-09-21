# tc_pmin_estimatefromR34kt.py

## Predict tropical cyclone Pmin from Vmax, R34kt, latitude, translation speed, and environmental pressure (or pressure of outermost closed isobar)
### Source: Chavas D.R., Knaff J.A. and P. Klotzbach  (2025). A Simple Model for Predicting Tropical Cyclone Minimum Central Pressure from Intensity and Size. Wea. For., 40(2), pp.333-346
### https://doi.org/10.1175/WAF-D-24-0031.1
### Simple: Eq 5 to estimate dP; use Penv (or Eq 6 from Poci) to estimate Pmin.

### Cite this code: Chavas, D. R. (2025). Chavas Knaff Klotzbach 2024 Pmin prediction model. Purdue University Research Repository. doi:10.4231/NKX4-MM81

import numpy as np


def predict_Pmin_from_R34kt(
    VmaxNHC_ms: float,
    R34ktmean_km: float,
    lat: float,
    Vtrans_ms: float,
    Penv_mb: float,
    verbose: bool = False
) -> float:
    
    """
    # Predict Pmin (mb) from NHC point‐max wind speed (m/s), translation speed (m/s),
    # R34kt mean radius (km), latitude (deg), and environmental pressure (mb)
    # following CKK25
    # verbose: if True, print a status message (default False: silent, so the function can be called in a loop)
    # """

    # Calculate final predictors
    R34ktmean_m = R34ktmean_km * 1000
    Vmaxmean_ms = VmaxNHC_ms - 0.55 * Vtrans_ms  #Eq 2 of CKK25 -- simple method to estimate azimuthal-mean Vmax from NHC point-max; factor originally from Lin and Chavas (2012)
    omeg = 7.292e-5  # [s-1]
    fcor = 2 * omeg * np.sin(np.radians(np.abs(lat)))
    
    # Define intermediate predictor
    halffcorR34kt_ms = 0.5 * fcor * R34ktmean_m
    
    # Eq. 5 of CKK25: Best fit model coefficients
    coefs = np.array([-6.60, -0.0127, -5.506, 109.013])
    dP_predict_mb = (
        coefs[0] + coefs[1] * Vmaxmean_ms**2 +
        coefs[2] * halffcorR34kt_ms +
        coefs[3] * halffcorR34kt_ms / Vmaxmean_ms
    )
    
    Pmin_predict_mb = Penv_mb + dP_predict_mb   #Eq 3 of CKK25
    # print("Vmaxmean_ms =", Vmaxmean_ms,' m/s')
    # print("R34ktmean_km =", R34ktmean_m/1000,'km =', R34ktmean_m/1000 / km_nautmi,'naut mi')
    # print("dP_predict_mb =", dP_predict_mb,' mb')
    # print("Pmin_predict_mb =", Pmin_predict_mb,' mb (Penv =', Penv_mb,' mb)')

    if verbose:
        print("Returning minimum central sea-level pressure [mb]")
    return Pmin_predict_mb, dP_predict_mb