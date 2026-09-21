# tc_rmax_estimatefromR34kt.py

## Predict tropical cyclone Rmax from Vmax, R34kt, and latitude
### Source: Chavas D.R. and J. A.. Knaff (2022). A simple model for predicting the tropical cyclone radius of maximum wind from outer size. Wea. For., 37(5), pp.563-579
### https://doi.org/10.1175/WAF-D-21-0103.1

### 3 simple steps: Eq 3+4+7 from the paper
### Plus optional final bias adjustment: Eq. 8 (applied only for Rmax>60 km; this adjustment is only significant for larger Rmax > 60 km  -- see Fig 8 of paper)

### Cite this code: Daniel Robert Chavas (2022). Code file for Rmax prediction model and relevant data used in the analyses (Chavas and Knaff 2022, Weather and Forecasting). (Version 1.1). Purdue University Research Repository. doi:10.4231/WMMS-XY76

import numpy as np

def predict_Rmax_from_R34kt(
    VmaxNHC_ms: float,
    R34ktmean_km: float,
    lat: float,
    verbose: bool = False
) -> float:
    """
    # Predict Rmax (km) from NHC point‐max wind speed (m/s),
    # R34kt mean radius (km), and latitude (deg), using CK22.
    # verbose: if True, print a status message (default False: silent, so the function can be called in a loop)
    # """
    
    R34ktmean_m = R34ktmean_km * 1000
    # Calculate M34kt (Eq 2 of CK22)
    omeg = 7.292e-5  # [s-1]
    fcor = 2 * omeg * np.sin(np.radians(np.abs(lat)))
    ms_kt = 0.5144444  # 1 kt = 0.514444 m/s
    V34kt = 34 * ms_kt  # [m/s]
    M34kt = R34ktmean_m * V34kt + 0.5 * np.abs(fcor) * R34ktmean_m**2
    
    # Define intermediate predictor
    halffcorR34kt_ms = 0.5 * fcor * R34ktmean_m
    
    # Estimate Mmax/M34kt (Eq 7 of CK22)
    # Model: Z = b*exp(beta_Vmax*(Vmax-V34kt)+beta_VfR*(Vmax-V34kt)*(0.5*f*R34kt))
    # Coefficients estimated from EBTK 2004-2020 (Eq 7 / Table 2 of CK22)
    coefs = {
        "b": 0.699,  # b
        "c": -0.00618,  # beta_Vmax
        "f": -0.00210  # beta_VfR
    }
    
    MmaxM34kt_predict = coefs["b"] * np.exp(
        coefs["c"] * (VmaxNHC_ms - V34kt) + coefs["f"] * (VmaxNHC_ms - V34kt) * halffcorR34kt_ms
    )
    
    # Solve for predicted Rmax (Mmax = Rmax*Vmax + 0.5*fcor*Rmax^2)
    # (Eq 3 of CK22)
    Mmax_predict = MmaxM34kt_predict * M34kt  # Mmax = MmaxM34kt * M34kt
    
    # (Eq 4 of CK22)
    Rmax_predict = (VmaxNHC_ms / fcor) * (
        np.sqrt(1 + (2 * fcor * Mmax_predict / (VmaxNHC_ms**2))) - 1
    )
    
    Rmax_predict_km = Rmax_predict / 1000
    km_nautmi = 1.852
    Rmax_predict_nautmi = Rmax_predict_km / km_nautmi
    # print("VmaxNHC_ms =", VmaxNHC_ms,'m/s')
    # print("R34ktmean_km =", R34ktmean_m/1000,'km =', R34ktmean_m/1000 / km_nautmi,'naut mi')
    # print('Rmax_predict_km=', Rmax_predict_km,'km =', Rmax_predict_nautmi,'naut mi')
    
    ## Apply optional final bias correction (Eq.)
    ## Do only if Rmax > 60 km; below that the uncertainties are same order as the value itself
    ## Also, for very small values (<9 km) it could give a negative Rmax which makes no sense
    if Rmax_predict_km > 60:
      Rmax_predict_km_biasadjgt60km = (1/0.76)*(Rmax_predict_km-9.02)
      Rmax_predict_nautmi_biasadjgt60km = Rmax_predict_km_biasadjgt60km / km_nautmi
      # print('Rmax_predict_km_biasadjgt60km=', Rmax_predict_km_biasadjgt60km,' km =', Rmax_predict_nautmi_biasadjgt60km,' naut mi')
      # print('Rmax estimation from R34kt: final bias adjustment applied for Rmax>60km, based on CK22')
    # else:
      # print('Rmax estimation from R34kt: no final bias adjustment (only done for Rmax>60km)')

    if verbose:
        print("Returning estimated radius of maximum wind [km]")
    return Rmax_predict_km
