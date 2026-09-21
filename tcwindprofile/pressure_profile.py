# tcwindprofile/pressure_profile.py
# Calculate the pressure profile from an input wind profile that matches the input minimum pressure

# Integrates the gradient wind balance equation, assumes density is constant.
# Allowing density to vary with pressure is also doable (since P decreases moving towards the center),
# but the math works out to be the same as simply rescaling the constant-density solution by a constant
# Default density value is 1.15, but exact value does NOT matter anyways since profile
# is rescaled by a constant in the end to match Pmin

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def pressure_profile_calcfromwindprofile(
    rr_km: np.ndarray,
    vv_ms: np.ndarray,
    Renv_km: float,
    Penv_mb: float,
    Pmin_mb: float,
    lat: float,
    rho_kgm3: float = 1.15,  # [kg/m³], default constant density
    plot=False,
    verbose: bool = False
) -> tuple[float, np.ndarray]:
    """
    Estimate minimum pressure using gradient wind balance (GWB).
    
    Parameters
    ----------
    rr_km : ndarray
        Radius [m]
    vv_ms : ndarray
        Azimuthal wind speed [m/s]
    Renv_km : float
        Environmental radius where Penv applies [m]
    Penv_mb : float
        Environmental pressure [mb/hPa]
    Pmin_mb : float
        Minimum central pressure [mb/hPa]
    lat : float
        Latitude [deg]
    rho_kgm3 : float
        Air density (constant) [kg/m³]
    verbose : bool
        If True, print a status message (default False: silent, so the function can be called in a loop)

    Returns
    -------
    Pres_mb : ndarray
        Full pressure profile [mb]
    """
    # Check input range
    if np.nanmax(rr_km[~np.isnan(vv_ms)]) < Renv_km:
        raise ValueError("Renv lies beyond the input radial profile")

    # Unit conversion to MKS
    rr_m = rr_km * 1000 #[m]
    Renv_m = Renv_km * 1000 #[m]
    
    # Constants
    omega = 7.292e-5  # [s^-1]
    fcor = 2 * omega * np.sin(np.radians(abs(lat)))

    # High-resolution radius grid
    rr_hr = np.arange(0, Renv_m + 10, 10)  # [m]

    # Interpolate V to high-res grid
    vv_hr = interp1d(rr_m, vv_ms, bounds_error=False, fill_value=np.nan)(rr_hr)

    # Compute dP/dr via gradient wind balance
    term1 = fcor * vv_hr
    with np.errstate(divide='ignore', invalid='ignore'):
        term2 = np.where(rr_hr>0, vv_hr**2/rr_hr, 0.0)

    # term2 = vv_hr**2 / rr_hr
    dPdr = rho_kgm3 * (term1 + term2)

    # Integrate inward to get pressure deficit (dP defined POSITIVE)
    dr = np.nanmean(np.diff(rr_hr))
    dP_hr = np.flip(np.nancumsum(np.flip(dPdr * dr)))

    # Linear extrapolation for first point (if needed)
    if len(dP_hr) >= 3:
        dP_hr[0] = dP_hr[1] + (dP_hr[1] - dP_hr[2])

    # Interpolate back to original radius vector
    dP_interp = interp1d(rr_hr, dP_hr, bounds_error=False, fill_value=np.nan)
    dP_mb = dP_interp(rr_m) / 100  # convert Pa → mb (a.k.a. hPa)

    # Rescale pressure deficit by constant to match input dP
    dP_center_mb = Penv_mb - Pmin_mb  #dP defined POSITIVE
    dP_center_prof_mb = np.nanmax(dP_mb)  #dP defined POSITIVE
    # print(dP_center_prof_mb)
    dP_mb = dP_mb*(dP_center_mb/dP_center_prof_mb)
    
    # Compute full pressure profile
    pp_mb = Penv_mb - dP_mb
    


    if plot:
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #%% Make plot
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        # Figure dimensions: original 30 cm x 30 cm, now half: 15 cm x 15 cm (converted to inches)
        fig, ax = plt.subplots(figsize=(15/2.54, 15/2.54))
        ax.plot(rr_m / 1000, pp_mb, 'm-', linewidth=3)
        ax.plot(0, Pmin_mb, 'k.', markersize=15)
        ax.plot([Renv_km, Renv_km], [0, 10*Penv_mb], 'k--', linewidth=1)
        ax.plot(Renv_km, Penv_mb, 'k.', markersize=15)
        ax.set_xlabel('radius [km]')
        ax.set_ylabel('sea-level pressure [mb]')
        ax.set_title('Complete pressure profile', fontsize=12)
        ax.axis([0, 1.1 * Renv_km, 0.98 * Pmin_mb, 1020])
        
        # Annotate the top right corner with
        annotation = (f"Pmin = {Pmin_mb:.1f} mb\n"
                      f"Penv = {Penv_mb:.1f} mb")
        ax.text(0.1, 0.2, annotation, transform=ax.transAxes, ha='left', va='top',
                fontsize=10, bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))
        

        plt.savefig('pressureprofile.jpg', format='jpeg')
        plt.show()
        print("Made plot of pressure profile from input wind profile!")

    # Return pressure profile
    if verbose:
        print("Returning sea-level pressure vector [mb] along input radius vector [km]")
    return pp_mb
