# windprofile.py

## Create a fast and robust radial profile of the tropical cyclone rotating wind from inputs Vmax, R34kt, Rmax, and latitude.

#### This code uses a modified‐Rankine vortex between Rmax and R34kt and the E04 model beyond R34kt (and a quadratic profile inside the eye). It is very similar to the full physics-based wind profile model of Chavas et al. (2015) ([code here](http://doi.org/10.4231/CZ4P-D448)), but is simpler and much faster.

#### It is designed to guarantee that the profile fits both Rmax and R34kt and will be very close to the true outer radius (R0) as estimated by the full E04 outer solution. Hence, it is very firmly grounded in the known physics of the tropical cyclone wind field while also matching the input data. It is also guaranteed to be very well‐behaved for basically any input parameter combination.

#### Model basis:
#### Modified Rankine profile between Rmax and R34kt was shown to compare very well against high-quality subset of Atlantic Best Track database -- see Fig 8 of [Klotzbach et al. (2022, JGR-A)](https://doi.org/10.1029/2022JD037030)
#### Physics-based non-convecting wind field profile beyond R34kt was shown to compare very well against entire QuikSCAT database -- see Fig 6 of [Chavas et al. (2015, JAS)](https://doi.org/10.1175/JAS-D-15-0014.1)
#### Quadratic in the eye (U-shape is common)

from tcwindprofile.tc_outer_radius_estimate import estimate_outer_radius
from tcwindprofile.tc_outer_windprofile import outer_windprofile


def generate_wind_profile(Vmaxmean_ms, Rmax_km, R34ktmean_km, lat, plot=False):
    

    import numpy as np
    import math
    import matplotlib.pyplot as plt
    from scipy.interpolate import interp1d
    
    # Set default plot properties: font size 12 (half of 24) and default font
    plt.rcParams.update({
        'font.size': 12,
        'lines.linewidth': 2
    })
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ### INPUTS: AZIMUTHAL-MEAN VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # If you have VmaxNHC (i.e. Best Track values): Vmaxmean_ms = VmaxNHC_ms - 0.55 * Vtrans_ms (Eq 2 of CKK25 -- simple method to estimate azimuthal-mean Vmax from NHC point-max; factor originally from Lin and Chavas (2012))
    # If you have R34ktNHCquadmax (i.e. Best Track values): R34ktmean_m = 0.85 * fac_R34ktNHCquadmax2mean (Eq 1 of CKK25 -- simple estimate of mean R34kt radius from NHC R34kt (which is maximum radius of 34kt); factor originally from DeMaria et al. (2009))
    
    #Vmaxmean_ms = 45.8            # [m/s]; default 45.8; AZIMUTHAL-MEAN maximum wind speed
    #Rmax_m = 38.1 * 1000          # [m]; default 38.1*1000; from bias-adjusted CK22
    #R34ktmean_m = 228.3 * 1000    # [m]; default R34kt = 228.3*1000; MEAN radius of 34kt wind speed
    #lat = 20                      # [degN]; default 20N; storm-center latitude;
    
    Rmax_m = Rmax_km * 1000
    R34ktmean_m = R34ktmean_km * 1000
    ## Default values: Vmaxmean_ms = 45.8 m/s, Rmax_m = 38.1 km, R34ktmean_m = 228.3 km, lat = 20 --> R0=1174.6 km (sanity check)
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ### CALCULATE WIND PROFILE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    # print("This code uses a modified‐Rankine vortex between Rmax and R34ktmean_m (default R34kt) and the E04 model beyond R34ktmean_m (and a quadratic profile inside the eye).")
    # print("It is designed to guarantee that the profile fits both Rmax and R34kt and will be very close to the true outer radius (R0) as estimated by the full E04 outer solution.")
    # print("It is also guaranteed to be very well‐behaved for basically any input parameter combination.")
    
    #%% Calculate some params
    ms_kt = 0.5144444             # 1 kt = 0.514444 m/s
    V34kt_ms = 34 * ms_kt            # [m/s]; outermost radius to calculate profile
    km_nautmi = 1.852
    omeg = 7.292e-5  # Earth's rotation rate
    fcor = 2 * omeg * math.sin(math.radians(abs(lat)))  # [s^-1]
    
    #%% Initial radius vector
    dr = 100  # [m]
    R0mean_initialtoolarge = 3000 * 1000
    rr_full = np.arange(0, R0mean_initialtoolarge + dr, dr)  # [m]
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%% STEP 1) Outer region wind profile (r>R34kt)
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%% STEP 1a) Very good estimate of outer edge of storm, R0mean (where v=0)
    #%% Analytic approximation of R0mean, from physical model of non-convecting wind profile (Emanuel 2004; Chavas et al. 2015 JAS)
    
    Cd = 1.5e-3  # [-]
    w_cool = 2 / 1000  # [m/s]
    beta = 1.35
    R0mean_dMdrcnstmod = estimate_outer_radius(
    R34ktmean_m=R34ktmean_m,
    V34kt_ms=V34kt_ms,
    fcor=fcor,
    Cd=Cd,
    w_cool=w_cool,
    beta=beta
    )
    
    R0_km = R0mean_dMdrcnstmod / 1000
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%% STEP 1b) Wind profile from R0mean --> R34kt
    #%% Exact solution to model of non-convecting wind profile (Emanuel 2004; Chavas et al. 2015 JAS)
    #%% Notes:
    #%% - It is a single integration inwards from R0. *Cannot* integrate out
    #%%   from R34kt, it may blow up (math is weird)
    #%% - Could do a simpler approximate profile instead, but this solution is
    #%%   just as fast doing some sort of curve fit, so might as well use exact
 
    # Returns:
    #   rr_full : radii [m]
    #   vv : wind speeds [m/s]
    rr_outer_m, vv_outer_ms = outer_windprofile(
    r0_m=R0mean_dMdrcnstmod,     # or your computed R0mean_dMdrcnstmod
    fcor=fcor,
    Cd=Cd,
    w_cool=w_cool,
    V34kt_ms=V34kt_ms,
    )
    
    # Interpolate to original radius vector, chop off NaNs beyond outer edge (which will extend a bit beyond r0)
    r0_plot = np.max(rr_outer_m)
    rr_full=rr_full[rr_full<=r0_plot]
    vv_E04approx = np.interp(rr_full, rr_outer_m, vv_outer_ms)

    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%% STEP 2) Inner region wind profile (r<=R34kt)
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    #%% Simple modified Rankine profile between R34kt and Rmax
    #%% With quadratic profile inside of Rmax
    alp_outer = np.log(Vmaxmean_ms / V34kt_ms) / np.log(Rmax_m / R34ktmean_m)
    vv_outsideRmax = np.full(rr_full.shape, np.nan)
    rr_outer_mask = rr_full > Rmax_m
    vv_outsideRmax[rr_outer_mask] = Vmaxmean_ms * (rr_full[rr_outer_mask] / Rmax_m)**alp_outer
    
    alp_eye = 1  # 1 = linear; 2 = quadratic profile in eye
    vv_eye = Vmaxmean_ms * (rr_full / Rmax_m)**alp_eye
    vv_MR = np.concatenate((vv_eye[rr_full <= Rmax_m], vv_outsideRmax[rr_full > Rmax_m]))
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%% STEP 3) Merge inner and outer wind profiles
    #%% Simple exponential smoother of (inner-outer) moving outwards from R34kt
    
    # No smoothing: direct merge at R34ktmean_m -- will *not* match dV/dr without smoothing
    vv_MR_E04_R34ktmean_nosmooth = np.concatenate((vv_MR[rr_full <= R34ktmean_m], vv_E04approx[rr_full > R34ktmean_m]))
    
    # Match dV/dr at R34kt: exponential smoother of (inner-outer) moving outwards from R34kt
    ii_temp = rr_full > R34ktmean_m
    rr_beyondR34ktmean = rr_full[ii_temp]
    vv_MR_beyondR34ktmean = vv_MR[ii_temp]
    vv_E04approx_beyondR34ktmean = vv_E04approx[ii_temp]
    R_transitionouteredge = R34ktmean_m + (R0_km*1000 - R34ktmean_m)* (3 / 3)   #set to R0; played with smaller radii but can give weird kinks
    def taper_cubic(r, Rin, Rout):
        """1 at Rin, 0 at Rout, derivative continuous at both, smooth cubic Hermite."""
        s = (r - Rin) / (Rout - Rin)
        s = np.clip(s, 0, 1)
        return 1 - (3*s*s - 2*s*s*s)
        #note: tried many other more complex options, but they only create additional "knees" in the curve that make things worse
        #tried: quintic (not better), cosine (nearly identical), cubic with s^p p<1 (sharper transition but now derivative may be discontinuous at Rin), beta with sharpness param (maintains derivative continuity, but is just bumpier)
    w = taper_cubic(rr_beyondR34ktmean, R34ktmean_m, R_transitionouteredge)

    vv_E04approx_adj = vv_MR_beyondR34ktmean * w + vv_E04approx_beyondR34ktmean * (1 - w)
    vv_MR_E04_R34ktmean = np.concatenate((vv_MR[rr_full <= R34ktmean_m], vv_E04approx_adj))
    
    # Interpolate to 1 km resolution
    rr_km = rr_full / 1000
    rr_km_interp = np.arange(0, rr_km.max(), 0.1)  # 0.1 km resolution
    vv_ms_interp = np.interp(rr_km_interp, rr_km, vv_MR_E04_R34ktmean)



    if plot:
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #%% Make plot
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Figure dimensions: original 30 cm x 30 cm, now half: 15 cm x 15 cm (converted to inches)
        fig, ax = plt.subplots(figsize=(15/2.54, 15/2.54))
        ax.plot(rr_km_interp, vv_ms_interp, 'm-', linewidth=3)
        ax.plot(rr_full / 1000, vv_MR_E04_R34ktmean_nosmooth, 'g:', linewidth=3)   #compare no smooth vs. smoothed
        ax.plot(Rmax_m / 1000, Vmaxmean_ms, 'k.', markersize=20)
        ax.plot(R34ktmean_m / 1000, V34kt_ms, 'k.', markersize=20)
        ax.plot(R0_km, 0, 'm.', markersize=20)
        ax.set_xlabel('radius [km]')
        ax.set_ylabel('azimuthal wind speed [m/s]')
        ax.axis([0, 1.1 * R0mean_dMdrcnstmod / 1000, 0, 1.1 * Vmaxmean_ms])
        # ax.axis([0, 600, 0, 1.1 * Vmaxmean_ms])
        ax.set_title('Complete wind profile', fontsize=12)
        
        # Annotate the top right corner with "Inputs:" and the values of Vmaxmean_ms, Rmax_m, (R34ktmean_m, V34kt_ms), and lat
        annotation = (f"Inputs:\n"
                      f"Vmax_mean = {Vmaxmean_ms:.1f} m/s\n"
                      f"Rmax = {Rmax_m/1000:.1f} km\n"
                      f"(R34kt_mean, V34kt_ms) = ({R34ktmean_m/1000:.1f} km, {V34kt_ms:.1f} m/s)\n"
                      f"lat = {lat:.1f}°N\n"
                      f"\n"
                      f"Output:\n"
                      f"R0_mean = {R0mean_dMdrcnstmod / 1000:.1f} km")
        ax.text(0.95, 0.95, annotation, transform=ax.transAxes, ha='right', va='top',
                fontsize=10, bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))
        
        plt.savefig('windprofile.jpg', format='jpeg')
        plt.show()
        print("Made plot of wind profile!")
               
    # Return radius [km], wind speed [m/s], true R0 [km], estimated R0 [km]
    print("Returning radius vector [km], wind speed vector [m/s], estimated outer radius [km]")
    return rr_km_interp, vv_ms_interp, R0_km


    
    
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--vmax_ms", type=float, default=45.8, help="Vmax in m/s")
    parser.add_argument("--rmax_km", type=float, default=38.1, help="Rmax in km")
    parser.add_argument("--r34_km", type=float, default=228.3, help="R34ktmean in km")
    parser.add_argument("--lat", type=float, default=20, help="Latitude in degrees")
    args = parser.parse_args()

    # Convert km to m
    rr_km, vv_ms, R0_km = generate_wind_profile(args.vmax_ms, args.rmax_km, args.r34_km, args.lat)

