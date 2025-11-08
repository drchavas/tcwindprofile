import numpy as np
import matplotlib.pyplot as plt

def generate_wind_profile(Vmaxmean_ms, Rmax_km, R34ktmean_km, lat, plot=False):
    """
    Model described in Tao et al. 2025 GRL

    Analytic wind profile combining segments for:
    1) eye: r<Rmax (linear model);
    2) inner-core: Rmax to R34kt (linear-M model; Tao+ 2023 GRL);
    3) intermediate radii: R34kt to transition radius (modified Rankine model; Tao+ 2023 GRL, Klotzbach+ 2022 JGRA); and
    4) large radii: transition radius to outer radius (Ekman suction model; Emanuel 2004; Chavas+ 2015/2016 JAS).

    Inputs:
    - Vmaxmean_ms: Maximum wind speed at Rmax [m/s]
    - Rmax_km: Radius of maximum wind [km]
    - R34ktmean_km: Radius of 34-kt wind [km]
    - lat: Latitude [degrees]
    - plot: Optional flag to plot the profile

    Returns:
    - rr_km_interp: Radius vector [km]
    - vv_ms_interp: Wind speed profile [m/s]
    - R0_km: Outer radius where wind goes to zero [km]
    - Raa_km, Rba_km: Transition radii [km]
    - ir0, ira, irb: Indices of transitions (SL→mR, mR→E04, outer edge)
    """

    # Constants
    ms_kt = 0.5144444
    V34kt_ms = 34 * ms_kt
    omega = 7.2921e-5   #[s^{-1]}]; Earth's rotation rate
    fcor = np.abs(2 * omega * np.sin(np.radians(lat)))  #[s^{-1]}]; absolute value of Coriolis parameter at storm center latitude

    # Convert to meters
    Rmax_m = Rmax_km * 1000
    R34ktmean_m = R34ktmean_km * 1000
    dR = 100

    # Angular momentum
    Mmax = Vmaxmean_ms * Rmax_m + 0.5 * fcor * Rmax_m**2
    M34kt = V34kt_ms * R34ktmean_m + 0.5 * fcor * R34ktmean_m**2
    SL34kt = (M34kt / Mmax - 1) / (R34ktmean_m / Rmax_m - 1)

    Aa = 0.5 / fcor * (SL34kt * Mmax / Rmax_m)**2 + Mmax * (1 - SL34kt)
    Raa_m = SL34kt * Mmax / fcor / Rmax_m
    sqrt_term = np.sqrt(2 * Aa / fcor)
    chi = 1.5 #[-]; physical outer profile slope parameter; set to 1.5 constant (see Tao+ 2025, forthcoming); (= 2*C_d/w_cool)
    Rba_m = 0.5 * np.sqrt(chi * Aa * sqrt_term + Aa / (2 * fcor)) - 0.25 * sqrt_term
    R0a_m = 0.5 * np.sqrt(chi * Aa * sqrt_term + Aa / (2 * fcor)) + 0.75 * sqrt_term

    # Radius vector
    rr_full = np.arange(0, R0a_m + dR, dR)

    # Wind profiles
    V_eye = Vmaxmean_ms * (rr_full / Rmax_m)
    V_eye[rr_full > Rmax_m] = 0

    V_SL = ((SL34kt * (rr_full / Rmax_m - 1) + 1) * Mmax - 0.5 * fcor * rr_full**2) / rr_full
    alpha = 1.0
    V_mR = Aa * rr_full**(-alpha)
    V_E04 = (Aa - 0.5 * fcor * Rba_m**2) / rr_full + fcor * Rba_m - 0.5 * fcor * rr_full

    # Region indices
    ira = np.where(rr_full <= Raa_m)[0][-1]
    irb = np.where(rr_full <= Rba_m)[0][-1]
    ir0 = np.where(rr_full <= R0a_m)[0][-1]

    # Composite profile
    vv_ms_interp = np.zeros_like(rr_full)
    vv_ms_interp[rr_full < Rmax_m] = V_eye[rr_full < Rmax_m]
    vv_ms_interp[(rr_full >= Rmax_m) & (rr_full <= Raa_m)] = V_SL[(rr_full >= Rmax_m) & (rr_full <= Raa_m)]
    vv_ms_interp[(rr_full > Raa_m) & (rr_full <= Rba_m)] = V_mR[(rr_full > Raa_m) & (rr_full <= Rba_m)]
    vv_ms_interp[(rr_full > Rba_m) & (rr_full <= R0a_m)] = V_E04[(rr_full > Rba_m) & (rr_full <= R0a_m)]

    # Convert outputs to km
    rr_km_interp = rr_full / 1000
    R0_km = R0a_m / 1000
    Raa_km = Raa_m / 1000
    Rba_km = Rba_m / 1000

    # Optional plot
    if plot:
        fig, ax = plt.subplots(figsize=(15/2.54, 15/2.54))
        ax.plot(rr_km_interp, vv_ms_interp, 'm-', linewidth=3)
        # ax.axvline(Raa_km, linestyle='--', color='b', label='Raa')
        # ax.axvline(Rba_km, linestyle='--', color='r', label='Rba')
        # ax.axvline(R0_km, linestyle='--', color='g', label='R0')
        ax.plot(Rmax_km, Vmaxmean_ms, 'ko', label='(Rmax, Vmax)')
        ax.plot(R34ktmean_km, V34kt_ms, 'ko', label='(R34kt, 34 kt)')
        ax.set_xlabel('Radius [km]')
        ax.set_ylabel('Azimuthal wind speed [m/s]')
        ax.set_title('Complete wind profile')
        ax.set_xlim(0, R0_km * 1.1)
        ax.set_ylim(0, Vmaxmean_ms * 1.1)
        ax.grid(True)
        # ax.legend()
        annotation = (f"Inputs:\n"
                      f"Vmax_mean = {Vmaxmean_ms:.1f} m/s\n"
                      f"Rmax = {Rmax_km:.1f} km\n"
                      f"(R34kt_mean, V34kt_ms) = ({R34ktmean_km:.1f} km, {V34kt_ms:.1f} m/s)\n"
                      f"lat = {lat:.1f}°N\n\n"
                      f"Output:\n"
                      f"R0_mean = {R0_km:.1f} km")
        ax.text(0.95, 0.95, annotation, transform=ax.transAxes, ha='right', va='top',
                fontsize=10, bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))
        plt.savefig('windprofile.jpg', format='jpeg')
        plt.show()
        print("Made plot of wind profile!")

    return rr_km_interp, vv_ms_interp, R0_km
