# windprofile.py

## Create a fast and robust radial profile of the tropical cyclone rotating wind from inputs Vmax, R34kt, optional Rmax, and latitude.


# tcwindprofile/full_profile.py

from typing import Optional

def run_full_wind_model(
    VmaxNHC_kt: float,
    Vtrans_kt: float,
    R34kt_quad_max_nautmi: float,
    lat: float,
    Penv_mb: float,
    plot: bool = False,
    Rmax_km: Optional[float] = None
):
    """
    Full modeling pipeline:
    - If no Rmax input: estimate Rmax from R34kt -- ref Chavas and Knaff 2022 WAF
    - Estimate R0 from R34kt: approximate version of outer model ref Emanuel 2004 / Chavas et al. 2015 JAS / Chavas and Lin 2016 JAS
    - Generate wind profile: merge simple inner + outer models, ref Klotzbach et al. 2022 JGR-A / Chavas and Lin 2016 JAS
    - Estimate Pmin: ref Chavas Knaff Klotzbach 2025 WAF
    - Generate pressure profile that matches Pmin: ref Chavas Knaff Klotzbach 2025 WAF
    """
    # Convert input units
    ms_per_kt = 0.5144444
    km_per_nautmi = 1.852
    VmaxNHC_ms = VmaxNHC_kt * ms_per_kt
    Vtrans_ms = Vtrans_kt * ms_per_kt
    Vmaxmean_ms = VmaxNHC_ms - 0.55 * Vtrans_ms  #Lin and Chavas (2012), Chavas Knaff Klotzbach (2025)
    R34ktmean_km = 0.85 * R34kt_quad_max_nautmi * km_per_nautmi
    
    V34kt_ms = 17.49
    
    if Vmaxmean_ms < V34kt_ms:
        raise ValueError("Vmaxmean_ms cannot be < 34 kt")
    ###########################
    # 1) Estimate Rmax (CK22) unless provided by user
    if Rmax_km is None:
        from tcwindprofile.tc_rmax_estimatefromR34kt import predict_Rmax_from_R34kt
        Rmax_km = predict_Rmax_from_R34kt(
            VmaxNHC_ms=VmaxNHC_ms,
            R34ktmean_km=R34ktmean_km,
            lat=lat
        )

    ###########################
    # 2) Estimate wind profile + R0 (Tao et al. 2025 GRL; approximation to CLE15 model)
    # print(Vmaxmean_ms,Rmax_km,R34ktmean_km,lat,plot)
    from tcwindprofile.windprofile import generate_wind_profile
    # from tcwindprofile.tc_outer_radius_estimate import estimate_outer_radius
    # from tcwindprofile.tc_outer_windprofile import outer_windprofile
    
    rr_km, vv_ms, R0_km = generate_wind_profile(
        Vmaxmean_ms=Vmaxmean_ms,
        Rmax_km=Rmax_km,
        R34ktmean_km=R34ktmean_km,
        lat=lat,
        # plot=True
        plot= plot 
    )

    ###########################
    # 3) Estimate Pmin (CKK25)
    # print(VmaxNHC_ms,R34ktmean_km,lat,Vtrans_ms,Penv_mb)
    from tcwindprofile.tc_pmin_estimatefromR34kt import predict_Pmin_from_R34kt
    
    Pmin_estimate_mb, dP_estimate_mb = predict_Pmin_from_R34kt(
        VmaxNHC_ms=VmaxNHC_ms,
        R34ktmean_km=R34ktmean_km,
        lat=lat,
        Vtrans_ms=Vtrans_ms,
        Penv_mb=Penv_mb
    )

    ###########################
    # 4) Calculate pressure profile from wind field and matches Pmin
    from tcwindprofile.pressure_profile import pressure_profile_calcfromwindprofile                        # [deg]
    
    pp_mb = pressure_profile_calcfromwindprofile(
        rr_km=rr_km,
        vv_ms=vv_ms,
        Renv_km=R0_km,
        Penv_mb=Penv_mb,
        Pmin_mb=Pmin_estimate_mb,
        lat=lat,
        plot=False)
    
    ###########################
    # 5) Return final results
    return {
        "rr_km": rr_km,
        "vv_ms": vv_ms,
        "pp_mb": pp_mb,
        "Vmaxmean_ms": Vmaxmean_ms,
        "Rmax_km": Rmax_km,
        "V34kt_ms": V34kt_ms,
        "R34ktmean_km": R34ktmean_km,
        "R0_km": R0_km,
        "lat": lat,
        "Penv_mb": Penv_mb,
        "Pmin_mb": Pmin_estimate_mb,
    }