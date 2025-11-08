# A fast, robust, physics-based model for the complete radial profile of tropical cyclone wind and pressure

#### Author: Dan Chavas (2025)

https://pypi.org/project/tcwindprofile/ : A package that creates a fast, robust, physics-based radial profile of the tropical cyclone rotating wind and pressure from input Vmax, R34kt, latitude, translation speed, and environmental pressure. Based on the latest observationally-validated science on the structure of the wind field and pressure.

Cite this package: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15442673.svg)](https://doi.org/10.5281/zenodo.15442673) .

Objective: To provide a faster and easier-to-use approximation to the physics-based wind model [code](https://doi.org/10.4231/CZ4P-D448) of Chavas+ (2015 JAS) that is also better fit to observations. That model is state-of-the-art physics, but slow and too rigid to best match observations. Our new model here is fast, analytic, and physics-inspired, but it is also better anchored to observed relationships between R34kt, Rmax, R0, Vmax, and Pmin.  The analytic wind model component is described in: Tao D., Nystrom R., Chavas D. R., and A. Avenas (2025). A fast analytical model for the complete radial structure of tropical cyclone low-level wind field. Geophys. Res. Lett., forthcoming. ([preprint](https://doi.org/10.22541/essoar.175376692.26220135/v1))

This code provides a very good analytic approximation to the wind field model of Chavas+ (2015), which has been extensively validated for real-world TCs in terms of both physics and hazards/impact applications:

Physics:
1. Reproduces characteristic TC wind structure from the entire QuikSCAT and HWIND databases ([Chavas+ 2015 JAS](https://doi.org/10.1175/JAS-D-15-0014.1))
2. Reproduces characteristic modes of TC wind field variability due to variations in intensity and outer size from the Extended Best Track database ([Chavas and Lin 2016 JAS](https://doi.org/10.1175/JAS-D-15-0185.1)).
3. Successfully predicts that wind field structure does not change significantly in a warmer world as seen in both climate-scale and storm-scale models ([Schenkel+ 2023](https://doi.org/10.1175/JCLI-D-22-0066.1))

Hazards/impacts:
1. When used as forcing for a surge model, it reproduces the historical record of U.S. peak storm surge remarkably well ([Gori+ 2023 JGR-A](https://doi.org/10.1029/2022JD037312)). It performs much better than the commonly-used Holland 1980 empirical wind field model [Wang+ 2022 JGR-A](https://doi.org/10.1029/2021JD036359)).
2. When used as forcing for a physics-based rainfall model, it reproduces the climatology of U.S. tropical cyclone inland rainfall remarkably well -- and dramatically better than existing empirical wind field models ([Xi+ 2020 J. Hydromet.](https://doi.org/10.1175/JHM-D-20-0035.1)).
3. When used to model all hazards (wind, coastal surge, inland flooding), predicts the county-level distribution of economic damage quite well ([Gori+ 2025 ERL](https://iopscience.iop.org/article/10.1088/1748-9326/add60d/meta)).

Full modeling pipeline:
1. Estimate Rmax from R34kt: ref [Chavas and Knaff 2022 WAF](https://doi.org/10.1175/WAF-D-21-0103.1) "A Simple Model for Predicting the Tropical Cyclone Radius of Maximum Wind from Outer Size"
2. Estimate R0 from R34kt: analytic approximate solution, from model of ref [Emanuel 2004](https://doi.org/10.1017/CBO9780511735035.010) ("Tropical cyclone energetics and structure") / [Chavas+ 2015 JAS](https://doi.org/10.1175/JAS-D-15-0014.1) "A model for the complete radial structure of the tropical cyclone wind field. Part I: Comparison with observed structure" / [Chavas and Lin 2016 JAS](https://doi.org/10.1175/JAS-D-15-0185.1) "Part II: Wind field variability"
3. Generate wind profile: Analytic complete wind profile: ref Tao et al. (2025, GRL, forthcoming) ([preprint](https://doi.org/10.22541/essoar.175376692.26220135/v1))
    1) eye: r<Rmax (linear model);
    2) inner-core: Rmax to R34kt (linear-M model; Tao+ 2023 GRL);
    3) intermediate radii: R34kt to transition radius (modified Rankine model; Tao+ 2023 GRL, Klotzbach+ 2022 JGRA); and
    4) large radii: transition radius to outer radius (Ekman suction model; Emanuel 2004; Chavas+ 2015/2016 JAS).
6. Estimate Pmin: ref [Chavas Knaff Klotzbach 2025 WAF](https://doi.org/10.1175/WAF-D-24-0031.1) ("A simple model for predicting tropical cyclone minimum central pressure from intensity and size")
7. Generate pressure profile that matches Pmin: same ref as previous

It is very similar to the full physics-based wind profile model of Chavas et al. (2015) ([code here](http://doi.org/10.4231/CZ4P-D448)), but is simpler and much faster, and also includes a more reasonable eye model.

The model starts from the radius of 34kt, which is the most robust measure of size we have: it has long been routinely-estimated operationally; it is at a low enough wind speed to be accurately estimated by satellites over the ocean (higher confidence in data); and it is less noisy because it is typically outside the convective inner-core of the storm. The model then encodes the latest science to estimate 1) Rmax from R34kt (+ Vmax, latitude), 2) the radius of vanishing wind R0 from R34kt (+ latitude, an environmental constant), and 3) the minimum pressure Pmin from Vmax, R34kt, latitude, translation speed, and environmental pressure. Hence, it is very firmly grounded in the known physics of the tropical cyclone wind field while also matching the input data. It is also guaranteed to be very well‐behaved for basically any input parameter combination.

