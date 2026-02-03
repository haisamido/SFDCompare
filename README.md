# The Merge: Space Flight Dynamics Tools Unification Project

This project aims to unify capabilities from major open-source space flight dynamics tools:
- **OreKit** - Java-based low-level library for space mechanics
- **GMAT** - NASA's General Mission Analysis Tool
- **42** - NASA Goddard's spacecraft attitude control simulation

---

## Feature Comparison: OreKit vs GMAT vs 42

### 1. Propagation

| Feature | OreKit | GMAT | 42 |
|---------|--------|------|-----|
| **Numerical Integrators** | Runge-Kutta (various), Adams-Bashforth, Adams-Moulton, Dormand-Prince | Runge-Kutta, PrinceDormand45/78, Adams-Bashforth-Moulton, Bulirsch-Stoer | 4th-order Runge-Kutta |
| **Analytical Propagators** | Kepler, Brouwer-Lyddane, Eckstein-Hechler | Keplerian (two-body) | Two-body, Three-body |
| **SGP4/SDP4 (TLE)** | ✅ Full support + TLE generation | ✅ TLE propagation (R2022a+) | ❌ |
| **DSST (Semi-analytical)** | ✅ Draper Semi-analytical Satellite Theory | ❌ | ❌ |
| **Ephemeris Propagation** | ✅ SP3, SPICE, tabulated | ✅ SPICE, Code500 ephemeris | ✅ Meeus algorithms |
| **Multi-spacecraft** | ✅ Parallel propagation | ✅ Coupled dynamics, synchronized epochs | ✅ Concurrent multi-spacecraft |
| **CR3BP/Libration Points** | ✅ Halo orbit propagation | ✅ Libration point missions | ✅ Three-body orbits |
| **Flexible Body Dynamics** | ❌ | ❌ | ✅ Rigid and flexible bodies |
| **Multi-body Dynamics** | ❌ | ❌ | ✅ Tree topology joints |

### 2. Force Models

| Feature | OreKit | GMAT | 42 |
|---------|--------|------|-----|
| **Gravity (Point Mass)** | ✅ | ✅ | ✅ |
| **Gravity (Spherical Harmonics)** | ✅ ICGEM, EGM, SHA formats | ✅ COF, GRV, GFC, TAB formats | ✅ EGM96 (Earth), GMM-2B (Mars), GLGM2 (Luna) up to 18x18 |
| **Max Gravity Degree/Order** | Configurable (70x70+) | Configurable (70x70+) | 18x18 |
| **Third Body** | ✅ Sun, Moon, planets | ✅ Sun, Moon, planets | ✅ All planets and major moons |
| **Atmospheric Drag** | ✅ DTM2000, JB2006/2008, NRLMSISE-00, Harris-Priester | ✅ Jacchia-Roberts, MSISE90, JB2008 | ✅ Jacchia-Roberts (Earth), Exponential (Mars) |
| **Solar Radiation Pressure** | ✅ With eclipse modeling | ✅ Basic + N-plate SRP (R2022a+) | ✅ |
| **Solid Tides** | ✅ | ✅ | ❌ |
| **Ocean Tides** | ✅ | ✅ | ❌ |
| **Relativistic Corrections** | ✅ General relativistic effects | ✅ | ❌ |
| **Albedo/IR Radiation** | ✅ Earth albedo and infrared | ❌ | ❌ |
| **Gravity Gradient Torque** | ✅ | ✅ | ✅ |
| **Aerodynamic Torque** | ❌ | ❌ | ✅ |
| **Magnetic Field** | ✅ WMM, IGRF | ❌ | ✅ Planetary magnetic field models |
| **Contact Forces** | ❌ | ❌ | ✅ Spacecraft-surface contact |

### 3. Coordinate Systems & Frames

| Feature | OreKit | GMAT | 42 |
|---------|--------|------|-----|
| **Inertial Frames** | EME2000, GCRF, ICRF, MOD, TOD, TEME | MJ2000Eq, MJ2000Ec, ICRF | J2000, Heliocentric |
| **Earth-Fixed** | ITRF (multiple versions), TIRF | BodyFixed, BodyInertial | Body-fixed for any body |
| **Local Orbital Frames** | LVLH, VNC, TNW, QSW | VNB, LVLH | LVLH, body frames |
| **Body-Centered** | Any celestial body | Earth, Moon, Sun, planets | All solar system bodies |
| **Topocentric** | ✅ Ground station frames | ✅ | ✅ |
| **Barycentric** | ✅ Solar system barycenter | ✅ | ✅ |
| **Libration Point Frames** | ✅ L1-L5 for any system | ✅ | ✅ |
| **User-Defined Frames** | ✅ Hierarchical frame trees | ✅ | ✅ |

### 4. Orbit Representation

| Feature | OreKit | GMAT | 42 |
|---------|--------|------|-----|
| **Cartesian** | ✅ Position/Velocity | ✅ X, Y, Z, VX, VY, VZ | ✅ |
| **Keplerian** | ✅ a, e, i, Ω, ω, ν/M/E | ✅ SMA, ECC, INC, RAAN, AOP, TA/MA | ✅ |
| **Circular** | ✅ For near-circular orbits | ❌ | ❌ |
| **Equinoctial** | ✅ Singularity-free | ✅ ModifiedEquinoctial | ❌ |
| **Spherical** | ❌ | ✅ RMAG, RA, DEC, VMAG, AZI, FPA | ❌ |
| **Two-Line Elements** | ✅ Parse and generate | ✅ Parse and propagate | ❌ |

### 5. Time Systems

| Feature | OreKit | GMAT | 42 |
|---------|--------|------|-----|
| **UTC** | ✅ With leap seconds | ✅ UTCGregorian, UTCModJulian | ✅ |
| **TAI** | ✅ | ✅ TAIGregorian, TAIModJulian | ❌ |
| **TT (Terrestrial Time)** | ✅ | ✅ TTGregorian, TTModJulian | ✅ |
| **TDB (Barycentric)** | ✅ | ✅ TDBGregorian, TDBModJulian | ✅ |
| **GPS Time** | ✅ | ✅ | ✅ |
| **UT1** | ✅ | ❌ | ❌ |
| **Julian Date** | ✅ | ✅ | ✅ |

### 6. Maneuvers

| Feature | OreKit | GMAT | 42 |
|---------|--------|------|-----|
| **Impulsive Burns** | ✅ | ✅ ImpulsiveBurn | ✅ |
| **Finite Burns** | ✅ Continuous thrust | ✅ FiniteBurn with thruster models | ✅ Thruster models |
| **Low-Thrust** | ✅ | ✅ | ✅ |
| **Propulsion Modeling** | ✅ User-defined | ✅ Tanks, Thrusters, ISP, thrust curves | ✅ Thrusters with fuel consumption |
| **Mass Decrement** | ✅ | ✅ | ✅ |
| **Thrust Direction** | ✅ Any frame | ✅ VNB, Body-fixed, inertial | ✅ Body-fixed |
| **Maneuver Triggers** | ✅ Event-based | ✅ Command-based | ✅ Flight software control |

### 7. Solvers & Optimization

| Feature | OreKit | GMAT | 42 |
|---------|--------|------|-----|
| **Differential Corrector** | ❌ (use external) | ✅ Newton-Raphson, Broyden, Modified Broyden | ❌ |
| **Batch Least Squares** | ✅ Levenberg-Marquardt, Gauss-Newton | ✅ Batch Estimator | ❌ |
| **Kalman Filters** | ✅ EKF, UKF, Semi-analytical | ✅ Extended Kalman Filter | ❌ |
| **Smoother** | ✅ RTS smoother | ✅ EKF Smoother (R2022a+) | ❌ |
| **Nonlinear Programming** | ❌ | ✅ VF13ad (SQP), fmincon (MATLAB) | ❌ |
| **Targeting** | ❌ | ✅ Target/Vary/Achieve commands | ❌ |
| **Trajectory Optimization** | ✅ Pontryagin/indirect methods | ✅ Optimize/Minimize commands | ❌ |

### 8. Event Detection

| Feature | OreKit | GMAT | 42 |
|---------|--------|------|-----|
| **Eclipse (Umbra/Penumbra)** | ✅ | ✅ EclipseLocator | ✅ |
| **Ground Station Visibility** | ✅ | ✅ ContactLocator | ✅ |
| **Apogee/Perigee** | ✅ | ✅ Periapsis/Apoapsis stop conditions | ✅ |
| **Node Crossings** | ✅ Ascending/Descending | ✅ | ❌ |
| **Altitude Crossing** | ✅ | ✅ | ❌ |
| **Inter-satellite LOS** | ✅ | ✅ | ✅ |
| **Angular Separation** | ✅ | ✅ | ❌ |
| **Surface Contact** | ❌ | ❌ | ✅ Lander/rover contact |

### 9. Orbit Determination

| Feature | OreKit | GMAT | 42 |
|---------|--------|------|-----|
| **Initial Orbit Determination** | ✅ Gibbs, Herrick-Gibbs, Gooding, Lambert, Gauss, Laplace | ✅ IOD capability (R2022a+) | ❌ |
| **Range Measurements** | ✅ One-way, two-way, TDRSS | ✅ | ❌ |
| **Range-Rate (Doppler)** | ✅ | ✅ | ❌ |
| **Angles (Az/El, RA/Dec)** | ✅ | ✅ | ❌ |
| **GNSS Measurements** | ✅ Code, carrier phase, ambiguity resolution | ❌ Limited | ✅ GPS receiver model |
| **TDOA/FDOA** | ✅ | ❌ | ❌ |
| **Covariance Propagation** | ✅ | ✅ (R2022a+) | ❌ |

### 10. Spacecraft Modeling

| Feature | OreKit | GMAT | 42 |
|---------|--------|------|-----|
| **Mass Properties** | ✅ Dry mass, fuel mass | ✅ DryMass, FuelMass | ✅ Mass, inertia tensor |
| **Drag Properties** | ✅ Cd, drag area | ✅ Cd, DragArea | ✅ Cd, drag area |
| **SRP Properties** | ✅ Cr, SRP area | ✅ Cr, SRPArea | ✅ Cr, SRP area |
| **Tanks** | ✅ Basic | ✅ ChemicalTank, ElectricTank | ✅ |
| **Thrusters** | ✅ Basic | ✅ ChemicalThruster, ElectricThruster | ✅ Multiple thruster types |
| **Power Systems** | ❌ | ✅ SolarPowerSystem, NuclearPowerSystem | ❌ |
| **Flexible Bodies** | ❌ | ❌ | ✅ Modal analysis |
| **Multi-body Joints** | ❌ | ❌ | ✅ Rotational/translational joints |
| **Formations** | ✅ Walker constellations | ✅ Formation object | ✅ Parent-child, peer-to-peer |

### 11. Attitude

| Feature | OreKit | GMAT | 42 |
|---------|--------|------|-----|
| **Attitude Dynamics** | ✅ Kinematic only | ✅ Limited | ✅ Full 6-DOF dynamics |
| **Attitude Laws** | ✅ Nadir, target tracking, yaw compensation, spin, inertial | ✅ CoordinateSystemFixed, Spinner, NadirPointing | ✅ Multiple pointing modes |
| **Euler Angles** | ✅ | ✅ All 12 sequences | ✅ |
| **Quaternions** | ✅ | ✅ | ✅ |
| **Direction Cosine Matrix** | ✅ | ✅ | ✅ |
| **Modified Rodrigues** | ❌ | ✅ | ❌ |
| **Attitude Control Laws** | ❌ | ❌ | ✅ PID, LQR, custom |
| **GNSS-Specific Attitudes** | ✅ GPS, GLONASS, Galileo, Beidou | ❌ | ❌ |

### 12. Sensors & Actuators (GNC Hardware)

| Feature | OreKit | GMAT | 42 |
|---------|--------|------|-----|
| **Gyroscopes** | ❌ | ❌ | ✅ With noise models |
| **Magnetometers** | ❌ | ❌ | ✅ 3-axis |
| **Sun Sensors** | ❌ | ❌ | ✅ Coarse and fine |
| **Star Trackers** | ❌ | ❌ | ✅ With noise models |
| **GPS Receivers** | ✅ Measurement modeling | ❌ | ✅ Position/velocity |
| **Accelerometers** | ❌ | ❌ | ✅ |
| **Reaction Wheels** | ❌ | ❌ | ✅ With momentum management |
| **Magnetic Torquers** | ❌ | ❌ | ✅ |
| **Control Moment Gyros** | ❌ | ❌ | ✅ |
| **Thrusters (ACS)** | ✅ | ✅ | ✅ |

### 13. File Formats

| Feature | OreKit | GMAT | 42 |
|---------|--------|------|-----|
| **TLE** | ✅ Read/Write | ✅ Read | ❌ |
| **SP3** | ✅ versions a-d | ❌ | ❌ |
| **CCSDS OEM** | ✅ | ✅ | ❌ |
| **CCSDS AEM** | ✅ | ✅ | ❌ |
| **SPICE SPK** | ✅ DE4xx, INPOP | ✅ | ❌ |
| **RINEX** | ✅ v2, v3, v4 | ❌ | ❌ |
| **Plain Text Config** | ❌ | ✅ Script files | ✅ Input files |
| **Socket IPC** | ❌ | ❌ | ✅ External app interface |

### 14. Visualization & Output

| Feature | OreKit | GMAT | 42 |
|---------|--------|------|-----|
| **3D Orbit View** | ❌ (external tools) | ✅ OrbitView with OpenGL | ✅ OpenGL visualization |
| **Spacecraft 3D Model** | ❌ | ✅ | ✅ Attitude visualization |
| **Ground Track Plot** | ❌ | ✅ GroundTrackPlot | ❌ |
| **XY Plots** | ❌ | ✅ XYPlot | ❌ |
| **Report Files** | ❌ | ✅ ReportFile | ✅ Text output |
| **Real-time Display** | ❌ | ❌ | ✅ |

### 15. Scripting & Integration

| Feature | OreKit | GMAT | 42 |
|---------|--------|------|-----|
| **Native Language** | Java | C++ | C |
| **Script Interface** | ❌ (API only) | ✅ MATLAB-like script language | ✅ Text input files |
| **GUI** | ❌ | ✅ Full GUI | ✅ OpenGL visualization |
| **Python Bindings** | ✅ via JCC or Orekit-Python wrapper | ✅ via SWIG (experimental) | ❌ |
| **MATLAB Integration** | ❌ | ✅ Native MATLAB function calls | ✅ MATLAB support |
| **Socket IPC** | ❌ | ❌ | ✅ Hardware-in-the-loop |
| **Julia Support** | ❌ | ❌ | ✅ |

### 16. Special Capabilities

| Feature | OreKit | GMAT | 42 |
|---------|--------|------|-----|
| **Collision Probability** | ✅ Multiple methods (Chan, Alfriend, Alfano, Patera) | ❌ | ❌ |
| **GNSS Multi-Constellation** | ✅ GPS, GLONASS, Galileo, Beidou, NavIC, QZSS | ❌ | ❌ |
| **Dilution of Precision** | ✅ GDOP, PDOP, TDOP | ❌ | ❌ |
| **Mission Sequence** | ❌ | ✅ Full mission scripting with control flow | ❌ |
| **Targeting Loops** | ❌ | ✅ Target/Vary/Achieve | ❌ |
| **Proximity Operations** | ❌ | ✅ | ✅ Rendezvous, servicing |
| **Formation Flying** | ✅ | ✅ | ✅ Precision formation |
| **Lander/Rover Ops** | ❌ | ❌ | ✅ Surface contact dynamics |
| **Flight Software Testing** | ❌ | ❌ | ✅ GNC algorithm validation |
| **Hardware-in-the-Loop** | ❌ | ❌ | ✅ Socket IPC |

---

## Summary Comparison

| Aspect | OreKit | GMAT | 42 |
|--------|--------|------|-----|
| **Primary Purpose** | Orbit determination, GNSS, precision ephemeris | Mission design, trajectory optimization | Attitude control system design & test |
| **Architecture** | Low-level library (Java) | Complete application (C++ with GUI/Script) | Simulation framework (C) |
| **Best For** | Library integration, OD, GNSS processing | Mission planning, visualization, analysis | GNC design, ADCS validation, HIL testing |
| **Propagation** | Many analytical options, DSST | Multi-spacecraft synchronization | Multi-body dynamics |
| **Force Models** | Most complete (albedo, tides, etc.) | Good coverage, simpler config | Attitude-relevant forces |
| **Attitude** | Kinematic only | Basic modeling | Full 6-DOF dynamics with sensors/actuators |
| **Sensors/Actuators** | None | None | Comprehensive GNC hardware models |
| **Estimation** | Superior GNSS support, IOD | Integrated targeting | None |
| **Solvers** | External optimization | Built-in DC, NLP, targeting | None |
| **Visualization** | None (external) | Integrated 3D/2D plotting | Real-time 3D attitude display |
| **Learning Curve** | Steeper (API-based) | Moderate (GUI + script) | Moderate (config files) |
| **Use Case** | Flight dynamics, navigation | Mission analysis | ADCS development |

---

## Variable Name Mapping

### Orbital Elements

| Concept | OreKit | GMAT | 42 | SFDaaS |
|---------|--------|------|-----|--------|
| Semi-major axis | `a` | `SMA` | `SMA` | N/A (Cartesian) |
| Eccentricity | `e` | `ECC` | `ecc` | N/A |
| Inclination | `i` | `INC` | `inc` | N/A |
| RAAN | `Ω` (omega) | `RAAN` | `RAAN` | N/A |
| Arg of Periapsis | `ω` (smallOmega) | `AOP` | `ArgP` | N/A |
| True Anomaly | `ν` (nu) | `TA` | `anom` | N/A |
| Mean Anomaly | `M` | `MA` | `MeanAnom` | N/A |
| Position X | `position.getX()` | `X` | `PosN[0]` | `r0[0]` |
| Position Y | `position.getY()` | `Y` | `PosN[1]` | `r0[1]` |
| Position Z | `position.getZ()` | `Z` | `PosN[2]` | `r0[2]` |
| Velocity X | `velocity.getX()` | `VX` | `VelN[0]` | `v0[0]` |
| Velocity Y | `velocity.getY()` | `VY` | `VelN[1]` | `v0[1]` |
| Velocity Z | `velocity.getZ()` | `VZ` | `VelN[2]` | `v0[2]` |

### Propagation Parameters

| Concept | OreKit | GMAT | 42 | SFDaaS |
|---------|--------|------|-----|--------|
| Step size | `stepSize` | `InitialStepSize` | `DT` | `stepSize` |
| Min step | `minStep` | `MinStep` | N/A | `minStep` |
| Max step | `maxStep` | `MaxStep` | N/A | `maxStep` |
| Position tolerance | `positionTolerance` | (part of `Accuracy`) | N/A | `positionTolerance` |
| Velocity tolerance | `velocityTolerance` | (part of `Accuracy`) | N/A | `velocityTolerance` |

### Integrator Types

| Description | OreKit | GMAT | 42 | SFDaaS |
|-------------|--------|------|-----|--------|
| Runge-Kutta 4th order | `ClassicalRungeKuttaIntegrator` | N/A | ✅ (default) | `rungekutta` |
| Runge-Kutta 8/9 | N/A | `RungeKutta89` | N/A | N/A |
| Dormand-Prince 8(5,3) | `DormandPrince853Integrator` | `PrinceDormand78` | N/A | `dormandprince` |
| Adams-Bashforth | `AdamsBashforthIntegrator` | `AdamsBashforthMoulton` | N/A | `adamsbashforth` |
| Adams-Moulton | `AdamsMoultonIntegrator` | (combined above) | N/A | `adamsmoulton` |

### Reference Frames

| Description | OreKit | GMAT | 42 | SFDaaS |
|-------------|--------|------|-----|--------|
| J2000 Equatorial | `FramesFactory.getEME2000()` | `EarthMJ2000Eq` | `J` frame | `eme2000` |
| J2000 Ecliptic | N/A | `EarthMJ2000Ec` | `H` frame | N/A |
| GCRF/ICRF | `FramesFactory.getGCRF()` | `EarthICRF` | N/A | `gcrf` |
| Earth-Fixed | `FramesFactory.getITRF()` | `EarthFixed` | `W` frame | `itrf` |
| Body Frame | N/A | N/A | `B` frame | N/A |
| LVLH | `LOFType.LVLH` | `LVLH` | `L` frame | N/A |

### Spacecraft Properties

| Concept | OreKit | GMAT | 42 | SFDaaS |
|---------|--------|------|-----|--------|
| Dry Mass | `SpacecraftState.getMass()` | `DryMass` | `mass` | N/A |
| Drag Coefficient | `IsotropicDrag.getCd()` | `Cd` | `Cd` | N/A |
| Reflectivity Coeff | `IsotropicRadiationSingleCoefficient.getCr()` | `Cr` | `Cr` | N/A |
| Drag Area | `IsotropicDrag.getCrossSection()` | `DragArea` | `DragArea` | N/A |
| SRP Area | `IsotropicRadiationSingleCoefficient.getCrossSection()` | `SRPArea` | `SRPArea` | N/A |
| Inertia Tensor | N/A | N/A | `I` (3x3) | N/A |

### Attitude Representations

| Concept | OreKit | GMAT | 42 |
|---------|--------|------|-----|
| Quaternion | `Rotation.getQ0/Q1/Q2/Q3()` | `Q1, Q2, Q3, Q4` | `q` (4-vector) |
| Euler Angles | `Rotation.getAngles()` | `EulerAngle1/2/3` | `Ang` (3-vector) |
| DCM | `Rotation.getMatrix()` | `DCM11...DCM33` | `C` (3x3) |
| Angular Velocity | `AngularCoordinates.getRotationRate()` | `EulerAngleRate1/2/3` | `wbn` (3-vector) |

### Sensors (42)

| Sensor Type | 42 Variable | Description |
|-------------|-------------|-------------|
| Gyroscope | `Gyro` | Angular rate sensor |
| Magnetometer | `MAG` | Magnetic field sensor |
| Sun Sensor | `CSS`, `FSS` | Coarse/Fine sun sensors |
| Star Tracker | `ST` | Star tracker |
| GPS Receiver | `GPS` | Position/velocity |
| Accelerometer | `Accel` | Linear acceleration |

### Actuators (42)

| Actuator Type | 42 Variable | Description |
|---------------|-------------|-------------|
| Reaction Wheel | `Whl` | Momentum wheel |
| Magnetic Torquer | `MTB` | Magnetic torque rod |
| Thruster | `Thr` | Thruster |
| Control Moment Gyro | `CMG` | CMG |

---

## Gravitational Parameters (μ in m³/s²)

| Body | OreKit | GMAT | 42 | SFDaaS |
|------|--------|------|-----|--------|
| Earth | `Constants.EGM96_EARTH_MU` (3.986004415e14) | 3.986004418e14 | 3.986004418e14 | 3.986004418e14 |
| Sun | `Constants.JPL_SSD_SUN_GM` | 1.32712440018e20 | 1.32712440018e20 | 1.32712440018e20 |
| Moon | `Constants.JPL_SSD_MOON_GM` | 4.9028e12 | 4.9028e12 | 4.9028e12 |
| Mars | via CelestialBodyFactory | 4.282837e13 | 4.282837e13 | 4.282837e13 |
| Jupiter | via CelestialBodyFactory | 1.26686534e17 | 1.26686534e17 | 1.26686534e17 |

---

## Default Values Comparison

### Propagator Defaults

| Parameter | OreKit | GMAT | 42 | SFDaaS |
|-----------|--------|------|-----|--------|
| Step Size | User-defined | 60 s | User-defined | 60 s |
| Min Step | User-defined | 0.001 s | N/A | 0.001 s |
| Max Step | User-defined | 2700 s | N/A | 1000 s |
| Position Tolerance | User-defined | (via Accuracy) | N/A | 10 m |
| Velocity Tolerance | User-defined | (via Accuracy) | N/A | 0.01 m/s |

### Force Model Defaults

| Property | OreKit | GMAT | 42 |
|----------|--------|------|-----|
| Central Body | User-defined | Earth | Earth |
| Gravity Degree | User-defined | 4 | 18 |
| Gravity Order | User-defined | 4 | 18 |
| Gravity Model | User-defined | JGM2 | EGM96 |
| Atmospheric Model | User-defined | Jacchia-Roberts | Jacchia-Roberts |
| Magnetic Field | IGRF/WMM | N/A | IGRF |

---

## References

- [OreKit 13.1.3 API Documentation](https://www.orekit.org/site-orekit-13.1.3/apidocs/index.html)
- [NASA GMAT Software Catalog](https://software.nasa.gov/software/GSC-19097-1)
- [GMAT Wiki](https://gmat.atlassian.net/wiki/spaces/GW/overview)
- [42 GitHub Repository](https://github.com/ericstoneking/42)
- [42 NASA Software Catalog](https://software.nasa.gov/software/GSC-16720-1)
- [42 Introduction to Simulation (PDF)](https://granasat.ugr.es/wp-content/uploads/2019/02/42-Intro-to-Simulation.pdf)
- [42 NASA Technical Report](https://ntrs.nasa.gov/citations/20180000954)
