# Earth–Mars Patched Conics Mission Calculator

A Python tool that sizes a simplified Earth → Mars mission with the **patched conics** method: heliocentric Hohmann transfer, hyperbolic departure from a circular low Earth orbit, and hyperbolic capture into a user-defined elliptical orbit around Mars. The heliocentric leg is also propagated numerically and compared against the analytical solution.

Given three inputs, the program returns the transfer geometry, the hyperbolic excess velocities, the Δv budget, and six plots.

## Inputs

| Input | Description |
|---|---|
| `h_leo` | Altitude of the initial circular parking orbit around Earth [km] |
| `h_p_mars` | Periapsis altitude of the final orbit around Mars [km] |
| `h_a_mars` | Apoapsis altitude of the final orbit around Mars [km] |

## Outputs

- Heliocentric Hohmann transfer: semi-major axis, time of flight, velocities at departure and arrival
- Ideal phase angle of Mars at launch
- Hyperbolic excess velocities $v_\infty$ at Earth and at Mars
- Departure Δv from LEO, capture Δv at Mars, total Δv
- Numerical checks on the propagated heliocentric arc (final radius and velocity errors, path length)
- Six plots (see [Figures](#figures))

## Method

**Model assumptions:** circular, coplanar planetary orbits; sphere-of-influence patching with instantaneous (impulsive) burns applied at periapsis.

1. **Heliocentric transfer.** Hohmann ellipse between the Earth and Mars orbital radii:

$$a_t = \frac{r_E + r_M}{2}, \qquad \text{TOF} = \pi\sqrt{\frac{a_t^3}{\mu_\odot}}$$

   The velocities at the two apsides follow from the vis-viva equation, $v = \sqrt{\mu\left(\frac{2}{r}-\frac{1}{a}\right)}$.

2. **Hyperbolic excess velocities.** Difference between the transfer velocity and the circular velocity of the planet:

$$v_{\infty,E} = \left|v_{t,1} - v_{E}\right|, \qquad v_{\infty,M} = \left|v_{t,2} - v_{M}\right|$$

3. **Departure from LEO.** Tangential burn at the periapsis of the escape hyperbola:

$$\Delta v_{dep} = \sqrt{v_{\infty,E}^2 + \frac{2\mu_E}{r_{LEO}}} - \sqrt{\frac{\mu_E}{r_{LEO}}}$$

4. **Capture at Mars.** Tangential burn at the periapsis of the arrival hyperbola, directly into the target ellipse $(r_p, r_a)$ with $a = (r_p + r_a)/2$:

$$\Delta v_{cap} = \sqrt{v_{\infty,M}^2 + \frac{2\mu_M}{r_p}} - \sqrt{\mu_M\left(\frac{2}{r_p}-\frac{1}{a}\right)}$$

5. **Launch phase angle.** Angle by which Mars must lead Earth at departure, with $n_M$ the mean motion of Mars:

$$\varphi = \pi - n_M\,\text{TOF}$$

6. **Numerical verification.** The planar two-body equations of motion are integrated with `scipy.integrate.solve_ivp` (adaptive Runge–Kutta, `rtol = 1e-10`) from the Earth's position over one transfer time of flight. The final radius and velocity are compared with the Mars orbit radius and the analytical arrival velocity, and the path length is computed with Simpson's rule.

## Installation and usage

```bash
git clone https://github.com/<your-username>/earth-mars-patched-conics.git
cd earth-mars-patched-conics
pip install numpy scipy matplotlib
python earth_mars_patched_conics.py
```

The script asks for the three altitudes interactively. To save the plots as PNG files, set `SAVE_FIGURES = True` at the top of the script.

**Requirements:** Python 3.8+, NumPy, SciPy, Matplotlib.

## Example

Input: `h_leo = 300 km`, Mars orbit of `300 × 20 000 km` altitude.

| Quantity | Value |
|---|---|
| Transfer semi-major axis | 1.262 AU |
| Time of flight | 258.9 days |
| Ideal Mars phase angle at launch | 44.35° |
| $v_\infty$ at Earth | 2.945 km/s |
| $v_\infty$ at Mars | 2.649 km/s |
| Departure Δv (from LEO) | 3.592 km/s |
| Capture Δv (into the target orbit) | 1.020 km/s |
| **Total Δv** | **4.612 km/s** |

These values match the classical Hohmann Earth–Mars figures (about 259 days, $v_\infty \approx$ 2.94 and 2.65 km/s, phase angle about 44°).

## Figures

The script produces:

| | |
|---|---|
| ![Heliocentric transfer](figures/transfer_heliocentric.png) | Heliocentric transfer orbit, with Earth and Mars at launch and Mars at arrival |
| ![Radius vs time](figures/radius_vs_time.png) | Heliocentric radius over time |
| ![Speed vs time](figures/speed_vs_time.png) | Heliocentric speed over time |
| ![Delta-v budget](figures/delta_v_budget.png) | Δv budget: departure, capture, total |
| ![Earth departure](figures/earth_departure.png) | Escape hyperbola from the circular parking orbit |
| ![Mars capture](figures/mars_capture.png) | Arrival hyperbola and insertion into the final elliptical orbit |

## Limitations and possible extensions

- Planetary orbits are circular and coplanar: Mars' eccentricity (about 0.093) and inclination (about 1.85°) are neglected, so no plane-change Δv is included.
- No ephemerides: the phase angle is an ideal value, not a real launch date.
- Impulsive burns only: no gravity or finite-burn losses, no mid-course corrections.
- Hohmann transfer only: no Lambert solver or porkchop plot for real launch windows.
- Mars' gravitational parameter is that of the Mars system (it includes Phobos and Deimos).

## Constants

Distances in km, times in s. $\mu_\odot = 1.32712440041279419 \times 10^{11}$ km³/s², $\mu_E = 398600.435507$ km³/s², $\mu_M = 42828.375816$ km³/s² (Mars system), $R_E = 6371.0084$ km, $R_M = 3389.50$ km, 1 AU = 149 597 870.700 km.

## Author

Andrea Bianchini

## License

MIT (see `LICENSE`).
