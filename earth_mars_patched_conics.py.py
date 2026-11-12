import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, simpson

# ============================================================
# COSTANTI (km, s)
# ============================================================
AU = 149597870.700
MU_SUN = 1.32712440041279419e11
MU_EARTH = 398600.435507
MU_MARS = 42828.375816   # Mars system
R_EARTH = 6371.0084
R_MARS = 3389.50

# Orbite eliocentriche circolari semplificate
R_EARTH_SUN = 1.00000018 * AU
R_MARS_SUN = 1.52371243 * AU

SAVE_FIGURES = False  # metti True se vuoi salvare i grafici in PNG


# ============================================================
# FUNZIONI DI BASE
# ============================================================
def circular_speed(mu, r):
    return np.sqrt(mu / r)


def vis_viva(mu, r, a):
    return np.sqrt(mu * (2.0 / r - 1.0 / a))


def hohmann_transfer(mu_central, r1, r2):
    a_t = 0.5 * (r1 + r2)
    tof = np.pi * np.sqrt(a_t**3 / mu_central)
    v_t1 = vis_viva(mu_central, r1, a_t)
    v_t2 = vis_viva(mu_central, r2, a_t)
    return a_t, tof, v_t1, v_t2


def departure_delta_v(mu_planet, r_parking, v_inf):
    v_circ = circular_speed(mu_planet, r_parking)
    v_per_hyp = np.sqrt(v_inf**2 + 2.0 * mu_planet / r_parking)
    dv = v_per_hyp - v_circ
    return dv, v_circ, v_per_hyp


def capture_delta_v(mu_planet, r_p, a_capture, v_inf):
    v_per_hyp = np.sqrt(v_inf**2 + 2.0 * mu_planet / r_p)
    v_per_ellipse = vis_viva(mu_planet, r_p, a_capture)
    dv = v_per_hyp - v_per_ellipse
    return dv, v_per_hyp, v_per_ellipse


def planar_two_body(t, Y, mu):
    x, y, vx, vy = Y
    r = np.sqrt(x**2 + y**2)

    dxdt = vx
    dydt = vy
    dvxdt = -mu * x / r**3
    dvydt = -mu * y / r**3

    return [dxdt, dydt, dvxdt, dvydt]


def normalize_angle_deg(angle_deg):
    return angle_deg % 360.0


def hyperbola_geometry(mu, r_p, v_inf, n_points=600):
    """
    Restituisce un ramo di iperbole con periapside a theta = 0.
    Utile solo per grafico locale.
    """
    a_abs = mu / v_inf**2
    e = 1.0 + r_p / a_abs
    theta_inf = np.arccos(-1.0 / e)
    theta = np.linspace(-0.92 * theta_inf, 0.92 * theta_inf, n_points)
    p = r_p * (1.0 + e)
    r = p / (1.0 + e * np.cos(theta))
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y, e, theta_inf


def ellipse_geometry(a, e, n_points=800):
    theta = np.linspace(0.0, 2.0 * np.pi, n_points)
    p = a * (1.0 - e**2)
    r = p / (1.0 + e * np.cos(theta))
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y


# ============================================================
# INPUT
# ============================================================
print("=== Earth -> Mars patched-conics mission ===")
h_leo = float(input("Quota orbita circolare iniziale attorno alla Terra [km]: "))
h_p_mars = float(input("Quota del periapside orbita finale attorno a Marte [km]: "))
h_a_mars = float(input("Quota dell'apoapside orbita finale attorno a Marte [km]: "))

if h_leo <= 0:
    raise ValueError("La quota LEO deve essere positiva.")

if h_p_mars <= 0 or h_a_mars <= 0:
    raise ValueError("Le quote finali attorno a Marte devono essere positive.")

if h_a_mars < h_p_mars:
    raise ValueError("L'apoapside finale deve essere >= del periapside finale.")


# ============================================================
# 1) ORBITA INIZIALE TERRESTRE
# ============================================================
r_leo = R_EARTH + h_leo
v_leo = circular_speed(MU_EARTH, r_leo)

# ============================================================
# 2) ORBITA FINALE MARZIANA
# ============================================================
r_p_mars = R_MARS + h_p_mars
r_a_mars = R_MARS + h_a_mars
a_capture = 0.5 * (r_p_mars + r_a_mars)
e_capture = (r_a_mars - r_p_mars) / (r_a_mars + r_p_mars)

# ============================================================
# 3) TRASFERIMENTO ELIOCENTRICO TERRA -> MARTE
# ============================================================
a_transfer, tof, v_t1, v_t2 = hohmann_transfer(MU_SUN, R_EARTH_SUN, R_MARS_SUN)

v_earth_helio = circular_speed(MU_SUN, R_EARTH_SUN)
v_mars_helio = circular_speed(MU_SUN, R_MARS_SUN)

v_inf_earth = abs(v_t1 - v_earth_helio)
v_inf_mars = abs(v_t2 - v_mars_helio)

dv_departure, v_circ_leo, v_per_hyp_earth = departure_delta_v(MU_EARTH, r_leo, v_inf_earth)
dv_capture, v_per_hyp_mars, v_per_ell_mars = capture_delta_v(MU_MARS, r_p_mars, a_capture, v_inf_mars)

dv_total = dv_departure + dv_capture

# angolo di fase ideale di Marte al lancio, rispetto alla Terra
n_mars = np.sqrt(MU_SUN / R_MARS_SUN**3)
phase_angle_rad = np.pi - n_mars * tof
phase_angle_deg = normalize_angle_deg(np.degrees(phase_angle_rad))

# ============================================================
# 4) PROPAGAZIONE NUMERICA DEL TRATTO ELIOCENTRICO
# ============================================================
# Stato iniziale: Terra su asse x, velocità tangenziale +y
Y0 = [R_EARTH_SUN, 0.0, 0.0, v_t1]

t_eval = np.linspace(0.0, tof, 2500)
sol = solve_ivp(
    planar_two_body,
    (0.0, tof),
    Y0,
    t_eval=t_eval,
    args=(MU_SUN,),
    rtol=1e-10,
    atol=1e-10
)

x = sol.y[0]
y = sol.y[1]
vx = sol.y[2]
vy = sol.y[3]
t = sol.t

r_num = np.sqrt(x**2 + y**2)
v_num = np.sqrt(vx**2 + vy**2)

path_length = simpson(v_num, x=t)

r_final_num = r_num[-1]
v_final_num = v_num[-1]

err_r_final = r_final_num - R_MARS_SUN
err_v_initial = v_num[0] - v_t1
err_v_final = v_final_num - v_t2

# ============================================================
# 5) STAMPA RISULTATI
# ============================================================
print("\n" + "=" * 72)
print("RIEPILOGO MISSIONE")
print("=" * 72)

print("\n--- Input missione ---")
print(f"Quota LEO iniziale                  = {h_leo:12.3f} km")
print(f"Quota periapside finale Marte      = {h_p_mars:12.3f} km")
print(f"Quota apoapside finale Marte       = {h_a_mars:12.3f} km")

print("\n--- Orbita iniziale Terra ---")
print(f"r_LEO                              = {r_leo:12.3f} km")
print(f"v_LEO                              = {v_leo:12.6f} km/s")

print("\n--- Orbita finale Marte ---")
print(f"r_p Marte                          = {r_p_mars:12.3f} km")
print(f"r_a Marte                          = {r_a_mars:12.3f} km")
print(f"a cattura                          = {a_capture:12.3f} km")
print(f"e cattura                          = {e_capture:12.6f}")

print("\n--- Trasferimento eliocentrico ---")
print(f"a trasferimento                    = {a_transfer:12.3f} km")
print(f"Tempo di volo                      = {tof:12.3f} s")
print(f"Tempo di volo                      = {tof/86400:12.6f} giorni")
print(f"Vel. Terra eliocentrica            = {v_earth_helio:12.6f} km/s")
print(f"Vel. Marte eliocentrica            = {v_mars_helio:12.6f} km/s")
print(f"Vel. trasferimento a partenza      = {v_t1:12.6f} km/s")
print(f"Vel. trasferimento ad arrivo       = {v_t2:12.6f} km/s")
print(f"v_inf uscita dalla Terra           = {v_inf_earth:12.6f} km/s")
print(f"v_inf arrivo a Marte               = {v_inf_mars:12.6f} km/s")
print(f"Angolo di fase ideale al lancio    = {phase_angle_deg:12.6f} deg")

print("\n--- Delta-v ---")
print(f"Delta-v partenza da LEO            = {dv_departure:12.6f} km/s")
print(f"Delta-v cattura a Marte            = {dv_capture:12.6f} km/s")
print(f"Delta-v totale                     = {dv_total:12.6f} km/s")

print("\n--- Velocità locali ---")
print(f"Vel. periapside iperbole Terra     = {v_per_hyp_earth:12.6f} km/s")
print(f"Vel. periapside iperbole Marte     = {v_per_hyp_mars:12.6f} km/s")
print(f"Vel. periapside orbita finale      = {v_per_ell_mars:12.6f} km/s")

print("\n--- Controlli numerici sul tratto eliocentrico ---")
print(f"Cammino numerico sul transfer      = {path_length:12.3f} km")
print(f"r finale numerico                  = {r_final_num:12.3f} km")
print(f"Errore su r finale                 = {err_r_final:12.6f} km")
print(f"Errore su v iniziale               = {err_v_initial:12.6e} km/s")
print(f"Errore su v finale                 = {err_v_final:12.6e} km/s")

# ============================================================
# 6) GRAFICI
# ============================================================
theta = np.linspace(0.0, 2.0 * np.pi, 800)

# orbite dei pianeti
x_earth_orbit = R_EARTH_SUN * np.cos(theta)
y_earth_orbit = R_EARTH_SUN * np.sin(theta)
x_mars_orbit = R_MARS_SUN * np.cos(theta)
y_mars_orbit = R_MARS_SUN * np.sin(theta)

# fase di Marte al tempo t=0
x_mars0 = R_MARS_SUN * np.cos(np.radians(phase_angle_deg))
y_mars0 = R_MARS_SUN * np.sin(np.radians(phase_angle_deg))

# posizione ideale di Marte all'arrivo
mars_arrival_angle = np.radians(phase_angle_deg) + n_mars * tof
x_mars_arr = R_MARS_SUN * np.cos(mars_arrival_angle)
y_mars_arr = R_MARS_SUN * np.sin(mars_arrival_angle)

# --- Grafico 1: trasferimento eliocentrico
plt.figure(figsize=(8, 8))
plt.plot(x_earth_orbit, y_earth_orbit, label="Orbita Terra")
plt.plot(x_mars_orbit, y_mars_orbit, label="Orbita Marte")
plt.plot(x, y, linewidth=2, label="Transfer orbit (numerico)")
plt.scatter([0], [0], s=80, label="Sole")
plt.scatter([R_EARTH_SUN], [0], s=40, label="Terra al lancio")
plt.scatter([x_mars0], [y_mars0], s=40, label="Marte al lancio")
plt.scatter([x_mars_arr], [y_mars_arr], s=40, label="Marte all'arrivo")
plt.xlabel("x [km]")
plt.ylabel("y [km]")
plt.title("Trasferimento eliocentrico Terra -> Marte")
plt.axis("equal")
plt.grid(True)
plt.legend()

if SAVE_FIGURES:
    plt.savefig("transfer_heliocentric.png", dpi=200, bbox_inches="tight")

# --- Grafico 2: raggio eliocentrico nel tempo
plt.figure(figsize=(9, 4.5))
plt.plot(t / 86400, r_num)
plt.axhline(R_EARTH_SUN, linestyle="--", label="Orbita Terra")
plt.axhline(R_MARS_SUN, linestyle="--", label="Orbita Marte")
plt.xlabel("Tempo [giorni]")
plt.ylabel("r eliocentrico [km]")
plt.title("Raggio eliocentrico nel tempo")
plt.grid(True)
plt.legend()

if SAVE_FIGURES:
    plt.savefig("radius_vs_time.png", dpi=200, bbox_inches="tight")

# --- Grafico 3: velocità eliocentrica nel tempo
plt.figure(figsize=(9, 4.5))
plt.plot(t / 86400, v_num)
plt.axhline(v_t1, linestyle="--", label="v trasferimento in partenza")
plt.axhline(v_t2, linestyle="--", label="v trasferimento in arrivo")
plt.xlabel("Tempo [giorni]")
plt.ylabel("Velocità eliocentrica [km/s]")
plt.title("Velocità eliocentrica nel tempo")
plt.grid(True)
plt.legend()

if SAVE_FIGURES:
    plt.savefig("speed_vs_time.png", dpi=200, bbox_inches="tight")

# --- Grafico 4: budget delta-v
plt.figure(figsize=(7, 4.5))
labels = ["Partenza", "Cattura", "Totale"]
values = [dv_departure, dv_capture, dv_total]
plt.bar(labels, values)
plt.ylabel("Delta-v [km/s]")
plt.title("Budget delta-v missione")
plt.grid(True, axis="y")

if SAVE_FIGURES:
    plt.savefig("delta_v_budget.png", dpi=200, bbox_inches="tight")

# --- Grafico 5: iperbole di partenza dalla Terra
x_hyp_e, y_hyp_e, e_hyp_e, theta_inf_e = hyperbola_geometry(MU_EARTH, r_leo, v_inf_earth)
x_parking = r_leo * np.cos(theta)
y_parking = r_leo * np.sin(theta)

plt.figure(figsize=(7, 7))
plt.plot(x_parking, y_parking, label="LEO circolare")
plt.plot(x_hyp_e, y_hyp_e, label="Iperbole di fuga")
earth_disk_x = R_EARTH * np.cos(theta)
earth_disk_y = R_EARTH * np.sin(theta)
plt.fill(earth_disk_x, earth_disk_y, alpha=0.3, label="Terra")
plt.scatter([r_leo], [0], s=40, label="Punto di iniezione")
plt.xlabel("x [km]")
plt.ylabel("y [km]")
plt.title("Partenza dalla Terra")
plt.axis("equal")
plt.grid(True)
plt.legend()

if SAVE_FIGURES:
    plt.savefig("earth_departure.png", dpi=200, bbox_inches="tight")

# --- Grafico 6: arrivo a Marte e cattura
x_hyp_m, y_hyp_m, e_hyp_m, theta_inf_m = hyperbola_geometry(MU_MARS, r_p_mars, v_inf_mars)
x_cap, y_cap = ellipse_geometry(a_capture, e_capture)
mars_disk_x = R_MARS * np.cos(theta)
mars_disk_y = R_MARS * np.sin(theta)

plt.figure(figsize=(7, 7))
plt.plot(x_hyp_m, y_hyp_m, label="Iperbole di arrivo")
plt.plot(x_cap, y_cap, label="Orbita finale ellittica")
plt.fill(mars_disk_x, mars_disk_y, alpha=0.3, label="Marte")
plt.scatter([r_p_mars], [0], s=40, label="Periapside di cattura")
plt.xlabel("x [km]")
plt.ylabel("y [km]")
plt.title("Arrivo a Marte e inserzione orbitale")
plt.axis("equal")
plt.grid(True)
plt.legend()

if SAVE_FIGURES:
    plt.savefig("mars_capture.png", dpi=200, bbox_inches="tight")

plt.show()