"""
Quarter-Car Suspension Simulator
--------------------------------
Simulates one wheel's share of a car body (the "quarter car") on a spring and
damper, and studies how damping affects ride comfort.

Model:   m*x'' + c*(x' - y') + k*(x - y) = 0
         x = body displacement, y = road height
         m = sprung mass, k = spring stiffness, c = damper coefficient

Run:     python quarter_car_suspension.py
Output:  five graphs saved in the "figures" folder + results printed in terminal
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------- parameters
M = 250.0                              # sprung mass (kg), one quarter of the body
K = 20000.0                            # spring stiffness (N/m)
C_VALUES = [300, 1500, 4472, 8000]     # damper values to compare (Ns/m)

W0 = np.sqrt(K / M)                    # natural frequency (rad/s)
C_CRIT = 2 * np.sqrt(K * M)            # critical damping (Ns/m)

FIG_DIR = "figures"
os.makedirs(FIG_DIR, exist_ok=True)


def damping_ratio(c):
    return c / C_CRIT


def save(name):
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, name), dpi=150)
    plt.close()


# ---------------------------------------------------------------- simulation
def simulate(c, road, t_end, n_points, x0=0.0, v0=0.0, max_step=np.inf):
    """Solve the quarter-car equation for a given road function road(t) -> (y, dy)."""
    t = np.linspace(0, t_end, n_points)

    def rhs(time, s):
        x, v = s
        y, dy = road(time)
        return [v, (-c * (v - dy) - K * (x - y)) / M]

    sol = solve_ivp(rhs, (0, t_end), [x0, v0], t_eval=t, max_step=max_step)
    x, v = sol.y
    y, dy = np.array([road(ti) for ti in t]).T
    acc = (-c * (v - dy) - K * (x - y)) / M       # body acceleration
    return t, x, v, acc, y


def flat_road(t):
    return 0.0, 0.0


def bump_road(h=0.08, length=2.0, speed=10.0):
    """A smooth half-sine speed bump of height h (m) and given length (m)."""
    tb = length / speed

    def road(t):
        if 0 <= t <= tb:
            return h * np.sin(np.pi * t / tb), h * (np.pi / tb) * np.cos(np.pi * t / tb)
        return 0.0, 0.0

    return road


def wavy_road(amplitude, freq_hz):
    w = 2 * np.pi * freq_hz

    def road(t):
        return amplitude * np.sin(w * t), amplitude * w * np.cos(w * t)

    return road


def transmissibility(freq_hz, c):
    """Theory: body motion / road motion for a sinusoidal road."""
    r = 2 * np.pi * freq_hz / W0
    z = damping_ratio(c)
    return np.sqrt(1 + (2 * z * r) ** 2) / np.sqrt((1 - r**2) ** 2 + (2 * z * r) ** 2)


# ---------------------------------------------------------------- experiments
def free_vibration(c=1500.0, x0=0.10):
    """Check the simulation against the analytical damped-oscillator solution."""
    z = damping_ratio(c)
    wd = W0 * np.sqrt(1 - z**2)
    t, x_num, *_ = simulate(c, flat_road, 3.0, 3000, x0=x0, max_step=0.001)
    x_exact = np.exp(-z * W0 * t) * (x0 * np.cos(wd * t) + z * W0 * x0 / wd * np.sin(wd * t))

    print("--- 1. Free vibration check ---")
    print(f"Natural frequency  = {W0 / (2 * np.pi):.2f} Hz")
    print(f"Damping ratio      = {z:.3f}")
    print(f"Max |numeric - analytical| = {np.max(np.abs(x_num - x_exact)):.2e} m\n")

    plt.plot(t, x_num * 100, label="Numerical")
    plt.plot(t, x_exact * 100, "--", label="Analytical")
    plt.xlabel("Time (s)")
    plt.ylabel("Body displacement (cm)")
    plt.title("Free vibration: simulation vs theory")
    plt.legend()
    plt.grid(True)
    save("1_free_vibration.png")


def damping_comparison(x0=0.10):
    for c in C_VALUES:
        t, x, *_ = simulate(c, flat_road, 3.0, 3000, x0=x0, max_step=0.001)
        plt.plot(t, x * 100, label=f"c = {c}")
    plt.xlabel("Time (s)")
    plt.ylabel("Body displacement (cm)")
    plt.title("Effect of damping (c = 4472 is critical damping)")
    plt.legend()
    plt.grid(True)
    save("2_damping_comparison.png")


def bump_response():
    road = bump_road()
    print("--- 2. Speed bump (8 cm high, 36 km/h) ---")
    for c in C_VALUES:
        t, x, v, acc, y = simulate(c, road, 3.0, 6000, max_step=0.001)
        print(f"c = {c:5d}: peak body height = {x.max() * 100:4.1f} cm, "
              f"peak acceleration = {np.abs(acc).max():5.1f} m/s^2")
        plt.plot(t, x * 100, label=f"c = {c}")
    plt.plot(t, y * 100, "k--", label="Road")
    plt.xlabel("Time (s)")
    plt.ylabel("Height (cm)")
    plt.title("Car body response to a speed bump")
    plt.legend()
    plt.grid(True)
    save("3_speed_bump.png")
    print()


def resonance_sweep(amplitude=0.02):
    freqs = np.linspace(0.3, 6, 60)
    print("--- 3. Resonance sweep (simulation vs theory) ---")
    for c in [300, 1500, 4472]:
        sim = []
        for f in freqs:
            t, x, *_ = simulate(c, wavy_road(amplitude, f), 30.0, 6000)
            sim.append(x[t > 20].max() / amplitude)       # ignore start-up transient
        sim = np.array(sim)
        theory = transmissibility(freqs, c)
        print(f"c = {c:5d}: peak ratio {sim.max():.2f} at {freqs[sim.argmax()]:.2f} Hz "
              f"(theory peak {theory.max():.2f})")
        line, = plt.plot(freqs, sim, label=f"c = {c} (simulation)")
        plt.plot(freqs, theory, "--", color=line.get_color(), label=f"c = {c} (theory)")
    plt.axhline(1, color="gray", linewidth=0.8)
    plt.xlabel("Road frequency (Hz)")
    plt.ylabel("Body motion / road motion")
    plt.title("Resonance of a quarter-car")
    plt.legend()
    plt.grid(True)
    save("4_resonance.png")
    print()


def speed_sweep(wavelength=5.0):
    speeds = np.linspace(5, 80, 200)                   # km/h
    f_hit = (speeds / 3.6) / wavelength                # Hz
    for c in [300, 1500, 4472]:
        plt.plot(speeds, transmissibility(f_hit, c), label=f"c = {c}")
    print("--- 4. Ride vs speed ---")
    print(f"Resonance speed on {wavelength:.0f} m waves: "
          f"{W0 / (2 * np.pi) * wavelength * 3.6:.0f} km/h\n")
    plt.axhline(1, color="gray", linewidth=0.8)
    plt.xlabel("Car speed (km/h)")
    plt.ylabel("Body motion / road motion")
    plt.title(f"Ride vs speed on a wavy road ({wavelength:.0f} m waves)")
    plt.legend()
    plt.grid(True)
    save("5_ride_vs_speed.png")


if __name__ == "__main__":
    free_vibration()
    damping_comparison()
    bump_response()
    resonance_sweep()
    speed_sweep()
    print(f"Done. Graphs saved in the '{FIG_DIR}' folder.")
