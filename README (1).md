# Quarter-Car Suspension Simulator

A Python simulation of one wheel's share of a car body (a "quarter car") on a spring and damper. It studies how damping affects ride comfort over a speed bump and over a wavy road, and finds the resonance of the system.

Built as a self-study project after my first-semester Applied Physics labs (damped and forced oscillators, spring constant) at VIT Vellore.

## Model

A single-degree-of-freedom spring-mass-damper driven by the road:

```
m x'' + c (x' - y') + k (x - y) = 0
```

| Symbol | Meaning | Value used |
|---|---|---|
| m | Sprung mass (one quarter of the body) | 250 kg |
| k | Spring stiffness | 20,000 N/m |
| c | Damper coefficient | 300 to 8,000 Ns/m |
| x | Body displacement | |
| y | Road height | |

Natural frequency = 1.42 Hz. Critical damping = 4,472 Ns/m. The base case (c = 1,500) has a damping ratio of 0.335, typical of a passenger car.

## What the code does

1. **Free vibration:** checks the numerical solution (`scipy.integrate.solve_ivp`) against the analytical damped-oscillator formula.
2. **Damping comparison:** under-damped, critically damped and over-damped responses.
3. **Speed bump:** an 8 cm bump at 36 km/h. Compares peak body height and peak acceleration (a measure of discomfort).
4. **Resonance sweep:** simulates a wavy road at many frequencies and compares the result with the theoretical transmissibility curve.
5. **Ride vs speed:** converts road frequency into car speed on a road with 5 m waves, and finds the worst speed.

## Results

| Damper c (Ns/m) | Peak body height on bump (cm) | Peak acceleration (m/s²) |
|---|---|---|
| 300 | 7.7 | 7.3 |
| 1,500 | 6.8 | 12.0 |
| 4,472 (critical) | 7.3 | 22.5 |
| 8,000 | 7.7 | 40.2 |

- The simulation matches theory almost exactly (free vibration error about 1e-15 m; resonance peaks agree to within 0.01).
- Resonance occurs near the natural frequency (about 1.4 Hz). More damping lowers the peak from about 6.9x to 1.8x to 1.15x.
- On 5 m road waves, the worst speed is about 26 km/h.
- Trade-off: a soft damper gives a gentle first impact but keeps bouncing; a stiff damper settles fast but passes harsh jolts to the body. Above about 2 Hz (1.42 x √2) the suspension isolates the body from the road, and there softer damping performs better.

![Resonance](figures/4_resonance.png)
![Speed bump](figures/3_speed_bump.png)

## How to run

```bash
pip install numpy matplotlib scipy
python quarter_car_suspension.py
```

Graphs are saved in the `figures/` folder, and results are printed in the terminal.

## Limitations

- Single-degree-of-freedom model: no tyre stiffness or wheel (unsprung) mass.
- Linear spring and damper; real dampers are non-linear.
- Idealised road (smooth half-sine bump, pure sine waves).

## Possible next steps

- Extend to a two-degree-of-freedom quarter-car model with tyre stiffness.
- Model the spring-damper assembly in SolidWorks and link its dimensions to the simulated stiffness.
- Try a random road profile and compute RMS acceleration.

## Author

Ratnesh, BTech Mechanical Engineering, VIT Vellore.
