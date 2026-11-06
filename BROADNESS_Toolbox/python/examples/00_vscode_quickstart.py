# %% 1. Import tools (run this cell first)
"""A small, self-contained BROADNESS lesson for VS Code's Interactive Window.

Run the cells from top to bottom using the `Run Cell` link above each section.
No participant data or MATLAB installation is needed.
"""

import numpy as np
import matplotlib.pyplot as plt

from broadness import BROADNESS_NetworkEstimation

print("BROADNESS is ready. Run cell 2 to make example data.")


# %% 2. Make small example data
# BROADNESS expects sources × time × conditions × participants.
# This synthetic dataset is tiny and is only for learning the workflow.
rng = np.random.default_rng(7)
time = np.linspace(-0.2, 0.6, 80)
source_pattern = rng.normal(size=(12, 1, 1, 1))
signal = np.sin(2 * np.pi * 5 * time)[None, :, None, None]
condition_strength = np.array([1.0, 1.4])[None, None, :, None]
data = source_pattern * signal * condition_strength
data = data + 0.15 * rng.normal(size=(12, 80, 2, 6))

print("Data shape (sources, time, conditions, participants):", data.shape)
print("Time range:", (time[0], time[-1]))


# %% 3. Estimate networks with PCA
pca = BROADNESS_NetworkEstimation(data, time, components_num=4)

print("Network time-series shape:", pca.TimeSeries_BrainNetworks.shape)
print("First four explained-variance percentages:", pca.Variance_BrainNetworks[:4])


# %% 4. Inspect a variable and make a plot
# Public component selections in BROADNESS start at 1, like MATLAB.
# NumPy array positions start at 0, so [:, 0, ...] means component 1.
component_1 = pca.TimeSeries_BrainNetworks[:, 0, :, :]
mean_component_1 = component_1.mean(axis=2)  # time × conditions

fig, ax = plt.subplots()
ax.plot(pca.Time, mean_component_1[:, 0], label="Condition 1")
ax.plot(pca.Time, mean_component_1[:, 1], label="Condition 2")
ax.set(xlabel="Time (s)", ylabel="Component 1 activity", title="First BROADNESS network")
ax.legend()
plt.show()

print("Open the Variables button in the Interactive Window to inspect data, pca, and fig.")
