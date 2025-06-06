# QS MedCyclone - Quasi-Stationary Mediterranean Cyclones
## Five Metrics to Classify Mediterranean Cyclones Based on Their Stationarity

QS MedCyclone contains Python scripts to identify the slowest-moving and/or most transient cyclone tracks in a dataset. It uses five metrics to define cyclone quasi-stationarity: one metric is based on propagation speed, while the other four are based on the distance traveled. These metrics consider either the entire cyclone life cycle (FT, full-track stationarity) or only part of it (AT, along-track stationarity). Together, these diverse metrics provide different perspectives on capturing the quasi-stationary properties of cyclones.

## The Quasi-Stationarity Metrics


- Median Speed (FT)
- Total Distance (FT)
- 12-hour Distance (AT)
- Radial Distance (AT)
- Circle Distance (AT)


![SketchNEW](https://github.com/user-attachments/assets/1039bd13-10c1-4464-8256-491f993829f6)
*Figure 1: Sketch of how the stationarity metrics (based on spatial distance) are calculated.*

## Tutorial
### 1 Load your Data

Use QS_setup.ipynb to upload the cyclone track composite file of Flaounas et al. (2023). In line x and x you can define the temporal and spatial extent of your analysis. The set of cyclone tracks we focus on are displayed as a table. Each track point (row index) of a cyclone (id) represents the minimum pressure at the cyclone centre (hPa) at a specific step in time (year, month, day, time) and space (lon, lat).

![image](https://github.com/user-attachments/assets/f3755185-2042-4e69-9580-8cfe96d092c4)

### 2 Calculate the Quasi-Stationarity Metrics

Use QS_metrics.ipynb to calculate the stationarity metrics for each cyclone.

### 3 Use the QS Table for your Quasi-Stationarity Analysis

QS_example.ipynb provides you possible ways of how the QS-table can be used as a tool to analyse the stationarity properties of the data.
