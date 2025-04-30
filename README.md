# QS MediCyclone - Quasi-Stationary Mediterranean Cyclones
## Five Metrics to Classify Mediterranean Cyclones Based on Their Stationarity

QS MediCyclone is a pyton script that allows you to detect the most persistent or stationary cyclone tracks in a dataset. Five metrics are presented that can be used to define the (quasi-) stationarity of a cyclone:

- median speed (FT)
- total distance (FT)
- 12-hours distance (AT)
- radial distance (AT)
- circle distance (AT)

Each metric is based on either propagation speed or spatial distance and taking into account either the total life cycle of a cyclone (FT, full-track stationarity) or only a part of it (AT, along-track stationarity). Through their variety, the five metrics thus allow different perspectives from which the (quasi-) stationarity properties of a cyclone can be captured.

![SketchNEW](https://github.com/user-attachments/assets/1039bd13-10c1-4464-8256-491f993829f6)
*Figure 1: Sketch of how the stationarity metrics (based on spatial distance) are calculated.*

## Tutorial
### 1 Load your Data

Use QS_setup.ipynb to upload your data. In line x and x you can define the temporal and spatial extent of your analysis. The set of cyclone tracks you focus on will finally be displayed as a table. Each track point (row index) of a cyclone (id) represents the minimum pressure at the cyclone center (hPa) at a specific step in time (year, month, day, time) and space (lon, lat).

![image](https://github.com/user-attachments/assets/f3755185-2042-4e69-9580-8cfe96d092c4)

### 2 Calculate the Quasi-Stationarity Metrics

Use QS_metrics.ipynb to calculate the stationarity metrics of your data.

### 3 Use the QS Table for your Quasi-Stationarity Analysis

QS_example.ipynb provides you possible ways of how the QS-table can be used as a tool to analyse the stationarity properties of your data.
