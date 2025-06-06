# QS MedCyclone - Quasi-Stationary Mediterranean Cyclones
## Five Metrics to Classify Mediterranean Cyclones Based on Their Stationarity

QS MedCyclone is a phyton script that allows you to identify the most persistent or stationary cyclone tracks in a dataset. Five metrics are presented that are used to define the (quasi-) stationarity of a cyclone. Each metric is based on either propagation speed or spatial distance and takes into account either the entire life cycle of a cyclone (FT, full-track stationarity) or only a part of it (AT, along-track stationarity). Due to their diversity, the five metrics offer different perspectives from which the (quasi-) stationarity properties of a cyclone can be captured.

## The (Quasi-) Stationarity Metrics


- Median Speed (FT)
- Total Distance (FT)
- 12-hour Distance (AT)
- Radial Distance (AT)
- Circle Distance (AT)


![SketchNEW](https://github.com/user-attachments/assets/1039bd13-10c1-4464-8256-491f993829f6)
*Figure 1: Sketch of how the stationarity metrics (based on spatial distance) are calculated.*

## Tutorial
### 1 Load your Data

Use QS_setup.ipynb to upload your data. In line x and x you can define the temporal and spatial extent of your analysis. The set of cyclone tracks you focus on will finally be displayed as a table. Each track point (row index) of a cyclone (id) represents the minimum pressure at the cyclone centre (hPa) at a specific step in time (year, month, day, time) and space (lon, lat).

![image](https://github.com/user-attachments/assets/f3755185-2042-4e69-9580-8cfe96d092c4)

### 2 Calculate the Quasi-Stationarity Metrics

Use QS_metrics.ipynb to calculate the stationarity metrics for each cyclone.

### 3 Use the QS Table for your Quasi-Stationarity Analysis

QS_example.ipynb provides you possible ways of how the QS-table can be used as a tool to analyse the stationarity properties of the data.
