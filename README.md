# QS Medi-Cyclone - Quasi-Stationary Mediterranean Cyclones
## Five Metrics to Classify Mediterranean Cyclones Based on Their Stationarity.

QS Medi-Cyclone is a pyton script that allows you to detect the most persistent or stationary cyclone tracks in a dataset. Five metrics are presented that can be used to define the (quasi-) stationarity of a cyclone:

- median speed (FT)
- total distance (FT)
- 12-hours distance (AT)
- radial distance (AT)
- circle distance (AT)

 Each metric is based on either propagation speed or spatial distance and taking into account either the total life cycle of a cyclone (FT, full-track stationarity) or only a part of it (AT, along-track stationarity). Through their variety, the five metrics thus allow different perspectives from which the (quasi-) stationarity of a cyclone can be captured. By setting specific thresholds, the cyclones finally can be assigned to different stationary classes and further be examined for their properties.

## Data Origin


## Instructions
### 1 Load your Data

Use ....ipynb to upload your data. In line x and x you can define the temporal and spatial extent of your analysis. The set of cyclone tracks you focus on will finally be displayed as a table. Each track point (row index) is represented as the cyclone's (id) minimum pressure at the center (hPa) at a specific step in time (year, month, day, time) and space (lon, lat).

![image](https://github.com/user-attachments/assets/f3755185-2042-4e69-9580-8cfe96d092c4)

### 2 Calculate the Stationarity Metrics

### 3 
