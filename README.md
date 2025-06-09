# QS MedCyclone - Quasi-Stationary Mediterranean Cyclones
## Five Metrics to Classify Mediterranean Cyclones Based on Their Stationarity

QS MedCyclone contains Python scripts to identify the slowest-moving and/or most transient cyclone tracks in a dataset. It uses five metrics to define cyclone quasi-stationarity: one metric is based on propagation speed, while the other four are based on the distance traveled. These metrics consider either the entire cyclone life cycle (FT, full-track stationarity) or only part of it (AT, along-track stationarity). Together, these diverse metrics provide different perspectives on capturing the quasi-stationary properties of cyclones.

---

## The Quasi-Stationarity Metrics

- Median Speed (FT)
- Total Distance (FT)
- 12-hour Distance (AT)
- Radial Distance (AT)
- Circle Distance (AT)

![SketchNEW](https://github.com/user-attachments/assets/1039bd13-10c1-4464-8256-491f993829f6)
*Figure 1: Sketch of how the stationarity metrics (based on spatial distance) are calculated.*

---

## More information
...work based on MSc Thesis at University of Bern
...how to cite this repository if use this code in your research: xyz

---

## Tutorial
### 01. `QS_setup.ipynb` — Load Dataset and Perform Filtering

Use this notebook to upload the cyclone track composite dataset from **Flaounas et al. (2023)**.  
You can define the temporal and spatial extent of your analysis in lines `x` and `x`.

The dataset consists of cyclone tracks represented as a table. Each row corresponds to a single time step of a cyclone, with the following fields:
- `id`: Cyclone identifier  
- `lon`, `lat`: Spatial coordinates
- `year`, `month`, `day`, `time`: Timestamp  
- `mslp`: Minimum sea level pressure at the cyclone center (hPa)

![image](https://github.com/user-attachments/assets/f3755185-2042-4e69-9580-8cfe96d092c4)

This script outputs an indexed dataset, `MedCrossers.mat`, which includes only those cyclone tracks that intersect grid points over the Mediterranean Sea. This heuristic filtering helps to exclude unwanted heat lows.


### 02. `QS_metrics.ipynb` — Calculate Quasi-Stationarity Metrics

This notebook takes `MedCrossers.mat` as input and computes stationarity metrics across all cyclones.  
The results are saved in `df_QS.csv`, a CSV file containing key metrics for further analysis.


### 03. `QS_example.ipynb` — Explore Quasi-Stationarity of Mediterranean Cyclones

This notebook provides examples of how to use `df_QS.csv` to analyse the quasi-stationary behavior of Mediterranean cyclones.  
It includes plotting routines and analysis workflows to support scientific exploration.

