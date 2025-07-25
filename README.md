# QS MedCyclones - Quasi-Stationary Mediterranean Cyclones
## Five Metrics to Classify Mediterranean Cyclones Based on Their Stationarity

QS MedCyclones contains Python scripts to identify the slowest-moving and/or most transient cyclone tracks in a dataset. It uses five metrics to define cyclone quasi-stationarity: one metric is based on propagation speed, while the other four are based on the distance traveled. These metrics consider either the entire cyclone life cycle (FT, full-track stationarity) or only part of it (AT, along-track stationarity). Together, these diverse metrics provide different perspectives on capturing the quasi-stationary properties of cyclones.

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
The code and results provided in this repository are based on work from my MSc thesis at the University of Bern, supervised by Prof. Dr. Olivia Romppainen-Martius, entitled [_Quasi-Stationary Mediterranean Cyclones_](https://github.com/milliganci/Quasi-Stationary-MediCyclones/blob/main/docs/MSc_Thesis_Unibe_MGanci.pdf).

If you use any part of this project, please cite it in your publication as:

- Ganci, M. B. (2025). Quasi-Stationary MedCyclones. GitHub. [DOI to follow soon].

---

## Tutorial
### 01. `QS_setup.ipynb` — Load Dataset and Perform Filtering

Use this notebook to upload the cyclone track composite dataset from **Flaounas et al. (2023)**.

The dataset consists of cyclone tracks represented as a table. Each row corresponds to a single time step of a cyclone, with the following fields:
- `id`: Cyclone identifier  
- `lon`, `lat`: Spatial coordinates
- `year`, `month`, `day`, `time`: Timestamp  
- `hPa`: Minimum sea level pressure at the cyclone center (in hectopascal units)

![image](https://github.com/user-attachments/assets/f3755185-2042-4e69-9580-8cfe96d092c4)

This script outputs an indexed dataset, `TRACKS_CL5_onlyMedcrossers.csv`, which includes only those cyclone tracks that intersect grid points over the Mediterranean Sea. This heuristic filtering helps to exclude unwanted heat lows.


### 02. `QS_metrics.ipynb` — Calculate Quasi-Stationarity Metrics

This notebook takes `TRACKS_CL5_onlyMedcrossers.csv` as input and computes stationarity metrics across all cyclones.  
The results are saved in `TRACKS_CL5_QS_Medcrossers.csv`, a CSV file containing key metrics for further analysis.


### 03. `QS_example.ipynb` — Explore Quasi-Stationarity of Mediterranean Cyclones

This notebook provides examples of how to use `TRACKS_CL5_QS_Medcrossers.csv` to visualise the quasi-stationary behaviour of Mediterranean cyclones.  
It includes plotting routines and some analysis workflows to support scientific exploration.

