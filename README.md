# QS Medi-Cyclone - Quasi-Stationary Mediterranean Cyclones
## Five ways to classify Mediterranean cyclones based on their stationarity.

QS Medi-Cyclone is a pyton script intended to detect the most persistent or stationary cyclone tracks in a dataset. Five metrics are used to define (quasi-) stationarity:

- median speed (FT)
- total distance (FT)
- 12-hours distance (AT)
- radial distance (AT)
- circle distance (AT)

 Each metric is based on either propagation speed or spatial distance and taking into account either the total life cycle of a cyclone (FT, full-track stationarity) or only a part of it (AT, along-track stationarity). By setting specific thresholds, the cyclones can further be assigned to different stationary classes and examined for their properties.
