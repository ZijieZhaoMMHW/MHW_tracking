# Temperature Extreme Tracking and Normalizing  
Methods to track and normalize temperature extreme objects. Specifically:  

## Tracking - `hwtrack.m`  
The tracking method to identify and track spatially contiguous extremes is adapted from Sun et al. (2023). The code is originally written by [Di Sun](https://github.com/cindyisok) and is modified here into a parameter-varying function that allows nouniform lon-lat grids.  
```
%Example
tracks = hwtrack(lon_used, lat_used, hw_ts, nums);
```
Here, `lon_used` and `lat_used` are 1D vectors representing the longitudes and latitudes of the input data, respectively. `hw_ts` is a 3D binary field [lon, lat, time], where a value of 1 indicates the presence of a temperature extreme event and 0 indicates its absence. `nums` defines the minimum number of pixels required for a detected object. The output `tracks` contains the spatiotemporal trajectories of tracks.

## Normalizing - `spn.m`  
Analyzing track entities as sequences of moving objects from a composite perspective is substantially more challenging than interpreting fixed-grid outputs, as each event exhibits distinct spatiotemporal dynamics, with evolving shapes, sizes, and durations. To overcome this issue, we develop a normalization method that scales the entities in space and time so that they can be more directly compared. To compare the spatiotemporal evolution of individual tracks, we normalize their shapes to fit within a unit circle (radius = 1) and standardize their durations from t=0 (onset) to t=1 (end). Specifically, we first compute the centroid of each track at every time step and identify the maximum distance (r) from these centroids to their outermost boundaries over the event's entire lifespan. Subsequently, we reposition each objects by translating all grid points so that the centroid aligns with the origin, and scale spatial dimensions uniformly by a factor 1/r to constrain the maximum radius within unity. Pixels covered by land remain blank in this and future composite process.  
```
%Example
tracks_normalized = spn(tracks, property, scale);
```
Where `tracks` is the output from `hwtrack.m`. `property` is the 3D dataset you want to project into the domain of moving tracks, which can be temperature, sensible heat flux, or any other variable. `scale` indicates the size of the normalized domain relative to the maximum extent of tracks. When `scale=1`, it is identical to the approach used in Holbrook et al. (submitted) and Zhao et al. (submitted); when it varies, it corresponds to the Cell Method proposed in Bian et al. (submitted).

Examples of an output  
![HW Morphological](https://github.com/ZijieZhaoMMHW/MHW_tracking/blob/main/mhwvsmhw_dist7.gif)

## Contact  
Zijie Zhao – Code author and maintainer; Framework Developer-Original 

zijiezhaomj@gmail.com; zijie.zhao@utas.edu.au

Department of Earth System Science, University of California Irvine


Ce Bian – Framework developer-Cell 

cebian@ldeo.columbia.edu

Department of Earth and Environmental Sciences, Columbia University  
