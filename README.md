# CUCI Detection Prototype

Official implementation of the bachelor's thesis prototype:

**Translating GPS Data into a Tool to Detect Candidate Unmapped Cycling Infrastructure in OpenStreetMap**

This repository contains an R Shiny prototype for detecting **Candidate Unmapped Cycling Infrastructure (CUCI)** by comparing observed GPS cycling traces with an OpenStreetMap-derived reference network. The tool was developed as a method-testing prototype for a controlled cycling route case study in Enschede, the Netherlands.

The prototype is intended to support visual inspection and prioritisation of possible map-data discrepancies. Detected candidates should not be interpreted as confirmed unmapped infrastructure without further validation using satellite imagery, OSM tag inspection, or fieldwork.

## Repository Structure

.
├── app.R
├── data/
│   ├── cleaned_gps.csv
│   └── overijssel.gpkg
├── outputs/
│   ├── cuci_candidate_summary.csv
│   ├── cuci_candidate_segments.geojson
│   └── analysed_gps_points.csv
└── README.md

Suggested file roles:

- `app.R` - R Shiny application for CUCI detection and visualisation.
- `data/cleaned_gps.csv` - cleaned GPS cycling traces.
- `data/overijssel.gpkg` - OpenStreetMap reference network GeoPackage.
- `outputs/` - exported candidate summaries, candidate geometries, and analysed GPS points.

## Preparing the Data

### GPS Cycling Data

The app expects a cleaned CSV file containing GPS cycling points. The following columns are recommended:

- `participant_id`
- `trip_id`
- `latitude`
- `longitude`
- `timestamp`
- `elevation_m`

The app automatically detects common alternatives for latitude, longitude, timestamp, trip ID, and participant ID columns.

### OpenStreetMap Reference Network

The prototype uses an OSM-derived GeoPackage as the reference network. For the thesis case study, the expected file is:

overijssel.gpkg


You can either upload this file through the Shiny interface or place it in the same folder as `app.R`.

## Running the Prototype

Install the required R packages:

install.packages(c(
  "shiny",
  "readr",
  "dplyr",
  "sf",
  "leaflet",
  "DT",
  "dtw"
))

Run the Shiny app:

shiny::runApp("app.R")


Then:

1. Upload the cleaned GPS CSV.
2. Upload `overijssel.gpkg`, or place it beside `app.R`.
3. Adjust preprocessing and detection settings if needed.
4. Click **Run detection and update map**.
5. Inspect candidate locations on the map and in the summary table.
6. Export candidate summaries and GeoJSON results for validation.

## Detection Settings

The default settings follow the thesis prototype assumptions:

| Parameter | Default | Purpose |
| --- | ---: | --- |
| Inactivity threshold | 900 s | Exclude trips with long time gaps |
| Maximum trip length | 10 km | Exclude trips outside the expected route scale |
| Maximum realistic speed | 50 km/h | Flag likely GPS jump segments |
| Altitude range | -20 to 150 m | Flag unrealistic elevation outliers |
| DTW interpolation spacing | 2 m | Resample raw and matched traces |
| DTW window size | 0 | Use unconstrained DTW by default |
| Local distance threshold | 10 m | Identify GPS-map deviations |
| Consecutive DTW points | 5 | Avoid isolated point errors |
| Minimum candidate gap length | 10 m | Require a continuous candidate segment |
| Repetition grouping distance | 30 m | Group nearby candidate detections |
| Minimum repeated trips | 2 | Retain candidates supported by multiple trips |
| Minimum repeated participants | 2 | Retain candidates supported by multiple participants, where available |
| OSM crop buffer | 500 m | Limit reference network to the GPS trace area |

## Outputs

The app can export:

- `cuci_candidate_summary_<date>.csv` - table of grouped candidate locations.
- `cuci_candidate_segments_<date>.geojson` - spatial candidate segments.
- `analysed_gps_points_<date>.csv` - GPS points with preprocessing and noise-filter attributes.

## Methodological Notes

The prototype uses nearest-network matching to create simplified map-matched traces. This is transparent and suitable for exploratory method testing, but it is not a full route-based map-matching algorithm. In complex areas, such as intersections, parallel paths, bridges, or dense campus networks, nearest-network matching may snap points to nearby features that were not actually used by cyclists.

CUCI outputs should therefore be interpreted as **candidate discrepancies** between observed cycling movement and the mapped reference network. Further validation is required before deciding whether a segment represents unmapped, misclassified, or otherwise incomplete cycling infrastructure in OpenStreetMap.

## Thesis Context

This prototype was developed for a bachelor's thesis in Industrial Design Engineering at the University of Twente. The case study used self-collected GPS cycling traces along a controlled origin-destination route between Enschede city centre and the University of Twente.

## License & Contact

This project is shared for academic and prototype demonstration purposes. For questions, please contact:

**Sabine van der Voorn**

## Citation

If you use or refer to this prototype, please cite:

```bibtex
@thesis{vanderVoorn2026cuci,
  title  = {Translating GPS Data into a Tool to Detect Candidate Unmapped Cycling Infrastructure in OpenStreetMap},
  author = {van der Voorn, Sabine},
  school = {University of Twente},
  year   = {2026},
  type   = {Bachelor thesis}
}
```
