# CUCI Detection Prototype

Official implementation of the bachelor's thesis prototype:

**Translating GPS Data into a Tool to Detect Candidate Unmapped Cycling Infrastructure in OpenStreetMap**

This repository contains an R Shiny prototype for detecting **Candidate Unmapped Cycling Infrastructure (CUCI)** by comparing observed GPS cycling traces with an OpenStreetMap-derived reference network. The tool was developed as a method-testing prototype for a controlled cycling route case study in Enschede, the Netherlands.

The prototype is intended to support visual inspection and prioritisation of possible map-data discrepancies. Detected candidates should not be interpreted as confirmed unmapped infrastructure without further validation using satellite imagery, OpenStreetMap tag inspection, or fieldwork.

## Repository Structure

```text
.
├── app.R
├── data/
│   ├── README.md
│   └── sample_gps.csv
├── outputs/
│   └── .gitkeep
├── .gitignore
├── LICENSE
└── README.md
```

Suggested file roles:

* `app.R` - R Shiny application for CUCI detection and visualisation.
* `data/sample_gps.csv` - small anonymised example GPS cycling dataset.
* `data/README.md` - instructions for preparing input data.
* `outputs/` - folder for exported candidate summaries, candidate geometries, and analysed GPS points.
* `.gitignore` - excludes local outputs, raw/private data, and large OpenStreetMap-derived files.
* `LICENSE` - license for the prototype code.

Large OpenStreetMap-derived files, such as `overijssel.gpkg`, are not included in this repository.

## Preparing the Data

### GPS Cycling Data

The app expects a cleaned CSV file containing GPS cycling points. The following columns are recommended:

* `participant_id`
* `trip_id`
* `latitude`
* `longitude`
* `timestamp`
* `elevation_m`

The app automatically detects common alternatives for latitude, longitude, timestamp, trip ID, and participant ID columns.

Do not upload raw or identifiable participant GPS data to a public repository. For demonstration purposes, only small anonymised sample files should be included.

### OpenStreetMap Reference Network

The prototype uses an OpenStreetMap-derived GeoPackage as the reference network. This file is not included in the repository because it is too large for normal GitHub storage and is derived from OpenStreetMap data.

Download the most recent Overijssel GeoPackage from Geofabrik:

```text
https://download.geofabrik.de/europe/netherlands/overijssel.html
```

Select the GeoPackage download:

```text
overijssel-latest-free.gpkg.zip
```

After downloading:

1. Unzip the file.
2. Place the extracted GeoPackage beside `app.R`, or upload it manually through the Shiny interface.
3. If needed, rename the extracted file to:

```text
overijssel.gpkg
```

Using the most recent Geofabrik download ensures that users can work with the latest available OpenStreetMap data. However, because OpenStreetMap is continuously updated, results may differ slightly depending on the date on which the GeoPackage was downloaded.

For reproducibility, users are encouraged to record the download date of the GeoPackage when reporting or comparing results.

The data is processed by Geofabrik and created by OpenStreetMap contributors. Please respect the OpenStreetMap attribution requirements and the Open Database License (ODbL 1.0).

## Running the Prototype

Install the required R packages:

```r
install.packages(c(
  "shiny",
  "readr",
  "dplyr",
  "sf",
  "leaflet",
  "DT",
  "dtw"
))
```

Run the Shiny app:

```r
shiny::runApp("app.R")
```

Then:

1. Upload the cleaned GPS CSV.
2. Upload the Overijssel GeoPackage, or place it beside `app.R`.
3. Adjust preprocessing and detection settings if needed.
4. Click **Run detection and update map**.
5. Inspect candidate locations on the map and in the summary table.
6. Export candidate summaries and GeoJSON results for validation.

## Detection Settings

The default settings follow the thesis prototype assumptions:

| Parameter                     |      Default | Purpose                                                               |
| ----------------------------- | -----------: | --------------------------------------------------------------------- |
| Inactivity threshold          |        900 s | Exclude trips with long time gaps                                     |
| Maximum trip length           |        10 km | Exclude trips outside the expected route scale                        |
| Maximum realistic speed       |      50 km/h | Flag likely GPS jump segments                                         |
| Altitude range                | -20 to 150 m | Flag unrealistic elevation outliers                                   |
| DTW interpolation spacing     |          2 m | Resample raw and matched traces                                       |
| DTW window size               |            0 | Use unconstrained DTW by default                                      |
| Local distance threshold      |         10 m | Identify GPS-map deviations                                           |
| Consecutive DTW points        |            5 | Avoid isolated point errors                                           |
| Minimum candidate gap length  |         10 m | Require a continuous candidate segment                                |
| Repetition grouping distance  |         30 m | Group nearby candidate detections                                     |
| Minimum repeated trips        |            2 | Retain candidates supported by multiple trips                         |
| Minimum repeated participants |            2 | Retain candidates supported by multiple participants, where available |
| OSM crop buffer               |        500 m | Limit reference network to the GPS trace area                         |

## Outputs

The app can export:

* `cuci_candidate_summary_<date>.csv` - table of grouped candidate locations.
* `cuci_candidate_segments_<date>.geojson` - spatial candidate segments.
* `analysed_gps_points_<date>.csv` - GPS points with preprocessing and noise-filter attributes.

Exported files are written locally and should not be committed to the repository unless they are anonymised and intended as example outputs.

## Methodological Notes

The prototype uses nearest-network matching to create simplified map-matched traces. This is transparent and suitable for exploratory method testing, but it is not a full route-based map-matching algorithm. In complex areas, such as intersections, parallel paths, bridges, or dense campus networks, nearest-network matching may snap points to nearby features that were not actually used by cyclists.

CUCI outputs should therefore be interpreted as **candidate discrepancies** between observed cycling movement and the mapped reference network. Further validation is required before deciding whether a segment represents unmapped, misclassified, or otherwise incomplete cycling infrastructure in OpenStreetMap.

## Thesis Context

This prototype was developed for a bachelor's thesis in Industrial Design Engineering at the University of Twente. The case study used self-collected GPS cycling traces along a controlled origin-destination route between Enschede city centre and the University of Twente.

## License

The R Shiny prototype code is shared under the license included in this repository.

This license applies to the code in this repository. It does not automatically apply to:

* OpenStreetMap-derived data;
* Geofabrik downloads;
* raw or processed participant GPS data;
* thesis text, figures, or university/client-owned materials.

OpenStreetMap-derived data is subject to OpenStreetMap attribution requirements and the Open Database License (ODbL 1.0).

## Contact

For questions, please contact:

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
