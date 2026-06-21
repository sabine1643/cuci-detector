library(shiny)
library(readr)
library(dplyr)
library(sf)
library(leaflet)
library(DT)
library(dtw)

make_lines_from_points <- function(pts_sf, id_col = ".trip_id", crs_out = 4326) {
  if (is.null(pts_sf) || nrow(pts_sf) < 2) return(NULL)
  
  pts_split <- split(pts_sf, pts_sf[[id_col]])
  
  lines <- lapply(names(pts_split), function(id) {
    trip_pts <- pts_split[[id]]
    if (nrow(trip_pts) < 2) return(NULL)
    
    if (".point_index" %in% names(trip_pts)) {
      trip_pts <- trip_pts %>% arrange(.point_index)
    } else if (".timestamp_order" %in% names(trip_pts)) {
      trip_pts <- trip_pts %>% arrange(.timestamp_order)
    }
    
    coords <- st_coordinates(trip_pts)
    line <- st_linestring(coords[, c("X", "Y")])
    
    participant_value <- if ("participant_id" %in% names(trip_pts)) {
      as.character(trip_pts$participant_id[1])
    } else {
      NA_character_
    }
    
    trip_value <- if (".trip_id" %in% names(trip_pts)) {
      as.character(trip_pts$.trip_id[1])
    } else {
      id
    }
    
    st_sf(
      trip_id = trip_value,
      detection_id = id,
      participant_id = participant_value,
      geometry = st_sfc(line, crs = st_crs(pts_sf))
    )
  })
  
  lines <- lines[!sapply(lines, is.null)]
  if (length(lines) == 0) return(NULL)
  
  st_transform(do.call(rbind, lines), crs_out)
}

calculate_trip_lengths_from_points <- function(pts_sf_m, id_col = ".trip_id") {
  if (is.null(pts_sf_m) || nrow(pts_sf_m) < 2) {
    return(data.frame(.trip_id = character(), .trip_length_m = numeric()))
  }
  
  pts_split <- split(pts_sf_m, pts_sf_m[[id_col]])
  
  length_list <- lapply(names(pts_split), function(id) {
    trip_pts <- pts_split[[id]]
    
    if (nrow(trip_pts) < 2) {
      return(data.frame(.trip_id = id, .trip_length_m = 0))
    }
    
    coords <- st_coordinates(trip_pts)[, c("X", "Y"), drop = FALSE]
    
    segment_lengths <- sqrt(
      diff(coords[, "X"])^2 +
        diff(coords[, "Y"])^2
    )
    
    data.frame(
      .trip_id = id,
      .trip_length_m = sum(segment_lengths, na.rm = TRUE)
    )
  })
  
  do.call(rbind, length_list)
}

interpolate_line_points <- function(line_sf, spacing_m = 2) {
  if (is.null(line_sf) || nrow(line_sf) == 0) return(NULL)
  
  line_m <- st_transform(line_sf, 28992)
  line_geom <- st_geometry(line_m)[[1]]
  line_length <- as.numeric(st_length(line_geom))
  
  if (is.na(line_length) || line_length <= 0) return(NULL)
  
  n_points <- max(2, ceiling(line_length / spacing_m))
  fractions <- seq(0, 1, length.out = n_points)
  
  sampled <- st_line_sample(line_geom, sample = fractions)
  sampled_points <- st_cast(sampled, "POINT")
  st_crs(sampled_points) <- 28992
  
  participant_value <- if ("participant_id" %in% names(line_sf)) {
    as.character(line_sf$participant_id[1])
  } else {
    NA_character_
  }
  
  st_sf(
    trip_id = line_sf$trip_id[1],
    participant_id = participant_value,
    point_index_interp = seq_along(sampled_points),
    distance_along_m = fractions * line_length,
    geometry = sampled_points
  )
}

dtw_path_basic <- function(raw_coords, matched_coords, window_size = NULL) {
  n <- nrow(raw_coords)
  m <- nrow(matched_coords)
  if (n < 2 || m < 2) return(NULL)
  
  dx <- outer(raw_coords[, 1], matched_coords[, 1], "-")
  dy <- outer(raw_coords[, 2], matched_coords[, 2], "-")
  dist_matrix <- sqrt(dx^2 + dy^2)
  
  use_window <- !is.null(window_size) && !is.na(window_size) && window_size > 0
  
  if (use_window) {
    safe_window_size <- max(as.integer(window_size), abs(n - m) + 2)
    
    alignment <- tryCatch(
      {
        dtw::dtw(
          dist_matrix,
          keep = TRUE,
          distance.only = FALSE,
          window.type = "sakoechiba",
          window.size = safe_window_size
        )
      },
      error = function(e) {
        dtw::dtw(
          dist_matrix,
          keep = TRUE,
          distance.only = FALSE
        )
      }
    )
  } else {
    alignment <- dtw::dtw(
      dist_matrix,
      keep = TRUE,
      distance.only = FALSE
    )
  }
  
  data.frame(
    raw_index = alignment$index1,
    matched_index = alignment$index2,
    local_distance_m = dist_matrix[cbind(alignment$index1, alignment$index2)]
  )
}

extract_gap_lines_from_dtw <- function(raw_interp, dtw_path, distance_threshold = 10, consecutive_threshold = 5, min_gap_length = 10) {
  if (is.null(raw_interp) || is.null(dtw_path) || nrow(dtw_path) == 0) return(NULL)
  
  dtw_path <- dtw_path %>%
    arrange(raw_index, matched_index) %>%
    mutate(is_gap = local_distance_m >= distance_threshold)
  
  if (!any(dtw_path$is_gap)) return(NULL)
  
  dtw_path <- dtw_path %>%
    mutate(run_id = cumsum(is_gap != lag(is_gap, default = first(is_gap))))
  
  gap_runs <- dtw_path %>%
    group_by(run_id, is_gap) %>%
    summarise(
      first_raw_index = min(raw_index),
      last_raw_index = max(raw_index),
      consecutive_dtw_points = n(),
      max_distance_m = max(local_distance_m, na.rm = TRUE),
      mean_distance_m = mean(local_distance_m, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(
      is_gap == TRUE,
      consecutive_dtw_points >= consecutive_threshold
    )
  
  if (nrow(gap_runs) == 0) return(NULL)
  
  raw_coords <- st_coordinates(raw_interp)
  raw_dist <- raw_interp$distance_along_m
  
  gap_lines <- lapply(seq_len(nrow(gap_runs)), function(k) {
    start_i <- gap_runs$first_raw_index[k]
    end_i <- gap_runs$last_raw_index[k]
    
    start_i <- max(1, min(start_i, nrow(raw_coords)))
    end_i <- max(1, min(end_i, nrow(raw_coords)))
    
    if (end_i <= start_i) return(NULL)
    
    gap_length <- raw_dist[end_i] - raw_dist[start_i]
    
    if (is.na(gap_length) || gap_length < min_gap_length) return(NULL)
    
    coords_gap <- raw_coords[start_i:end_i, c("X", "Y")]
    if (nrow(coords_gap) < 2) return(NULL)
    
    direction_info <- calculate_gap_direction(
      raw_coords = raw_coords,
      start_i = start_i,
      lookahead_points = 5L
    )
    
    participant_value <- if ("participant_id" %in% names(raw_interp)) {
      as.character(raw_interp$participant_id[1])
    } else {
      NA_character_
    }
    
    st_sf(
      trip_id = raw_interp$trip_id[1],
      participant_id = participant_value,
      gap_id = k,
      first_raw_index = start_i,
      last_raw_index = end_i,
      travel_direction = direction_info$travel_direction,
      bearing_deg = direction_info$bearing_deg,
      consecutive_dtw_points = gap_runs$consecutive_dtw_points[k],
      gap_length_m = round(gap_length, 2),
      max_distance_m = round(gap_runs$max_distance_m[k], 2),
      mean_distance_m = round(gap_runs$mean_distance_m[k], 2),
      geometry = st_sfc(st_linestring(coords_gap), crs = 28992)
    )
  })
  
  gap_lines <- gap_lines[!sapply(gap_lines, is.null)]
  if (length(gap_lines) == 0) return(NULL)
  
  do.call(rbind, gap_lines)
}

add_repetition_counts_to_gaps <- function(gaps_sf, repeat_distance_m = 15) {
  if (is.null(gaps_sf) || nrow(gaps_sf) == 0) return(NULL)
  
  gaps_m <- st_transform(gaps_sf, 28992)
  nearby <- st_is_within_distance(gaps_m, gaps_m, dist = repeat_distance_m)
  
  n <- nrow(gaps_m)
  visited <- rep(FALSE, n)
  cluster_id <- rep(NA_integer_, n)
  current_cluster <- 0L
  
  for (i in seq_len(n)) {
    if (visited[i]) next
    
    current_cluster <- current_cluster + 1L
    queue <- i
    visited[i] <- TRUE
    cluster_id[i] <- current_cluster
    
    while (length(queue) > 0) {
      node <- queue[1]
      queue <- queue[-1]
      neighbours <- nearby[[node]]
      
      for (j in neighbours) {
        if (!visited[j]) {
          visited[j] <- TRUE
          cluster_id[j] <- current_cluster
          queue <- c(queue, j)
        }
      }
    }
  }
  
  gaps_m$repeat_cluster_id <- cluster_id
  
  has_participants <- "participant_id" %in% names(gaps_m) &&
    any(!is.na(gaps_m$participant_id) & gaps_m$participant_id != "")
  
  if (has_participants) {
    repeat_summary <- gaps_m %>%
      st_drop_geometry() %>%
      group_by(repeat_cluster_id) %>%
      summarise(
        repeated_trip_count = n_distinct(trip_id),
        repeated_participant_count = n_distinct(participant_id[!is.na(participant_id) & participant_id != ""]),
        repeated_gap_count = n(),
        repeated_trip_ids = paste(sort(unique(trip_id)), collapse = ", "),
        repeated_participant_ids = paste(sort(unique(participant_id[!is.na(participant_id) & participant_id != ""])), collapse = ", "),
        .groups = "drop"
      )
  } else {
    repeat_summary <- gaps_m %>%
      st_drop_geometry() %>%
      group_by(repeat_cluster_id) %>%
      summarise(
        repeated_trip_count = n_distinct(trip_id),
        repeated_participant_count = NA_integer_,
        repeated_gap_count = n(),
        repeated_trip_ids = paste(sort(unique(trip_id)), collapse = ", "),
        repeated_participant_ids = NA_character_,
        .groups = "drop"
      )
  }
  
  gaps_m %>%
    left_join(repeat_summary, by = "repeat_cluster_id") %>%
    st_transform(st_crs(gaps_sf))
}

info_icon <- function(text) {
  tags$span(
    class = "info-icon",
    title = text,
    HTML("i")
  )
}

bearing_to_direction <- function(bearing_deg) {
  if (is.na(bearing_deg)) return(NA_character_)
  directions <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
  index <- floor(((bearing_deg + 22.5) %% 360) / 45) + 1
  directions[index]
}

calculate_bearing_deg <- function(start_xy, end_xy) {
  dx <- end_xy["X"] - start_xy["X"]
  dy <- end_xy["Y"] - start_xy["Y"]
  if (is.na(dx) || is.na(dy) || (dx == 0 && dy == 0)) return(NA_real_)
  bearing <- (atan2(dx, dy) * 180 / pi + 360) %% 360
  round(as.numeric(bearing), 1)
}

calculate_gap_direction <- function(raw_coords, start_i, lookahead_points = 5L) {
  if (is.null(raw_coords) || nrow(raw_coords) < 2) {
    return(list(bearing_deg = NA_real_, travel_direction = NA_character_))
  }
  start_i <- max(1L, min(as.integer(start_i), nrow(raw_coords)))
  end_i <- min(start_i + as.integer(lookahead_points), nrow(raw_coords))
  if (end_i == start_i && start_i > 1L) {
    end_i <- start_i
    start_i <- max(1L, start_i - as.integer(lookahead_points))
  }
  bearing <- calculate_bearing_deg(
    start_xy = raw_coords[start_i, c("X", "Y")],
    end_xy = raw_coords[end_i, c("X", "Y")]
  )
  list(
    bearing_deg = bearing,
    travel_direction = bearing_to_direction(bearing)
  )
}

ui <- fluidPage(
  tags$head(
    uiOutput("accessibility_font_css"),
    tags$style(HTML("
      body { 
        background-color: #f7f7f7;
        font-size: var(--app-font-size, 14px);
      }

      h2 {
        font-size: calc(var(--app-font-size, 14px) * 1.8);
      }

      h3 {
        font-size: calc(var(--app-font-size, 14px) * 1.45);
      }

      label,
      .control-label,
      .checkbox,
      .radio,
      .help-block,
      .small-note,
      .app-subtitle,
      .status-box,
      .heuristic-note,
      .btn,
      input,
      select,
      textarea {
        font-size: var(--app-font-size, 14px) !important;
      }

      details.control-card summary {
        font-size: calc(var(--app-font-size, 14px) * 1.15) !important;
      }

      .form-control {
        height: auto;
        min-height: calc(var(--app-font-size, 14px) * 2.4);
        padding: 8px 10px;
      }

      .leaflet-popup-content,
      .leaflet-control,
      .leaflet-control-layers,
      .legend {
        font-size: var(--app-font-size, 14px) !important;
        line-height: 1.45 !important;
      }

      .leaflet-popup-content {
        min-width: 220px;
      }

      .dataTables_wrapper,
      .dataTables_wrapper table,
      .dataTables_wrapper .dataTables_info,
      .dataTables_wrapper .dataTables_paginate,
      .dataTables_wrapper .dataTables_length,
      .dataTables_wrapper .dataTables_filter,
      table.dataTable {
        font-size: var(--app-font-size, 14px) !important;
      }

      table.dataTable td,
      table.dataTable th {
        padding: 8px 10px !important;
      }

      .app-title { margin-bottom: 4px; }
      .app-subtitle { color: #555; margin-bottom: 18px; }
      .control-card {
        background: white;
        border: 1px solid #ddd;
        border-radius: 8px;
        margin-bottom: 12px;
        box-shadow: 0 1px 2px rgba(0,0,0,0.04);
        overflow: hidden;
      }
      details.control-card summary {
        cursor: pointer;
        list-style: none;
        padding: 13px 14px;
        font-weight: 700;
        background: #ffffff;
        border-bottom: 1px solid #eeeeee;
      }
      details.control-card summary::-webkit-details-marker { display: none; }
      details.control-card summary::after {
        content: '\\25B8';
        float: right;
        color: #337ab7;
        transition: transform 0.15s ease-in-out;
      }
      details.control-card[open] summary::after {
        transform: rotate(90deg);
      }
      .section-body { padding: 14px; }
      .small-note { color: #666; }
      .primary-button { width: 100%; font-weight: 700; }
      .download-row .btn { width: 100%; margin-bottom: 8px; }
      .status-box {
        background: #fbfbfb;
        border-left: 4px solid #337ab7;
        padding: 10px;
        white-space: pre-wrap;
      }
      .heuristic-note {
        background: #eef7ff;
        border-left: 4px solid #5bc0de;
        padding: 10px;
        margin-bottom: 12px;
      }

      .info-icon {
        display: inline-block;
        width: calc(var(--app-font-size, 14px) * 1.35);
        height: calc(var(--app-font-size, 14px) * 1.35);
        line-height: calc(var(--app-font-size, 14px) * 1.35);
        text-align: center;
        border-radius: 50%;
        background: #337ab7;
        color: white;
        font-size: calc(var(--app-font-size, 14px) * 0.85);
        font-weight: 700;
        margin-left: 6px;
        cursor: help;
      }
      .input-label-with-info {
        font-weight: 700;
        margin-bottom: 4px;
      }
    "))
  ),
  
  div(class = "app-title", h2("Candidate Unmapped Cycling Infrastructure Detection Prototype")),
  div(
    class = "app-subtitle",
    "Compare GPS cycling traces with the Open Street Map (OSM) reference network and export candidate unmapped cycling infrastructure (CUCI) results."
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      
      div(
        class = "heuristic-note",
        strong("Workflow status: "),
        textOutput("workflow_status", inline = TRUE)
      ),
      
      tags$details(
        class = "control-card",
        tags$summary("Accessibility settings"),
        div(
          class = "section-body",
          selectInput(
            "font_size_choice",
            "Interface font size",
            choices = c(
              "Normal" = "14",
              "Large" = "16",
              "Extra large" = "18",
              "Very large" = "20"
            ),
            selected = "14"
          ),
          div(
            class = "small-note",
            "Increasing the font size enlarges interface text, table text, map popups, legends, and layer controls."
          )
        )
      ),
      
      tags$details(
        class = "control-card",
        open = NA,
        tags$summary("1. Upload input data"),
        div(
          class = "section-body",
          helpText("Upload the cleaned GPS CSV. The app will automatically detect common latitude, longitude, timestamp, trip ID and participant ID column names."),
          fileInput("file", "GPS CSV", accept = ".csv"),
          fileInput("gpkg_file", "OSM network GeoPackage: overijssel.gpkg", accept = c(".gpkg")),
          helpText("If overijssel.gpkg is saved in the same folder as app.R, uploading it here is optional.")
        )
      ),
      
      tags$details(
        class = "control-card",
        tags$summary("2. Preprocessing filters"),
        div(
          class = "section-body",
          helpText("These filters remove whole trips before CUCI detection. They are useful for excluding traces that may distort the analysis."),
          
          checkboxInput(
            "apply_inactivity_filter",
            "Exclude trips with long inactivity/time gaps",
            value = TRUE
          ),
          numericInput(
            "max_time_gap_s",
            "Maximum allowed time gap between consecutive GPS points, in seconds",
            value = 900,
            min = 1,
            step = 30
          ),
          div(class = "small-note", "Only applied when a timestamp column is available."),
          br(),
          
          checkboxInput(
            "apply_max_trip_length_filter",
            "Exclude trips longer than a maximum route length",
            value = TRUE
          ),
          numericInput(
            "max_trip_length_km",
            "Maximum allowed trip length, in kilometres",
            value = 10,
            min = 0.1,
            step = 0.5
          ),
          div(class = "small-note", "Trip length is calculated from consecutive GPS points after ordering."),
          br(),
          
          checkboxInput(
            "apply_speed_noise_filter",
            "Ignore GPS segments with unrealistic speeds during CUCI detection",
            value = TRUE
          ),
          numericInput(
            "max_realistic_speed_kmh",
            "Maximum realistic cycling speed, in km/h",
            value = 50,
            min = 5,
            step = 5
          ),
          div(class = "small-note", "This does not remove the whole trip. It prevents GPS jump segments from becoming CUCI candidates. Only applied when timestamps are available."),
          br(),
          
          checkboxInput(
            "apply_altitude_noise_filter",
            "Ignore GPS points with unrealistic altitudes during CUCI detection",
            value = TRUE
          ),
          numericInput(
            "min_realistic_altitude_m",
            "Minimum realistic altitude, in metres",
            value = -20,
            step = 5
          ),
          numericInput(
            "max_realistic_altitude_m",
            "Maximum realistic altitude, in metres",
            value = 150,
            step = 5
          ),
          div(class = "small-note", "This does not remove the whole trip. It prevents altitude outliers from becoming CUCI candidates. Only applied when an altitude/elevation column is available.")
        )
      ),
      
      tags$details(
        class = "control-card",
        tags$summary("3. Choose visible map layers"),
        div(
          class = "section-body",
          checkboxInput("show_points", "GPS points", TRUE),
          checkboxInput("show_raw_trace", "Raw GPS traces", TRUE),
          checkboxInput("show_map_matched", "Matched traces", TRUE),
          checkboxInput("show_gap_segments", "CUCI candidate gaps", TRUE),
          checkboxInput("show_osm_network", "OSM network layer", TRUE),
          div(class = "small-note", "Layer names avoid overclaiming: orange segments are candidates, not confirmed infrastructure.")
        )
      ),
      
      tags$details(
        class = "control-card",
        tags$summary("4. Detection settings"),
        div(
          class = "section-body",
          helpText("Default values follow the prototype assumptions. Hover over the information icons for an explanation of each setting."),
          div(
            class = "input-label-with-info",
            "Interpolation spacing for DTW, in metres",
            info_icon("Distance between sampled points along the raw and matched traces before DTW comparison. Smaller values are more detailed but slower; larger values are faster but less precise.")
          ),
          numericInput("dtw_spacing", label = NULL, value = 2, min = 1, max = 20),
          div(
            class = "input-label-with-info",
            "DTW window size, in interpolated points",
            info_icon("Limits how far the DTW alignment may warp away from the diagonal. Use 0 for unconstrained DTW. Smaller windows are faster, but too small a window can miss strong deviations.")
          ),
          numericInput("dtw_window_size", label = NULL, value = 0, min = 0, step = 5),
          div(
            class = "input-label-with-info",
            "Local distance threshold, in metres",
            info_icon("Minimum local mismatch between the raw GPS trace and the map-matched trace before a point pair is marked as a deviation.")
          ),
          numericInput("distance_threshold", label = NULL, value = 10, min = 1),
          div(
            class = "input-label-with-info",
            "Minimum consecutive DTW points",
            info_icon("Minimum number of consecutive deviating DTW point pairs needed before the deviation can become a CUCI candidate. This reduces isolated GPS-noise detections.")
          ),
          numericInput("consecutive_threshold", label = NULL, value = 5, min = 1),
          div(
            class = "input-label-with-info",
            "Minimum candidate gap length, in metres",
            info_icon("Minimum physical length of a deviation segment. With the default 2 m interpolation and 5 consecutive DTW points, 10 m matches the minimum continuous candidate signal.")
          ),
          numericInput("min_gap_length", label = NULL, value = 10, min = 1),
          div(
            class = "input-label-with-info",
            "Minimum repeated detections, distinct trips",
            info_icon("Minimum number of different trips that must support the same spatial candidate cluster before it is shown in the final results.")
          ),
          numericInput("min_repeated_trips", label = NULL, value = 2, min = 1, step = 1),
          div(
            class = "input-label-with-info",
            "Minimum repeated detections, distinct participants",
            info_icon("Minimum number of distinct participants that must support a candidate. This is only applied when a participant ID column is available.")
          ),
          numericInput("min_repeated_participants", label = NULL, value = 2, min = 1, step = 1),
          div(
            class = "input-label-with-info",
            "Distance for grouping repeated candidates, in metres",
            info_icon("Candidate gaps within this distance are grouped into the same repeat cluster. Larger values merge nearby detections more easily; smaller values keep them separate.")
          ),
          numericInput("repeat_distance_m", label = NULL, value = 30, min = 1, step = 1),
          div(
            class = "input-label-with-info",
            "OSM crop buffer around GPS traces, in metres",
            info_icon("Only OSM network features within this buffer around the uploaded GPS traces are used. A larger buffer is safer but slower; a smaller buffer is faster but may remove relevant network segments.")
          ),
          numericInput("osm_crop_buffer_m", label = NULL, value = 500, min = 50, step = 50),
          div(class = "small-note", "The participant filter is only applied when a participant ID column is available. The OSM crop buffer limits map-matching to the area around the uploaded GPS traces.")
        )
      ),
      
      tags$details(
        class = "control-card",
        open = NA,
        tags$summary("5. Run and export"),
        div(
          class = "section-body",
          actionButton("draw", "Run detection and update map", class = "btn-primary primary-button"),
          br(), br(),
          div(class = "download-row",
              downloadButton("download_gap_csv", "Export candidate summary as CSV"),
              downloadButton("download_gap_geojson", "Export candidate segments as GeoJSON"),
              downloadButton("download_cleaned_gps_csv", "Export analysed GPS points as CSV")
          ),
          div(class = "small-note", "Exports use the current trip exclusions, preprocessing filters, and detection settings.")
        )
      ),
      
      tags$details(
        class = "control-card",
        tags$summary("Status and feedback"),
        div(
          class = "section-body",
          div(class = "status-box", verbatimTextOutput("status"))
        )
      )
    ),
    
    mainPanel(
      width = 8,
      leafletOutput("map", height = 700),
      br(),
      h3("Detected CUCI candidate locations"),
      helpText("This table summarises candidate locations where cyclists deviated from the map-matched OSM network. Higher evidence candidates are more suitable for visual or field validation."),
      DTOutput("gap_table")
    )
  )
)

server <- function(input, output, session) {
  
  output$accessibility_font_css <- renderUI({
    font_size <- input$font_size_choice
    
    if (is.null(font_size) || is.na(font_size)) {
      font_size <- "14"
    }
    
    tags$style(HTML(paste0(
      ":root { --app-font-size: ", font_size, "px; }"
    )))
  })
  
  status_text <- reactiveVal("Step 1: upload a GPS CSV and overijssel.gpkg, then run detection.")
  
  output$workflow_status <- renderText({
    if (is.null(input$file)) {
      return("waiting for GPS CSV")
    }
    if (is.null(input$draw) || input$draw == 0) {
      return("data uploaded, detection not yet run")
    }
    candidates <- candidate_summary_table()
    paste("detection complete;", nrow(candidates), "candidate location(s) after filtering")
  })
  
  data_in <- reactive({
    req(input$file)
    
    notification_id <- showNotification(
      "Loading GPS CSV and reading trip information...",
      type = "message",
      duration = NULL,
      closeButton = FALSE
    )
    on.exit(removeNotification(notification_id), add = TRUE)
    
    status_text("Loading GPS CSV and detecting available columns...")
    
    df <- read_csv(input$file$datapath, show_col_types = FALSE)
    
    validate(
      need(nrow(df) > 0, "The uploaded CSV appears to be empty.")
    )
    
    status_text(
      paste0(
        "GPS CSV loaded successfully.",
        "\nRows detected: ", nrow(df),
        "\nColumns detected: ", ncol(df),
        "\nNext: adjust detection settings, then run detection."
      )
    )
    
    df
  })
  
  lon_col <- reactive({
    df <- data_in()
    possible_lon_names <- c("longitude", "lon", "long", "lng", "x", "X", "Longitude", "LONGITUDE", "Lon", "LON", "Long", "LONG", "Lng", "LNG")
    match <- names(df)[names(df) %in% possible_lon_names]
    if (length(match) == 0) stop("No longitude column found.")
    match[1]
  })
  
  lat_col <- reactive({
    df <- data_in()
    possible_lat_names <- c("latitude", "lat", "y", "Y", "Latitude", "LATITUDE", "Lat", "LAT")
    match <- names(df)[names(df) %in% possible_lat_names]
    if (length(match) == 0) stop("No latitude column found.")
    match[1]
  })
  
  time_col <- reactive({
    df <- data_in()
    possible_time_names <- c("timestamp", "time", "datetime", "date_time", "Timestamp", "TIME", "DateTime", "DATETIME")
    match <- names(df)[names(df) %in% possible_time_names]
    if (length(match) == 0) return(NULL)
    match[1]
  })
  
  altitude_col <- reactive({
    df <- data_in()
    possible_altitude_names <- c(
      "altitude", "alt", "elevation", "height",
      "altitude_m", "elevation_m", "Altitude", "ALTITUDE",
      "Elevation", "ELEVATION", "Alt", "ALT"
    )
    match <- names(df)[names(df) %in% possible_altitude_names]
    if (length(match) == 0) return(NULL)
    match[1]
  })
  
  trip_col <- reactive({
    df <- data_in()
    possible_trip_names <- c("unique_trip_id", "source_file", "trip_id", "trip", "route_id", "track_id", "Trip_ID", "TRIP_ID", "TripId", "tripid")
    match <- names(df)[names(df) %in% possible_trip_names]
    if (length(match) == 0) return(NULL)
    match[1]
  })
  
  participant_col <- reactive({
    df <- data_in()
    possible_participant_names <- c(
      "participant_id", "participant", "user_id", "user", "rider_id", "rider",
      "Participant_ID", "PARTICIPANT_ID", "ParticipantId", "participantid",
      "User_ID", "USER_ID", "UserId", "userid"
    )
    match <- names(df)[names(df) %in% possible_participant_names]
    if (length(match) == 0) return(NULL)
    match[1]
  })
  
  points_sf <- eventReactive(input$draw, {
    req(data_in())
    
    df <- data_in()
    lon <- lon_col()
    lat <- lat_col()
    time <- time_col()
    trip <- trip_col()
    participant <- participant_col()
    altitude <- altitude_col()
    
    original_point_count <- nrow(df)
    
    df <- df %>%
      filter(!is.na(.data[[lon]]), !is.na(.data[[lat]]))
    
    coordinate_filtered_point_count <- nrow(df)
    
    if (!is.null(trip)) {
      df <- df %>% mutate(.trip_id = as.character(.data[[trip]]))
    } else {
      df <- df %>% mutate(.trip_id = "trip_1")
    }
    
    if (!is.null(participant)) {
      df <- df %>% mutate(participant_id = as.character(.data[[participant]]))
    } else if (!"participant_id" %in% names(df)) {
      df <- df %>% mutate(participant_id = NA_character_)
    } else {
      df <- df %>% mutate(participant_id = as.character(participant_id))
    }
    
    original_trip_ids <- sort(unique(df$.trip_id))
    original_trip_count <- length(original_trip_ids)
    
    validate(need(nrow(df) > 1, "After missing-coordinate filtering, fewer than 2 GPS points remain."))
    
    if (!is.null(time)) {
      df <- df %>%
        mutate(.timestamp_order = as.POSIXct(.data[[time]], tz = "UTC")) %>%
        arrange(.trip_id, .timestamp_order)
      
      if (isTRUE(input$apply_inactivity_filter)) {
        df <- df %>%
          group_by(.trip_id) %>%
          mutate(
            .time_gap_s = as.numeric(difftime(.timestamp_order, lag(.timestamp_order), units = "secs")),
            .has_long_inactivity = any(.time_gap_s > input$max_time_gap_s, na.rm = TRUE)
          ) %>%
          ungroup()
        
        inactivity_removed_trips <- df %>%
          filter(.has_long_inactivity == TRUE) %>%
          distinct(.trip_id) %>%
          pull(.trip_id)
        
        df <- df %>%
          filter(.has_long_inactivity == FALSE)
      } else {
        inactivity_removed_trips <- character()
        df <- df %>%
          mutate(
            .time_gap_s = NA_real_,
            .has_long_inactivity = FALSE
          )
      }
    } else {
      df <- df %>%
        mutate(.row_order = row_number()) %>%
        arrange(.trip_id, .row_order)
      
      inactivity_removed_trips <- character()
    }
    
    validate(need(nrow(df) > 1, "After inactivity filtering, fewer than 2 GPS points remain."))
    
    pts_for_length <- st_as_sf(df, coords = c(lon, lat), crs = 4326, remove = FALSE) %>%
      st_transform(28992)
    
    trip_lengths <- calculate_trip_lengths_from_points(pts_for_length, id_col = ".trip_id")
    
    df <- df %>%
      left_join(trip_lengths, by = ".trip_id")
    
    if (isTRUE(input$apply_max_trip_length_filter)) {
      max_trip_length_m <- input$max_trip_length_km * 1000
      
      length_removed_trips <- df %>%
        filter(.trip_length_m > max_trip_length_m) %>%
        distinct(.trip_id) %>%
        pull(.trip_id)
      
      df <- df %>%
        filter(.trip_length_m <= max_trip_length_m)
    } else {
      length_removed_trips <- character()
    }
    
    validate(need(nrow(df) > 1, "After trip length filtering, fewer than 2 GPS points remain."))
    
    df_for_noise <- st_as_sf(df, coords = c(lon, lat), crs = 4326, remove = FALSE) %>%
      st_transform(28992)
    
    noise_coords <- st_coordinates(df_for_noise)[, c("X", "Y"), drop = FALSE]
    df$.x_m_for_noise <- noise_coords[, "X"]
    df$.y_m_for_noise <- noise_coords[, "Y"]
    
    df <- df %>%
      group_by(.trip_id) %>%
      mutate(
        .segment_distance_m = sqrt(
          (.x_m_for_noise - lag(.x_m_for_noise))^2 +
            (.y_m_for_noise - lag(.y_m_for_noise))^2
        ),
        .segment_time_s = if (".timestamp_order" %in% names(.)) {
          as.numeric(difftime(.timestamp_order, lag(.timestamp_order), units = "secs"))
        } else {
          NA_real_
        },
        .segment_speed_kmh = ifelse(
          !is.na(.segment_time_s) & .segment_time_s > 0,
          (.segment_distance_m / .segment_time_s) * 3.6,
          NA_real_
        )
      ) %>%
      ungroup()
    
    if (!is.null(time) && isTRUE(input$apply_speed_noise_filter)) {
      df <- df %>%
        group_by(.trip_id) %>%
        mutate(
          .speed_noise_point = (
            (!is.na(.segment_speed_kmh) & .segment_speed_kmh > input$max_realistic_speed_kmh) |
              lead(!is.na(.segment_speed_kmh) & .segment_speed_kmh > input$max_realistic_speed_kmh, default = FALSE)
          )
        ) %>%
        ungroup()
    } else {
      df <- df %>%
        mutate(.speed_noise_point = FALSE)
    }
    
    if (!is.null(altitude) && isTRUE(input$apply_altitude_noise_filter)) {
      df <- df %>%
        mutate(
          .altitude_m_for_noise = suppressWarnings(as.numeric(.data[[altitude]])),
          .altitude_noise_point = !is.na(.altitude_m_for_noise) &
            (
              .altitude_m_for_noise < input$min_realistic_altitude_m |
                .altitude_m_for_noise > input$max_realistic_altitude_m
            )
        )
    } else {
      df <- df %>%
        mutate(
          .altitude_m_for_noise = NA_real_,
          .altitude_noise_point = FALSE
        )
    }
    
    df <- df %>%
      mutate(
        .exclude_from_cuci = .speed_noise_point | .altitude_noise_point
      ) %>%
      group_by(.trip_id) %>%
      mutate(
        .cuci_segment_number = cumsum(lag(.exclude_from_cuci, default = FALSE)) + 1L,
        .cuci_segment_id = paste0(.trip_id, "_S", .cuci_segment_number)
      ) %>%
      ungroup()
    
    speed_noise_point_count <- sum(df$.speed_noise_point, na.rm = TRUE)
    altitude_noise_point_count <- sum(df$.altitude_noise_point, na.rm = TRUE)
    excluded_cuci_point_count <- sum(df$.exclude_from_cuci, na.rm = TRUE)
    
    df <- df %>%
      group_by(.trip_id) %>%
      mutate(.point_index = row_number()) %>%
      ungroup()
    
    pts <- st_as_sf(df, coords = c(lon, lat), crs = 4326, remove = FALSE)
    
    final_trip_ids <- sort(unique(pts$.trip_id))
    final_participant_ids <- sort(unique(pts$participant_id[!is.na(pts$participant_id) & pts$participant_id != ""]))
    
    inactivity_removed_count <- length(unique(inactivity_removed_trips))
    length_removed_count <- length(unique(length_removed_trips))
    
    participant_status <- if (!is.null(participant)) {
      paste0("\nParticipant ID column found: ", participant, ".")
    } else {
      "\nParticipant ID column not found; participant recurrence filter will not be applied."
    }
    
    inactivity_status <- if (!is.null(time) && isTRUE(input$apply_inactivity_filter)) {
      paste0(
        "\nInactivity filter: removed ",
        inactivity_removed_count,
        " trip(s) with gaps > ",
        input$max_time_gap_s,
        " seconds."
      )
    } else if (is.null(time) && isTRUE(input$apply_inactivity_filter)) {
      "\nInactivity filter: not applied because no timestamp column was found."
    } else {
      "\nInactivity filter: disabled."
    }
    
    length_status <- if (isTRUE(input$apply_max_trip_length_filter)) {
      paste0(
        "\nTrip length filter: removed ",
        length_removed_count,
        " trip(s) longer than ",
        input$max_trip_length_km,
        " km."
      )
    } else {
      "\nTrip length filter: disabled."
    }
    
    speed_noise_status <- if (!is.null(time) && isTRUE(input$apply_speed_noise_filter)) {
      paste0(
        "\nSpeed noise filter: flagged ",
        speed_noise_point_count,
        " point(s) as part of unrealistic-speed segment(s) > ",
        input$max_realistic_speed_kmh,
        " km/h. These points are ignored for CUCI detection only."
      )
    } else if (is.null(time) && isTRUE(input$apply_speed_noise_filter)) {
      "\nSpeed noise filter: not applied because no timestamp column was found."
    } else {
      "\nSpeed noise filter: disabled."
    }
    
    altitude_noise_status <- if (!is.null(altitude) && isTRUE(input$apply_altitude_noise_filter)) {
      paste0(
        "\nAltitude noise filter: flagged ",
        altitude_noise_point_count,
        " point(s) outside ",
        input$min_realistic_altitude_m,
        " to ",
        input$max_realistic_altitude_m,
        " m. These points are ignored for CUCI detection only."
      )
    } else if (is.null(altitude) && isTRUE(input$apply_altitude_noise_filter)) {
      "\nAltitude noise filter: not applied because no altitude/elevation column was found."
    } else {
      "\nAltitude noise filter: disabled."
    }
    
    cuci_noise_status <- paste0(
      "\nTotal points ignored for CUCI detection due to point-level noise filters: ",
      excluded_cuci_point_count
    )
    
    status_text(
      paste0(
        "GPS CSV loaded.",
        "\nOriginal points: ", original_point_count,
        "\nPoints after missing-coordinate filtering: ", coordinate_filtered_point_count,
        "\nOriginal trips: ", original_trip_count,
        participant_status,
        inactivity_status,
        length_status,
        speed_noise_status,
        altitude_noise_status,
        cuci_noise_status,
        "\nFinal analysed points: ", nrow(pts),
        "\nFinal analysed trips: ", length(final_trip_ids),
        "\nFinal distinct participants: ", length(final_participant_ids)
      )
    )
    
    pts
  })
  
  points_for_cuci_sf <- reactive({
    pts <- points_sf()
    req(pts)
    
    if (!".exclude_from_cuci" %in% names(pts)) {
      pts$.exclude_from_cuci <- FALSE
    }
    
    if (!".cuci_segment_id" %in% names(pts)) {
      pts$.cuci_segment_id <- pts$.trip_id
    }
    
    pts_clean <- pts %>%
      filter(!.exclude_from_cuci)
    
    if (is.null(pts_clean) || nrow(pts_clean) < 2) return(NULL)
    
    valid_segment_ids <- pts_clean %>%
      st_drop_geometry() %>%
      count(.cuci_segment_id, name = ".segment_point_count") %>%
      filter(.segment_point_count >= 2) %>%
      pull(.cuci_segment_id)
    
    pts_clean <- pts_clean %>%
      filter(.cuci_segment_id %in% valid_segment_ids)
    
    if (nrow(pts_clean) < 2) return(NULL)
    
    pts_clean
  })
  
  osm_network <- eventReactive(input$draw, {
    req(points_sf())
    
    gpkg_path <- NULL
    
    if (!is.null(input$gpkg_file)) {
      gpkg_path_original <- input$gpkg_file$datapath
      gpkg_path <- tempfile(fileext = ".gpkg")
      file.copy(gpkg_path_original, gpkg_path, overwrite = TRUE)
    } else if (file.exists("overijssel.gpkg")) {
      gpkg_path <- "overijssel.gpkg"
    } else {
      msg <- paste(
        "No GeoPackage found. Upload overijssel.gpkg or place it in the same folder as app.R. Current folder:",
        getwd()
      )
      status_text(msg)
      showNotification(msg, type = "error", duration = 15)
      return(NULL)
    }
    
    layers_info <- tryCatch(st_layers(gpkg_path), error = function(e) NULL)
    if (is.null(layers_info)) return(NULL)
    
    layer_names <- layers_info$name
    preferred_layers <- c(
      "gis_osm_roads_free_1",
      "gis_osm_roads_free",
      "roads",
      "road",
      "lines",
      "transport",
      "highways",
      "highway"
    )
    selected_layer <- preferred_layers[preferred_layers %in% layer_names][1]
    
    if (length(selected_layer) == 0 || is.na(selected_layer)) {
      selected_layer <- NULL
      
      for (layer_name in layer_names) {
        test_layer <- tryCatch(
          st_read(gpkg_path, layer = layer_name, quiet = TRUE),
          error = function(e) NULL
        )
        if (is.null(test_layer) || nrow(test_layer) == 0) next
        
        geom_types <- unique(as.character(st_geometry_type(test_layer)))
        
        if (any(geom_types %in% c("LINESTRING", "MULTILINESTRING", "GEOMETRY"))) {
          selected_layer <- layer_name
          break
        }
      }
    }
    
    if (is.null(selected_layer) || is.na(selected_layer)) {
      status_text(
        paste(
          "No line layer found in the GeoPackage. Available layers:",
          paste(layer_names, collapse = ", ")
        )
      )
      return(NULL)
    }
    
    network <- tryCatch(
      st_read(gpkg_path, layer = selected_layer, quiet = TRUE),
      error = function(e) NULL
    )
    if (is.null(network) || nrow(network) == 0) return(NULL)
    
    original_network_feature_count <- nrow(network)
    
    network_lines <- network %>%
      st_make_valid()
    
    geom_types <- unique(as.character(st_geometry_type(network_lines)))
    
    if (any(geom_types %in% c("GEOMETRYCOLLECTION", "MULTILINESTRING"))) {
      network_lines <- network_lines %>%
        st_collection_extract("LINESTRING", warn = FALSE)
    }
    
    network_lines <- network_lines %>%
      st_transform(28992)
    
    pts_m <- points_sf() %>%
      st_transform(28992)
    
    gps_bbox <- st_bbox(pts_m)
    gps_bbox_sfc <- st_as_sfc(gps_bbox)
    
    crop_buffer_m <- input$osm_crop_buffer_m
    if (is.null(crop_buffer_m) || is.na(crop_buffer_m)) {
      crop_buffer_m <- 500
    }
    
    crop_area <- gps_bbox_sfc %>%
      st_buffer(dist = crop_buffer_m)
    
    network_cropped <- tryCatch(
      {
        st_crop(network_lines, st_bbox(crop_area))
      },
      error = function(e) {
        NULL
      }
    )
    
    if (is.null(network_cropped) || nrow(network_cropped) == 0) {
      warning_msg <- paste0(
        "Spatial cropping returned no OSM features. The app will use the full network instead. ",
        "Try increasing the OSM crop buffer."
      )
      
      showNotification(warning_msg, type = "warning", duration = 10)
      
      status_text(
        paste(
          status_text(),
          "\nOSM spatial crop warning:",
          warning_msg
        )
      )
      
      network_cropped <- network_lines
    }
    
    cropped_network_feature_count <- nrow(network_cropped)
    
    status_text(
      paste(
        status_text(),
        "\nUsing GeoPackage layer for map-matching:",
        selected_layer,
        "\nOriginal OSM network features:",
        original_network_feature_count,
        "\nOSM features after spatial crop:",
        cropped_network_feature_count,
        "\nOSM crop buffer:",
        crop_buffer_m,
        "m"
      )
    )
    
    network_cropped
  })
  
  map_matched_points_sf <- eventReactive(input$draw, {
    req(points_sf())
    
    pts <- points_for_cuci_sf()
    network <- osm_network()
    
    if (is.null(pts) || is.null(network) || nrow(network) == 0 || nrow(pts) < 2) return(NULL)
    
    pts_m <- st_transform(pts, 28992)
    nearest_line_index <- st_nearest_feature(pts_m, network)
    
    snapped_points <- lapply(seq_len(nrow(pts_m)), function(i) {
      nearest_line <- network[nearest_line_index[i], ]
      
      nearest_pair <- st_nearest_points(
        st_geometry(pts_m[i, ]),
        st_geometry(nearest_line)
      )
      
      pair_coords <- st_coordinates(nearest_pair)[, c("X", "Y")]
      st_point(pair_coords[2, ])
    })
    
    snapped_sfc <- st_sfc(snapped_points, crs = 28992)
    
    matched_attributes <- pts_m %>%
      st_drop_geometry() %>%
      mutate(matched_osm_feature_index = nearest_line_index)
    
    matched_sf <- st_sf(matched_attributes, geometry = snapped_sfc, crs = 28992)
    
    st_transform(matched_sf, 4326)
  })
  
  raw_display_trace_sf <- reactive({
    req(points_sf())
    make_lines_from_points(points_sf(), id_col = ".trip_id", crs_out = 4326)
  })
  
  raw_trace_sf <- reactive({
    pts_clean <- points_for_cuci_sf()
    if (is.null(pts_clean)) return(NULL)
    make_lines_from_points(pts_clean, id_col = ".cuci_segment_id", crs_out = 4326)
  })
  
  map_matched_trace_sf <- reactive({
    matched_pts <- map_matched_points_sf()
    if (is.null(matched_pts)) return(NULL)
    make_lines_from_points(matched_pts, id_col = ".cuci_segment_id", crs_out = 4326)
  })
  
  all_gap_segments_sf <- reactive({
    raw_lines <- raw_trace_sf()
    matched_lines <- map_matched_trace_sf()
    
    if (is.null(raw_lines) || is.null(matched_lines)) return(NULL)
    
    detection_ids <- intersect(raw_lines$detection_id, matched_lines$detection_id)
    if (length(detection_ids) == 0) return(NULL)
    
    gap_list <- lapply(detection_ids, function(id) {
      raw_line <- raw_lines %>% filter(detection_id == id)
      matched_line <- matched_lines %>% filter(detection_id == id)
      
      raw_interp <- interpolate_line_points(raw_line, spacing_m = input$dtw_spacing)
      matched_interp <- interpolate_line_points(matched_line, spacing_m = input$dtw_spacing)
      
      if (is.null(raw_interp) || is.null(matched_interp)) return(NULL)
      
      raw_coords <- st_coordinates(raw_interp)[, c("X", "Y")]
      matched_coords <- st_coordinates(matched_interp)[, c("X", "Y")]
      
      path <- dtw_path_basic(
        raw_coords = raw_coords,
        matched_coords = matched_coords,
        window_size = input$dtw_window_size
      )
      
      extract_gap_lines_from_dtw(
        raw_interp = raw_interp,
        dtw_path = path,
        distance_threshold = input$distance_threshold,
        consecutive_threshold = input$consecutive_threshold,
        min_gap_length = input$min_gap_length
      )
    })
    
    gap_list <- gap_list[!sapply(gap_list, is.null)]
    if (length(gap_list) == 0) return(NULL)
    
    do.call(rbind, gap_list) %>%
      st_transform(4326)
  })
  
  gap_segments_sf <- reactive({
    gaps <- all_gap_segments_sf()
    if (is.null(gaps) || nrow(gaps) == 0) return(NULL)
    
    gaps_with_repetition <- add_repetition_counts_to_gaps(
      gaps_sf = gaps,
      repeat_distance_m = input$repeat_distance_m
    )
    
    gaps_filtered <- gaps_with_repetition %>%
      filter(repeated_trip_count >= input$min_repeated_trips)
    
    participant_filter_available <- "repeated_participant_count" %in% names(gaps_filtered) &&
      any(!is.na(gaps_filtered$repeated_participant_count))
    
    if (participant_filter_available) {
      gaps_filtered <- gaps_filtered %>%
        filter(repeated_participant_count >= input$min_repeated_participants)
    }
    
    gaps_filtered
  })
  
  detailed_gap_table <- reactive({
    gaps <- gap_segments_sf()
    
    if (is.null(gaps) || nrow(gaps) == 0) {
      return(data.frame(
        trip_id = character(),
        participant_id = character(),
        gap_id = integer(),
        travel_direction = character(),
        bearing_deg = numeric(),
        gap_length_m = numeric(),
        max_distance_m = numeric(),
        mean_distance_m = numeric(),
        consecutive_dtw_points = integer(),
        repeat_cluster_id = integer(),
        repeated_trip_count = integer(),
        repeated_participant_count = integer(),
        repeated_gap_count = integer(),
        repeated_trip_ids = character(),
        repeated_participant_ids = character()
      ))
    }
    
    gaps %>%
      st_drop_geometry() %>%
      select(any_of(c(
        "trip_id",
        "participant_id",
        "gap_id",
        "travel_direction",
        "bearing_deg",
        "gap_length_m",
        "max_distance_m",
        "mean_distance_m",
        "consecutive_dtw_points",
        "repeat_cluster_id",
        "repeated_trip_count",
        "repeated_participant_count",
        "repeated_gap_count",
        "repeated_trip_ids",
        "repeated_participant_ids",
        "first_raw_index",
        "last_raw_index"
      )))
  })
  
  candidate_summary_table <- reactive({
    gaps <- gap_segments_sf()
    
    if (is.null(gaps) || nrow(gaps) == 0) {
      return(data.frame(
        repeat_cluster_id = integer(),
        candidate_id = character(),
        recurrence_level = character(),
        suggested_action = character(),
        dominant_direction = character(),
        candidate_length_m = numeric(),
        trips_supporting = integer(),
        participants_supporting = integer(),
        strongest_mismatch_m = numeric(),
        average_mismatch_m = numeric(),
        source_trips = character()
      ))
    }
    
    gaps %>%
      st_drop_geometry() %>%
      group_by(repeat_cluster_id) %>%
      summarise(
        candidate_id = paste0("CUCI-", sprintf("%02d", repeat_cluster_id[1])),
        candidate_length_m = round(max(gap_length_m, na.rm = TRUE), 1),
        trips_supporting = max(repeated_trip_count, na.rm = TRUE),
        participants_supporting = ifelse(
          all(is.na(repeated_participant_count)),
          NA_integer_,
          max(repeated_participant_count, na.rm = TRUE)
        ),
        strongest_mismatch_m = round(max(max_distance_m, na.rm = TRUE), 1),
        average_mismatch_m = round(mean(mean_distance_m, na.rm = TRUE), 1),
        dominant_direction = {
          dirs <- travel_direction[!is.na(travel_direction) & travel_direction != ""]
          if (length(dirs) == 0) NA_character_ else names(sort(table(dirs), decreasing = TRUE))[1]
        },
        source_trips = paste(sort(unique(trip_id)), collapse = ", "),
        .groups = "drop"
      ) %>%
      mutate(
        recurrence_level = case_when(
          trips_supporting >= 3 & !is.na(participants_supporting) & participants_supporting >= 2 ~ "High",
          trips_supporting >= 2 ~ "Medium",
          TRUE ~ "Low"
        ),
        suggested_action = case_when(
          recurrence_level == "High" ~ "Prioritise for field or satellite validation",
          recurrence_level == "Medium" ~ "Review visually before validation",
          TRUE ~ "Treat cautiously; possible GPS noise"
        )
      ) %>%
      select(
        repeat_cluster_id,
        candidate_id,
        recurrence_level,
        suggested_action,
        dominant_direction,
        candidate_length_m,
        trips_supporting,
        participants_supporting,
        strongest_mismatch_m,
        average_mismatch_m,
        source_trips
      )
  })
  
  output$status <- renderText({
    status_text()
  })
  
  output$map <- renderLeaflet({
    leaflet(options = leafletOptions(preferCanvas = TRUE)) %>%
      addTiles(group = "OpenStreetMap") %>%
      addProviderTiles(
        providers$Esri.WorldImagery,
        group = "Satellite imagery"
      ) %>%
      addProviderTiles(
        providers$CartoDB.Positron,
        group = "Light map"
      ) %>%
      setView(lng = 6.8937, lat = 52.2215, zoom = 13) %>%
      addLayersControl(
        baseGroups = c("OpenStreetMap", "Satellite imagery", "Light map"),
        overlayGroups = c(
          "OSM network",
          "GPS points",
          "Raw GPS traces",
          "Map-matched trace",
          "Detected CUCI candidates",
          "Selected candidate"
        ),
        options = layersControlOptions(collapsed = FALSE)
      )
  })
  
  observeEvent(input$draw, {
    
    withProgress(
      message = "Running CUCI detection",
      detail = "Loading GPS CSV and applying selected preprocessing filters...",
      value = 0,
      {
        setProgress(value = 0.10, detail = "Loading GPS CSV and applying selected preprocessing filters...")
        pts <- points_sf()
        req(pts)
        
        setProgress(value = 0.25, detail = "Loading OSM GeoPackage and preparing the network...")
        network <- osm_network()
        
        setProgress(value = 0.40, detail = "Creating raw GPS polylines...")
        raw_trace <- raw_display_trace_sf()
        
        setProgress(value = 0.55, detail = "Creating map-matched traces...")
        matched_trace <- map_matched_trace_sf()
        
        setProgress(value = 0.70, detail = "Comparing traces with DTW and detecting candidate gaps...")
        gaps <- gap_segments_sf()
        
        setProgress(value = 0.85, detail = "Updating map layers and candidate table...")
        
        m <- leafletProxy("map")
        
        m %>%
          clearGroup("OSM network") %>%
          clearGroup("GPS points") %>%
          clearGroup("Raw GPS traces") %>%
          clearGroup("Map-matched trace") %>%
          clearGroup("Detected CUCI candidates") %>%
          clearGroup("Selected candidate") %>%
          removeControl("trace_legend")
        
        if (!is.null(network)) {
          network_display <- st_transform(network, 4326)
          
          m %>%
            addPolylines(
              data = network_display,
              color = "#555555",
              weight = 2,
              opacity = 0.65,
              group = "OSM network"
            )
          
          if (!isTRUE(input$show_osm_network)) {
            m %>% hideGroup("OSM network")
          }
        }
        
        if (isTRUE(input$show_points)) {
          m %>%
            addCircleMarkers(
              data = pts,
              radius = 4,
              stroke = FALSE,
              fillOpacity = 0.8,
              group = "GPS points",
              popup = ~paste(
                "Trip:", .trip_id,
                "<br>Participant:", participant_id,
                "<br>Point:", .point_index,
                "<br>Trip length:", round(.trip_length_m, 1), "m"
              )
            )
        }
        
        if (isTRUE(input$show_raw_trace) && !is.null(raw_trace)) {
          m %>%
            addPolylines(
              data = raw_trace,
              color = "blue",
              weight = 4,
              opacity = 0.8,
              group = "Raw GPS traces",
              popup = ~paste(
                "Raw GPS trace",
                "<br>Trip:", trip_id,
                "<br>Participant:", participant_id
              )
            )
        }
        
        if (isTRUE(input$show_map_matched) && !is.null(matched_trace)) {
          m %>%
            addPolylines(
              data = matched_trace,
              color = "red",
              weight = 4,
              opacity = 0.8,
              group = "Map-matched trace",
              popup = ~paste(
                "Map-matched trace",
                "<br>Trip:", trip_id,
                "<br>Participant:", participant_id
              )
            )
        }
        
        if (isTRUE(input$show_gap_segments) && !is.null(gaps) && nrow(gaps) > 0) {
          m %>%
            addPolylines(
              data = gaps,
              color = "orange",
              weight = 7,
              opacity = 0.95,
              group = "Detected CUCI candidates",
              popup = ~paste(
                "Detected CUCI candidate gap",
                "<br>Trip:", trip_id,
                "<br>Participant:", participant_id,
                "<br>Candidate group:", paste0("CUCI-", sprintf("%02d", repeat_cluster_id)),
                "<br>Travel direction:", travel_direction, "(", bearing_deg, "°)",
                "<br>Gap length:", gap_length_m, "m",
                "<br>Largest GPS-map mismatch:", max_distance_m, "m",
                "<br>Average GPS-map mismatch:", mean_distance_m, "m",
                "<br>Supporting trips:", repeated_trip_count,
                "<br>Supporting participants:", repeated_participant_count,
                "<br>Source trips:", repeated_trip_ids
              )
            )
        }
        
        m %>%
          addLegend(
            layerId = "trace_legend",
            position = "bottomright",
            colors = c("blue", "red", "orange"),
            labels = c(
              "GPS polylines",
              "Map-matched trace",
              "Detected CUCI candidate"
            ),
            title = "Trace interpretation",
            opacity = 1
          )
        
        bbox <- st_bbox(pts)
        
        m %>%
          fitBounds(
            lng1 = bbox["xmin"],
            lat1 = bbox["ymin"],
            lng2 = bbox["xmax"],
            lat2 = bbox["ymax"]
          )
        
        candidate_summary <- candidate_summary_table()
        detailed_summary <- detailed_gap_table()
        
        status_text(
          paste(
            status_text(),
            "\nDetection complete.",
            "\nCandidate locations after repetition filters:",
            nrow(candidate_summary),
            "\nDetailed gap segments within these candidate locations:",
            nrow(detailed_summary)
          )
        )
        
        setProgress(value = 1, detail = "Done.")
      }
    )
  })
  
  output$gap_table <- renderDT({
    table_to_show <- candidate_summary_table() %>%
      rename(
        "Cluster ID" = repeat_cluster_id,
        "Candidate ID" = candidate_id,
        "Recurrence level" = recurrence_level,
        "Suggested action" = suggested_action,
        "Dominant travel direction" = dominant_direction,
        "Length of candidate (m)" = candidate_length_m,
        "Supporting trips" = trips_supporting,
        "Supporting participants" = participants_supporting,
        "Largest GPS-map mismatch (m)" = strongest_mismatch_m,
        "Average GPS-map mismatch (m)" = average_mismatch_m,
        "Source trips" = source_trips
      )
    
    datatable(
      table_to_show,
      selection = "single",
      rownames = FALSE,
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        columnDefs = list(
          list(visible = FALSE, targets = 0)
        )
      )
    )
  })
  
  observeEvent(input$gap_table_rows_selected, {
    selected_row <- input$gap_table_rows_selected
    
    req(length(selected_row) == 1)
    
    summary_table <- candidate_summary_table()
    req(nrow(summary_table) >= selected_row)
    
    selected_cluster_id <- summary_table$repeat_cluster_id[selected_row]
    
    gaps <- gap_segments_sf()
    req(!is.null(gaps), nrow(gaps) > 0)
    
    selected_candidate <- gaps %>%
      filter(repeat_cluster_id == selected_cluster_id)
    
    req(nrow(selected_candidate) > 0)
    
    selected_candidate_4326 <- st_transform(selected_candidate, 4326)
    
    selected_bbox_geom <- selected_candidate_4326 %>%
      st_transform(28992) %>%
      st_union() %>%
      st_buffer(75) %>%
      st_transform(4326)
    
    bbox <- st_bbox(selected_bbox_geom)
    
    leafletProxy("map") %>%
      clearGroup("Selected candidate") %>%
      addPolylines(
        data = selected_candidate_4326,
        color = "black",
        weight = 10,
        opacity = 0.85,
        group = "Selected candidate",
        popup = ~paste(
          "Selected CUCI candidate",
          "<br>Candidate group:", paste0("CUCI-", sprintf("%02d", repeat_cluster_id)),
          "<br>Trip:", trip_id,
          "<br>Travel direction:", travel_direction, "(", bearing_deg, "°)",
          "<br>Gap length:", gap_length_m, "m",
          "<br>Supporting trips:", repeated_trip_count,
          "<br>Supporting participants:", repeated_participant_count
        )
      ) %>%
      fitBounds(
        lng1 = bbox["xmin"],
        lat1 = bbox["ymin"],
        lng2 = bbox["xmax"],
        lat2 = bbox["ymax"]
      )
  })
  
  output$download_gap_csv <- downloadHandler(
    filename = function() {
      paste0("cuci_candidate_summary_", Sys.Date(), ".csv")
    },
    content = function(file) {
      summary <- candidate_summary_table()
      readr::write_csv(summary, file)
    }
  )
  
  output$download_gap_geojson <- downloadHandler(
    filename = function() {
      paste0("cuci_candidate_segments_", Sys.Date(), ".geojson")
    },
    content = function(file) {
      gaps <- gap_segments_sf()
      if (is.null(gaps) || nrow(gaps) == 0) {
        empty_sf <- st_sf(
          trip_id = character(),
          participant_id = character(),
          gap_id = integer(),
          travel_direction = character(),
          bearing_deg = numeric(),
          gap_length_m = numeric(),
          max_distance_m = numeric(),
          mean_distance_m = numeric(),
          consecutive_dtw_points = integer(),
          repeat_cluster_id = integer(),
          repeated_trip_count = integer(),
          repeated_participant_count = integer(),
          repeated_gap_count = integer(),
          repeated_trip_ids = character(),
          repeated_participant_ids = character(),
          geometry = st_sfc(crs = 4326)
        )
        st_write(empty_sf, file, driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
      } else {
        st_write(gaps, file, driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
      }
    }
  )
  
  output$download_cleaned_gps_csv <- downloadHandler(
    filename = function() {
      paste0("analysed_gps_points_", Sys.Date(), ".csv")
    },
    content = function(file) {
      pts <- points_sf()
      if (is.null(pts) || nrow(pts) == 0) {
        readr::write_csv(data.frame(), file)
      } else {
        pts_out <- pts %>%
          st_drop_geometry() %>%
          select(any_of(c(
            ".trip_id", "participant_id", ".point_index", ".timestamp_order", ".row_order",
            ".time_gap_s", ".has_long_inactivity", ".trip_length_m",
            ".segment_distance_m", ".segment_time_s", ".segment_speed_kmh",
            ".speed_noise_point", ".altitude_m_for_noise", ".altitude_noise_point",
            ".exclude_from_cuci", ".cuci_segment_id",
            "trip_id", "unique_trip_id", "source_file",
            "timestamp", "time", "latitude", "longitude", "lat", "lon", "elevation_m"
          )))
        readr::write_csv(pts_out, file)
      }
    }
  )
}

shinyApp(ui, server)
