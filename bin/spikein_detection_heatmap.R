# Load necessary libraries
library(reshape2)
library(scales)
library(argparse)
library(tidyverse)

parser <- ArgumentParser()
parser$add_argument("--input", type = "character", nargs="+", help = "Input Spikein Count CSV file of all samples. Columns are SampleID, SpikeinID and Count.")
parser$add_argument("--expected", type = "character", help = "CSV containing location information for a SampleID. Each row should contain a SampleID, Plate and Well. Columns are SampleID, Plate and Well.")
parser$add_argument("--spikein-info", type = "character", help = "CSV containing what spike-in is expected at which location. Columns are SpikeinID and Well.")
parser$add_argument("--output", type = "character", help = "Output PNG file")

# Parsing arguments
args <- parser$parse_args()

# TEST DATA
# WHEN TESTING, DO NOT FORGET TO MELT THE DATA.
# counts_data <- tibble(
#   SampleID = c("SampleID_1", "SampleID_1", "SampleID_2", "SampleID_2", "SampleID_3", "SampleID_3", "SampleID_3", "SampleID_4", "SampleID_4"),
#   SpikeinID = c("SpikeinID_1", "SpikeinID_2", "SpikeinID_1", "SpikeinID_2", "SpikeinID_3", "SpikeinID_2", "SpikeinID_19", "SpikeinID_4", "SpikeinID_4"),
#   Count = c(100, 200, 300, 400, 500, 600, 1100, 700, 800)
# )
#
# expected_data <- tibble(
#   SampleID = c("SampleID_1", "SampleID_2", "SampleID_3", "SampleID_4"),
#   Plate = c("Plate_1"),
#   Well = c("A01", "B02", "C03", "D01")
# )
#
# spikein_info <- tibble(
#   SpikeinID = sprintf("SpikeinID_%d", 1:96),
#   Well = expand.grid(LETTERS[1:8], str_pad(1:12, 2, pad = "0")) |>
#     apply(1, paste, collapse = "") |>
#     as.character()
# )
# WHEN TESTING, DO NOT FORGET TO MELT THE DATA.


validate_data <- function(counts_data, expected_data, spikein_info) {
  if (nrow(counts_data) == 0) {
    stop("No spikein count data was detected.")
  }

  # Check that the files have the correct columns
  if (!all(c("SampleID", "SpikeinID", "Count") %in% colnames(counts_data))) {
    stop("The spikein count data does not contain the correct columns.")
  }

  if (!all(c("SampleID", "Plate", "Well") %in% colnames(expected_data))) {
    stop("The expected data does not contain the correct columns.")
  }

  if (!all(c("SpikeinID", "Well") %in% colnames(spikein_info))) {
    stop("The spikein info does not contain the correct columns.")
  }
}

melt_data <- function(counts_data) {

  # Transform the data to logaritmic (NOTE: use raw read counts for now)
  # counts_data$Count <- log10(counts_data$Count + 1)

  # Transform the data from long to wide format
  wide_data <- dcast(counts_data, SampleID ~ SpikeinID, value.var = "Count", fun.aggregate = sum)

  # Melt the data for ggplot2
  melted_data <- melt(wide_data, id.vars = "SampleID")

  # Return melted data
  return (melted_data)
}

plot_spikein_detection_heatmap_by_sampleid <- function(melted_data, expected_data, spikein_info) {

  # Create a 96x96 grid based on Well positions
  sample_ids <- spikein_info$Well
  spikein_ids <- spikein_info$Well

  # Create a complete grid of all 96x96 combinations of Well positions
  full_grid <- expand.grid(ExpectedSpikeinID = spikein_ids, variable = sample_ids) %>%
    left_join(spikein_info, by = c("ExpectedSpikeinID" = "Well")) %>%
    left_join(spikein_info, by = c("variable" = "Well"), suffix = c(".expected", ".variable"))

  # Join the spikein_info to associate each Well with its corresponding SpikeinID
  joined_data <- expected_data %>%
    left_join(spikein_info, by = c("Well" = "Well")) %>%
    dplyr::rename(ExpectedSpikeinID = SpikeinID)

  # Transform the melted data by joining with joined_data
  transformed_data <- melted_data %>%
    left_join(joined_data, by = c("SampleID")) %>%
    group_by(ExpectedSpikeinID) %>%
    mutate(total = sum(value),
           value = if_else(total == 0, 0, (value / total * 100))) %>%
    select(-total)  # Remove the temporary column if not needed

  # Merge transformed data with full grid to ensure complete 96x96 layout
  transformed_data <- full_grid %>%
    left_join(transformed_data, by = c(
      "SpikeinID.expected" = "ExpectedSpikeinID",
      "SpikeinID.variable" = "variable")
    )

  # Use 'Well' positions as labels for both axes
  transformed_data$variable <- factor(transformed_data$variable, levels = rev(sample_ids))
  transformed_data$ExpectedSpikeinID <- factor(transformed_data$ExpectedSpikeinID, levels = spikein_ids)

  # Generate the heatmap with swapped x and y axes
  g <- ggplot(transformed_data, aes(x = ExpectedSpikeinID, y = variable, fill = value)) +
    geom_tile() +
    scale_fill_gradient(
      name = "SDSI (%)",
      low = "white",
      high = "#0072B2",
      limits = c(0, 100),
      na.value = "white"
    ) +
    labs(
      x = "SDSI 1 → 96\nExpected Synthetic DNA spike-in",
      y = "Synthetic DNA spike-in\nSDSI 96 → 1"
    ) +
    theme_minimal(base_size = 8) +
    theme(
      # Customize axis and legend to match the reference style
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.length = unit(0.2, "cm"),
      axis.ticks = element_line(linewidth = 0.2),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),  # Black border around the plot
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.title.x = element_text(vjust = -1),
      axis.title.y = element_text(angle = 90, vjust = 1),
      legend.key.size = unit(0.5, "cm")
    ) +
    guides(
      fill = guide_colorbar(
        title.position = "bottom",
        title.hjust = 0.5,
        barwidth = unit(0.5, "cm"),
        barheight = unit(15, "cm")
      )
    ) +
    coord_fixed(ratio = 1)  # Ensure a square aspect ratio for the 96x96 layout

  # Add a black border around wells with detected spikein (value > 0 & value < 2)
  # IMPORTANT: This should be tun-able.
  g <- g + geom_tile(
    data = transformed_data %>% filter(value > 0 & value < 2 & SpikeinID.expected != SpikeinID.variable),
    aes(x = ExpectedSpikeinID, y = variable),
    fill = NA, color = "black", size = 0.25
  )

  return(g)
}


# This is a helper function for the functions below. It will create
# a matrix that represents a plate, and plot the metric value for each well
# given the metric calculation function.
create_plate_matrix_with_metric_value <- function(melted_data, expected_data, spikein_info, plate_id) {

  # Create the matrix
  plate_matrix <- matrix(NA, nrow = 8, ncol = 12)
  rownames(plate_matrix) <- LETTERS[1:8]
  colnames(plate_matrix) <- str_pad(1:12, 2, pad = "0")

  # Track whether there was expected or unexpected spikein at this well.
  unexpected_count <- list()

  # Fill the matrix with the calculated metric for each sample,
  # based on the well that the sample was in on the plate. The melted
  # data should be used to determine how many spikeins were detected and
  # at what quantity. Spikein_info will tell us what spikein is expected
  # at the well, and expected_data will tell us what sample was in the well
  # for the plate.
  expected_data_by_plate <- expected_data[expected_data$Plate == plate_id,]
  for (i in 1:nrow(expected_data_by_plate)) {
    sample_id <- expected_data_by_plate$SampleID[i]
    well <- expected_data_by_plate$Well[i]
    spikein_id <- spikein_info$SpikeinID[spikein_info$Well == well]
    row_id <- match(substr(well, 1, 1), rownames(plate_matrix))
    col_id <- match(substr(well, 2, 3), colnames(plate_matrix))
    # Get the expected count. This should be the number of spikeins we want for a SampleID.
    expected_count <- melted_data$value[melted_data$SampleID == sample_id & melted_data$variable == spikein_id]

    # If the expected count is empty (length = 0), then we have no spike-ins for this well.
    if (length(expected_count) == 0) {
      expected_count <- 0
    }

    # The the unexpected count. This is the number of Spikeins detected that do not match the expected spikein.
    unexpected_count[[well]] <- sum(melted_data$value[melted_data$SampleID == sample_id & melted_data$variable != spikein_id])

    # IMPORTANT: Calculate the metric!
    # 100% means we _only_ had "unexpected" spike-in
    # 0% means we _only_ had "expected" spike-in
    metric <- if (unexpected_count[[well]] == 0 && expected_count == 0) { 0 } else {
      unexpected_count[[well]] / (unexpected_count[[well]] + expected_count) * 100.
    }

    # Apply the value to the well in the plate matrix
    plate_matrix[row_id, col_id] <- metric
  }

  return (plate_matrix)
}

# This is a helper function for the functions below. It will create
# a matrix that represents a plate, and plot the metric value for each well
# given the metric calculation function.
create_plate_matrix_with_contaminant_tracing_for_sampleID <- function(melted_data, expected_data, spikein_info, sample_id) {

  # Create the matrix
  plate_matrix <- matrix(NA, nrow = 8, ncol = 12)
  rownames(plate_matrix) <- LETTERS[1:8]
  colnames(plate_matrix) <- str_pad(1:12, 2, pad = "0")

  # Track whether there was expected or unexpected spikein at this well.
  unexpected_count <- list()

  # Get the plate and well for the sample
  plate_id <- expected_data$Plate[expected_data$SampleID == sample_id]
  sample_id_well <- expected_data$Well[expected_data$SampleID == sample_id]
  sample_specific_melted_data <- melted_data[melted_data$SampleID == sample_id,]

  # Get the unexpected spike-in IDs that were found in the sample and
  # match the ID against the spikein_info to get the well it is supposed
  # to be in.
  expected_spikein_id <- spikein_info$SpikeinID[spikein_info$Well == sample_id_well]
  unexpected_spikein_ids <- spikein_info$SpikeinID[spikein_info$SpikeinID != expected_spikein_id]
  unexpected_spikein_ids_in_sample_id_well <- unexpected_spikein_ids[
    unexpected_spikein_ids %in% sample_specific_melted_data$variable]

  # Get the wells for the unexpected spike-ins
  unexpected_wells <- spikein_info$Well[spikein_info$SpikeinID %in% unexpected_spikein_ids_in_sample_id_well]

  # Fill the matrix with the proportion of spike-in reads
  # that were found in the `sample_id` well, that should
  # have been found in the `unexpected_wells` well. The metric
  # should summarize the proportion of spike-in from the wrong
  # well that made up the amount of spike-in in the `sample_id` well.
  #
  # The metric should be the proportion of the spike-in reads that
  # made up the total unexpected spike-in reads in the `sample_id` well.
  sample_specific_melted_data_unexpected <- sample_specific_melted_data %>%
    filter(variable %in% unexpected_spikein_ids_in_sample_id_well)
  total_spikein_in_sample_id_well <- sum(sample_specific_melted_data_unexpected$value)

  # Include a table of the unexpected spike-ins in the sample
  unexpected_spikein_table <- NULL

  for (well in unexpected_wells) {
    spikein_id <- spikein_info$SpikeinID[spikein_info$Well == well]
    count_of_unexpected_spikein_in_sample_id_well <- sample_specific_melted_data$value[
      sample_specific_melted_data$variable == spikein_id]

    row_id <- match(substr(well, 1, 1), rownames(plate_matrix))
    col_id <- match(substr(well, 2, 3), colnames(plate_matrix))

    # IMPORTANT: Calculate the metric!
    # 100% means we _only_ had "unexpected" spike-in
    # 0% means we _only_ had "expected" spike-in
    metric <- (count_of_unexpected_spikein_in_sample_id_well / total_spikein_in_sample_id_well) * 100.

    # Apply the value to the well in the plate matrix
    plate_matrix[row_id, col_id] <- metric

    # Spikein tibble
    unexpected_spikein_table <- bind_rows(
      unexpected_spikein_table,
      tibble(Well = well, Count = count_of_unexpected_spikein_in_sample_id_well, Metric = metric)
    )
  }

  # Add a -1 where the sample id well is located
  row_id <- match(substr(sample_id_well, 1, 1), rownames(plate_matrix))
  col_id <- match(substr(sample_id_well, 2, 3), colnames(plate_matrix))
  plate_matrix[row_id, col_id] <- -Inf

  result_list <- list(plate_matrix = plate_matrix, unexpected_spikein_table = unexpected_spikein_table)
  return (result_list)
}

# Plot of the plate where presence / absence of contamination outlined by green / red around a well, respectively.
# The metric used will be: unexpected spike-in / total spike-in (unexpected + expected)
# The border will appear red after surpassing a percent amount of unexpected spikein defined by `border_threshold` (default: `1%`).
plot_spikein_detection_plate_heatmap <- function(melted_data, expected_data, spikein_info, plate_id, border_threshold = 1.) {

  # Create a plot of a 96-well plate that shows the presence / absence of contamination
  # using a green or red border around the well, respectively.
  # The metric used will be: unexpected spike-in / total spike-in (unexpected + expected).
  # The plate_id is the ID of the plate to plot.
  plate_matrix <- create_plate_matrix_with_metric_value(melted_data, expected_data, spikein_info, plate_id)

  # Transform the plate_matrix into a dataframe to plot
  plate_df <- as_tibble(plate_matrix)
  plate_df$row <- rownames(plate_matrix)

  # Arrange the data.frame so that we can plot the metrics appropriately.
  # This will need to incorporate the unexpected_counts as the border
  # is going to be determined by the presence of unexpected spikeins.
  # Use pivot_longer to turn into long form
  plate_df_long <- plate_df |>
    pivot_longer(cols = -row, names_to = "col", values_to = "value") %>%
    mutate(row = factor(row, levels = rev(LETTERS[1:8])),
           col = factor(col, levels = str_pad(1:12, 2, pad = "0")))


  # Create the heatmap of the 96-well plate.
  # Use a green border for contamination and a red border for no contamination.
  # Color each well based on the metric calculated above with white close to 0 and orange closer to 1.
  # Column names (x-axis) should display on top of the plot.
  # Row names (y-axis) should display on the left side going from A-H, top to bottom.
  # The plot should have a white background.
  # The plot should fill by `value` and range from 0 - 1, white to orange.
  g <- ggplot(plate_df_long, aes(x = col, y = row)) +
    geom_tile() +
    geom_tile(aes(fill = value), colour = "black") +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 100)) +
    theme_minimal(base_size = 10) +
    scale_x_discrete(position = "top") +
    theme(panel.background = element_rect(fill = "white", colour = "white"),
          plot.background = element_rect(fill = "white", colour = "white")) +
    labs(x = element_blank(), y = element_blank(), fill = "SDSI (%)", title = plate_id)

  # Make x and y axis labels bold and larger
  g <- g + theme(axis.text.x = element_text(face = "bold", size = 12),
                 axis.text.y = element_text(face = "bold", size = 12))

  # add a border around the wells, green if it is clean, red if it is contaminated.
  # make sure that NA values are not plotted.
  # The red border will only appear after surpassing `threshold`
  sample_found <- !(plate_df_long$value %>% is.na())
  g <- g + geom_rect(data = plate_df_long[sample_found & plate_df_long$value <= border_threshold,], aes(xmin = as.numeric(col) - 0.5, xmax = as.numeric(col) + 0.5, ymin = as.numeric(row) - 0.5, ymax = as.numeric(row) + 0.5), fill = NA, colour = "green", size = 0.75) +
    geom_rect(data = plate_df_long[sample_found & plate_df_long$value > border_threshold,], aes(xmin = as.numeric(col) - 0.5, xmax = as.numeric(col) + 0.5, ymin = as.numeric(row) - 0.5, ymax = as.numeric(row) + 0.5), fill = NA, colour = "red", size = 0.75)

  g <- g + guides(
    fill = guide_colorbar(
      title.position = "bottom",
      title.hjust = 0.5,
      barwidth = unit(0.5, "cm"),
      barheight = unit(14, "cm")
    )
  )

  return (g)
}

# This function will plot a heatmap of the plate where the contamination originated for a given sample.
plot_contamination_origin_by_sample <- function(melted_data, expected_data, spikein_info, sample_id) {

  # Create the plate matrix for the sample
  contaminant_tracing_results <- create_plate_matrix_with_contaminant_tracing_for_sampleID(melted_data, expected_data, spikein_info, sample_id)
  plate_matrix <- contaminant_tracing_results$plate_matrix

  # Transform the plate_matrix into a dataframe to plot
  plate_df <- as_tibble(plate_matrix)
  plate_df$row <- rownames(plate_matrix)

  # Arrange the data.frame so that we can plot the metrics appropriately.
  # This will need to incorporate the unexpected_counts as the border
  # is going to be determined by the presence of unexpected spikeins.
  # Use pivot_longer to turn into long form
  plate_df_long <- plate_df |>
    pivot_longer(cols = -row, names_to = "col", values_to = "value") %>%
    mutate(row = factor(row, levels = rev(LETTERS[1:8])),
           col = factor(col, levels = str_pad(1:12, 2, pad = "0")))

  # Create the heatmap of the 96-well plate.
  # If there is a -Inf, that is the well of the sample id. This should
  # be filled in with black.
  g <- ggplot(plate_df_long, aes(x = col, y = row)) +
    geom_tile(aes(fill = value), colour = "black") +
    scale_fill_gradientn(
      colors = c("black", "white", "red"),
      values = scales::rescale(c(-Inf, 0, 100)),
      limits = c(0, 100)
    ) +
    theme_minimal(base_size = 10) +
    scale_x_discrete(position = "top") +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.background = element_rect(fill = "white", colour = "white"),
      axis.text.x = element_text(face = "bold", size = 12),
      axis.text.y = element_text(face = "bold", size = 12)
    ) +
    labs(x = element_blank(), y = element_blank(), fill = "SDSI (%)", title = sample_id, subtitle = plate_id)

  # Update the legend with a larger color bar
  g <- g + guides(
    fill = guide_colorbar(
      title.position = "bottom",
      title.hjust = 0.5,
      barwidth = unit(0.5, "cm"),
      barheight = unit(14, "cm")
    )
  )

  # Build return list
  unexpected_spikein_table <- contaminant_tracing_results$unexpected_spikein_table
  result_list <- list(plate_graphic = g, unexpected_spikein_table = unexpected_spikein_table)

  return (result_list)
}

args <- list()

args$input <- "~/Documents/perfect_case.csv"
args$input <- "~/Documents/contamination_case.csv"
args$input <- "~/Documents/low_contamination_case.csv"
args$expected <- "~/Documents/expected_data.csv"
args$spikein_info <- "~/Documents/spikein_info.csv"

# Load the data
counts_data <- args$input |> map_dfr(read_csv)
expected_data <- read_csv(args$expected)
spikein_info <- read_csv(args$spikein_info)

# Write concatenated file
# write.csv(counts_data, file = "spikein_counts_data.csv", quote = FALSE, row.names = FALSE)
validate_data(counts_data, expected_data, spikein_info)
melted_data <- melt_data(counts_data)

# Save the heatmap with a white background
g <- plot_spikein_detection_heatmap_by_sampleid(melted_data, expected_data, spikein_info)
ggsave(args$output, plot = g, bg = "white", width = 10, height = 8, dpi = 300)

g <- plot_spikein_detection_plate_heatmap(melted_data, expected_data, spikein_info, "Plate_1")
ggsave(args$output, plot = g, bg = "white", width = 10, height = 8, dpi = 300)

result_list <- plot_contamination_origin_by_sample(melted_data, expected_data, spikein_info, "SampleID_2")
ggsave(args$output, plot = result_list$plate_graphic, bg = "white", width = 10, height = 8, dpi = 300)
write_csv(result_list$unexpected_spikein_table, file.path("~/Documents", "unexpected_spikein_table.csv"))

# Save the heatmap with a white background
plot_contamination_origin_by_sample(melted_data, expected_data, spikein_info, "SampleID_3")
plot_spikein_detection_plate_heatmap(melted_data, expected_data, spikein_info, "Plate_1")

# DEBUGGING
melted_data %>% filter(SampleID == "SampleID_3")
spikein_info %>% filter(Well == "C03")
