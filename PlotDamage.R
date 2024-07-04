#!/usr/bin/env Rscript

# Plot damage for a given taxa
# Part of BamDam, written by Bianca De Sanctis, July 3 2024

library(optparse)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(grid)
pdf(NULL)

# Define command-line arguments
option_list <- list(
  make_option(c("-f", "--subs_file"), type = "character", default = NULL,
              help = "Path to the substitutions file", metavar = "character"),
  make_option(c("-t", "--tax"), type = "character", default = NULL,
              help = "Tax name or id. Accepts (unique) tax names in theory, but it's much safer to use tax ids.", 
              metavar = "character"),
  make_option(c("-s", "--stranded"), type="character", default = NULL,
              help = "ds for double stranded or ss for single stranded.",
              metavar = "character"),
  make_option(c("-o", "--output_plot"), type = "character", default = "damage_plot.png",
              help = "Output file for the plot ending in .png or .pdf", metavar = "character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the mandatory arguments are provided
if (is.null(opt$subs_file) || is.null(opt$tax)) {
  # If not, print the help message and stop execution
  print_help(opt_parser)
  stop("Both --subs_file and --tax arguments must be provided", call. = FALSE)
}

if (is.null(opt$stranded) || ! opt$stranded %in% c("ss","ds")){
  print_help(opt_parser)
  stop("Stranded can only be ds or ss.")
}

# Read the file
subs_file <- opt$subs_file
tax_strings <- unlist(strsplit(opt$tax, ","))

# Display the provided arguments
cat("Substitutions file:", subs_file, "\n")
cat("Tax strings:", paste(tax_strings, collapse = ", "), "\n")

# Modify numeric tax strings
modified_tax_strings <- sapply(tax_strings, function(tax) {
  if (grepl("^[0-9]+$", tax)) {
    paste0("^", tax, "\t")
  } else {
    tax
  }
})

# Read the content of the file
file_content <- readLines(subs_file)

# Grep the tax strings
matched_lines <- lapply(modified_tax_strings, function(tax) {
  grep(tax, file_content, value = TRUE)
})

# Combine all matched lines into one vector
matched_lines <- unlist(matched_lines)

# Check if any matched lines exist
if (length(matched_lines) ==0 ) {
  stop("No matching lines found for the given tax strings. Check your subs file, please.", call. = FALSE)
}
if (length(matched_lines) >1 ) {
  stop("More than one matching line found for the given tax strings. Please be more specific.", call. = FALSE)
}

# Now extract the information to make the plot 
extract_value <- function(item) {
  as.numeric(sub(".*:(.*)", "\\1", item))
}

split_line <- strsplit(matched_lines, "\t")[[1]]
title <- paste0("Damage plot for ", split_line[2], " (tax id: ", split_line[1], ")")
data_part <- split_line[3]
data_items <- unlist(strsplit(data_part, " "))
positions <- c(-15:-1, 1:15)
if (opt$stranded == "ds") {
  df <- data.frame(Pos = positions, O = 0, CT = 0, GA = 0)
} else {
  df <- data.frame(Pos = positions, O = 0, CT = 0)
}

for (item in data_items) {
  splititem <- strsplit(item, split = ":")[[1]]
  if (grepl("^O", item)) {
    pos <- as.numeric(gsub("O", "", splititem[1]))
    df$O[df$Pos == pos] <- splititem[2]
  } else if (grepl("^CT", item)) {
    pos <- as.numeric(gsub("CT", "", splititem[1]))
    df$CT[df$Pos == pos] <- splititem[2]
  } else if (grepl("^GA", item) && opt$stranded == "ds") {
    pos <- as.numeric(gsub("GA", "", splititem[1]))
    df$GA[df$Pos == pos] <- splititem[2]
  }
}

# Max y axis value 
max_y <- max(as.numeric(df$O), as.numeric(df$CT), if (opt$stranded == "ds") as.numeric(df$GA) else NA, na.rm = TRUE) * 1.1

# Convert the data frame to long format for ggplot
long_df <- pivot_longer(df, cols = c("O", "CT", if (opt$stranded == "ds") "GA"), names_to = "Type", values_to = "Value")

# Define color palette
color_palette <- c("CT" = "#F8766D", "GA" = "#56B4E9", "O" = "#009E73")

# Function to define x-axis breaks
x_breaks <- function(x) {
  breaks <- seq(min(x), max(x), by = 2)
  if (1 %in% x) breaks <- unique(c(breaks, 1))
  if (-1 %in% x) breaks <- unique(c(breaks, -1))
  breaks
}

# Create two plots and combine them
p1 <- ggplot(subset(long_df, Pos > 0), aes(x = as.factor(Pos), y = as.numeric(Value), color = Type, group = Type)) +
  geom_line(linewidth = 1) +  # Thicker lines using linewidth
  labs(x = "Position", y = "Frequency") +
  scale_color_manual(values = color_palette, labels = c("CT" = "C->T", "GA" = if (opt$stranded == "ds") "G->A" else NULL, "O" = "Other")) +
  scale_x_discrete(breaks = x_breaks(unique(subset(long_df, Pos > 0)$Pos))) +
  ylim(0, max_y) +
  theme_minimal() +
  theme(legend.position = "none")

p2 <- ggplot(subset(long_df, Pos < 0), aes(x = as.factor(Pos), y = as.numeric(Value), color = Type, group = Type)) +
  geom_line(linewidth = 1) +  # Thicker lines using linewidth
  labs(x = "Position", y = "") +
  scale_color_manual(values = color_palette, labels = c("CT" = "C->T", "GA" = if (opt$stranded == "ds") "G->A" else NULL, "O" = "Other")) +
  scale_x_discrete(breaks = x_breaks(unique(subset(long_df, Pos < 0)$Pos))) +
  ylim(0, max_y) +
  theme_minimal() +
  theme(legend.position = "right", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Arrange the plots
combined_plot <- arrangeGrob(
  p1,
  p2,
  ncol = 2
)

# Add the title to the arranged plots
outplot <- arrangeGrob(
  combined_plot,
  top = textGrob(title, gp = gpar(fontsize = 15, fontface = "bold"))
)

# Save the combined plot
ggsave(opt$output_plot, outplot, width = 7.5, height = 3.75)









