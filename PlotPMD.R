#!/usr/bin/env Rscript

# Easy plotting script for PMD score distributions

# Load necessary libraries
library(ggplot2)
library(optparse)

pdf(NULL)

# Define command line options
option_list <- list(
  make_option(c("-o", "--output"), type="character", default="histogram.png",
              help="Output file for histogram [default = %default]", metavar="file"),
  make_option(c("-t", "--title"), type="character", default="PMD scores",
              help="Title prefix for the plot [default = %default]", metavar="title")
)

# Parse command line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Open the connection to standard input
con <- file("stdin")
open(con)

# Read numbers from standard input
numbers <- readLines(con)

# Close the connection explicitly
close(con)

# Convert to numeric and filter out non-finite values and values outside the x-axis range
numbers <- as.numeric(numbers)
filtered_numbers <- numbers[is.finite(numbers) & !is.na(numbers) & numbers > -2 & numbers < 10]

# Determine the plot title
plot_title <- paste(opt$title, "PMD scores")

# Create a data frame from the filtered numbers
data <- data.frame(numbers = filtered_numbers)

# Check if there are any values left after filtering
if (length(filtered_numbers) == 0) {
  stop("No valid PMD scores to plot. Ensure that your input contains values within the range -2 to 10.")
}

# Calculate the maximum count for the histogram
max_count <- max(ggplot_build(ggplot(data, aes(x=numbers)) +
                                geom_histogram(binwidth=0.1))$data[[1]]$count)

# Set y-axis limit 10% higher than the maximum count
y_max <- max_count * 1.1

# Create histogram with white background
p <- ggplot(data, aes(x=numbers)) +
  geom_histogram(aes(y=after_stat(count)), binwidth=0.1, fill="orchid3", alpha=1) +
  labs(title=plot_title,
       x="PMD scores",
       y="Count") +
  scale_x_continuous(breaks=c(-2, 0, 2, 4, 6, 8, 10), limits=c(-2, 10)) +
  theme_bw() +
  ylim(0, y_max)

# Save plot to file
suppressWarnings(ggsave(filename = opt$output, plot = p, width = 8, height = 3))

print(paste("Histogram saved to", opt$output))
