get_sample_coverage_by_locus <- function(dat) {
  dat |>
    dplyr::select(sample_id, locus) |>
    dplyr::distinct() |>
    dplyr::group_by(locus) |>
    dplyr::summarize(num_samples = dplyr::n())
}

get_locus_coverage_by_sample <- function(dat) {
  dat |>
    dplyr::select(sample_id, locus) |>
    dplyr::distinct() |>
    dplyr::group_by(sample_id) |>
    dplyr::summarize(num_loci = dplyr::n())
}

dat <- readr::read_rds("processed/allele_dat.rds")

# min. number of reads for an allele to be called
min_reads <- 20
# min. number of samples that a locus has data for
min_sample_coverage <- 100
# min. number of loci that a sample has data for
min_locus_coverage <- 50

# remove alleles that do not have at least 20 reads
filtered_dat <- dat |>
  dplyr::filter(reads >= min_reads)

# identify loci with data for at least min_sample_coverage samples
valid_loci <- get_sample_coverage_by_locus(filtered_dat) |>
  dplyr::filter(num_samples > min_sample_coverage) |>
  dplyr::pull(locus)

# remove loci without enough coverage
filtered_dat <- filtered_dat |>
  dplyr::filter(locus %in% valid_loci)

# identify samples with data for at least min_locus_coverage loci
valid_samples <- get_locus_coverage_by_sample(filtered_dat) |>
  dplyr::filter(num_loci > min_locus_coverage) |>
  dplyr::pull(sample_id)

# remove samples that do not have enough coverage
filtered_dat <- filtered_dat |>
  dplyr::filter(sample_id %in% valid_samples)

# remove all uninformative loci, i.e. loci with only 1 allele present
uninformative_loci <- filtered_dat |>
  dplyr::select(locus, allele) |>
  dplyr::distinct() |>
  dplyr::group_by(locus) |>
  dplyr::summarise(n = dplyr::n()) |>
  dplyr::filter(n == 1) |>
  dplyr::pull(locus)

filtered_dat <- filtered_dat |>
  dplyr::filter(!(locus %in% uninformative_loci))

readr::write_rds(filtered_dat, file = "processed/filtered_dat.rds")
