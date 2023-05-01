sequence_to_alleles <- function(df) {
  df |>
    dplyr::select(locus, asv) |>
    dplyr::distinct() |>
    dplyr::group_by(locus) |>
    dplyr::mutate(allele = seq(1, dplyr::n())) |>
    dplyr::ungroup() |>
    dplyr::mutate(allele = paste0(locus, ".", allele)) |>
    dplyr::rename(sequence = asv)
}

dna_meta <- readr::read_rds("data/DNA_meta.rds")
dd_epi <- readr::read_rds("data/dd_epi.rds")
allele_data <- readr::read_rds("data/allele_data.rds")


#### Load and join sequencing data

allele_sequences <- sequence_to_alleles(allele_data)

# remove controls, combine data based on ASV, rename samples
processed_dat <- allele_data |>
  dplyr::filter(!grepl("Control", sampleID)) |>
  dplyr::left_join(allele_sequences |> dplyr::select(-locus), by = c("asv" = "sequence")) |>
  dplyr::group_by(sampleID, locus, allele) |>
  dplyr::mutate(norm.reads.allele = reads / sum(reads)) |>
  dplyr::group_by(sampleID, locus) |>
  dplyr::mutate(norm.reads.locus = reads / sum(reads)) |>
  dplyr::mutate(n.alleles = dplyr::n()) |>
  dplyr::ungroup() |>
  dplyr::rename(sample_id = sampleID) |>
  dplyr::mutate(
    barcode =
      purrr::map_chr(
        stringr::str_split(sample_id, "-"),
        ~ stringr::str_split(.x[[8]], "_", simplify = TRUE)[[1]]
      )
  ) |>
  dplyr::left_join(dna_meta |> dplyr::select(individual, barcode)) |>
  dplyr::distinct() |>
  dplyr::mutate(individual = dplyr::coalesce(individual, barcode)) |>
  dplyr::select(-sample_id, -barcode) |>
  dplyr::rename(sample_id = individual) |>
  dplyr::group_by(sample_id, locus, allele, asv) |>
  dplyr::summarise(reads = sum(reads)) |>
  dplyr::ungroup() |>
  dplyr::filter(grepl("-1A$", locus))

dir.create("processed", showWarnings = F)
dir.create("figures", showWarnings = F)

readr::write_rds(processed_dat, file = "processed/allele_dat.rds")
