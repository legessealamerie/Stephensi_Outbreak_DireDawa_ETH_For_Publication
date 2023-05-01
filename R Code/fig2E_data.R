################################################################################
# Code to generate data to produce fig 2E
################################################################################

source("cluster_analysis.R")

iconGroups <- expand.grid(
    site = c("City", "University"),
    cluster = 1:4,
    cluster_status = c("High Relatedness", "Medium Relatedness", "Low Relatedness", "Unrelated")
) |>
    dplyr::mutate(
        icon = dplyr::case_when(
            site == "City" ~ "circle",
            site == "University" ~ "triangle"
        ),
        color = cluster_colors[cluster],
        edge_color = dplyr::case_when(
            cluster_status == "High Relatedness" ~ "black",
            TRUE ~ cluster_colors[cluster]
        ),
        cluster_status = as.character(cluster_status),
        entry = dplyr::row_number()
    )

to_map <- clustered_data |>
    dplyr::left_join(epi_clustered_df |> dplyr::select(sample_id, max_out_weight)) |>
    dplyr::left_join(iconGroups) |>
    dplyr::select(sample_id, icon, color, edge_color, cluster_status, cluster)

readr::write_csv(to_map, "processed/cluster_geo_data.csv")
