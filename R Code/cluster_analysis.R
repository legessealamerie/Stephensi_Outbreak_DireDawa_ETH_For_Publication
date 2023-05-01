library(tidygraph)
library(ggplot2)
library(ggraph)


matrix_to_longformat <- function(mat) {
    as.data.frame(mat) |>
        tibble::rownames_to_column("sample.x") |>
        tidyr::pivot_longer(
            -sample.x,
            names_to = "sample.y", values_to = "value"
        )
}

dd_epi <- readr::read_rds("data/dd_epi.rds")
filtered_dat <- readr::read_rds("processed/filtered_dat.rds")
dres <- readr::read_rds("processed/dcifer_results.rds")


mydsmp <- dcifer::formatDat(filtered_dat, "sample_id", "locus", "allele")
coi <- dcifer::getCOI(mydsmp)
coi_data <- data.frame(
    sample_coi = coi,
    sample_id = sapply(stringr::str_split(names(coi), "\\."), function(x) x[1])
)


dmat <- dres[, , "estimate"]
dmat_full <- dmat
dmat_full[upper.tri(dmat)] <- t(dmat)[upper.tri(t(dmat))]

dmat_CI <- dres[, , "CI_lower"]
dmat_CI_full <- dmat_CI
dmat_CI_full[upper.tri(dmat_CI)] <- t(dmat_CI)[upper.tri(t(dmat_CI))]


samples_to_keep <- colnames(dmat_CI)

# get pairwise relatedness where lower CI > 0.1
ci_edge_df <- matrix_to_longformat(dmat_CI) |>
    dplyr::filter(!is.na(value), value > 0.1) |>
    dplyr::rename(from = sample.x, to = sample.y, lower = value)

rel_edge_df <- matrix_to_longformat(dmat) |>
    dplyr::rename(from = sample.x, to = sample.y, weight = value)

edge_df <- ci_edge_df |>
    dplyr::left_join(rel_edge_df, by = c("from", "to")) |>
    dplyr::select(-lower)

edge_df_mirror <- edge_df |> dplyr::rename(to = from, from = to)

# classify relatedness into bins
full_edge_df <- rbind(edge_df, edge_df_mirror) |>
    dplyr::mutate(relatedness = dplyr::case_when(
        weight > .9 ~ "High Relatedness",
        weight > .5 ~ "Medium Relatedness",
        weight > .2 ~ "Low Relatedness",
        TRUE ~ "Unrelated"
    )) |>
    dplyr::mutate(relatedness_order = dplyr::case_when(
        weight > .9 ~ 1,
        weight > .5 ~ 2,
        weight > .2 ~ 3,
        TRUE ~ 4
    ))

out_weight_summary <- full_edge_df |>
    dplyr::group_by(from) |>
    dplyr::summarise(
        max_out_weight = max(weight[weight > 0])
    ) |>
    dplyr::rename(sample_id = from)

cluster_status_colors <- RColorBrewer::brewer.pal(4, "Set2")

epi_clustered_df <- as.data.frame(dd_epi) |>
    dplyr::relocate(sample_id) |>
    dplyr::left_join(out_weight_summary) |>
    dplyr::filter(sample_id %in% samples_to_keep) |>
    dplyr::mutate(
        cluster_status = dplyr::case_when(
            max_out_weight > .9 ~ "High Relatedness",
            max_out_weight > .5 ~ "Medium Relatedness",
            max_out_weight > .2 ~ "Low Relatedness",
            TRUE ~ "Unrelated"
        ),
        cluster_status_order = dplyr::case_when(
            max_out_weight > .9 ~ 1,
            max_out_weight > .5 ~ 2,
            max_out_weight > .25 ~ 3,
            TRUE ~ 4
        ),
        cluster_status_color = dplyr::case_when(
            max_out_weight > .9 ~ cluster_status_colors[1],
            max_out_weight > .5 ~ cluster_status_colors[2],
            max_out_weight > .25 ~ cluster_status_colors[3],
            TRUE ~ cluster_status_colors[4]
        )
    )


graph <- igraph::simplify(
    igraph::graph_from_data_frame(
        full_edge_df,
        directed = FALSE, vertices = epi_clustered_df
    )
)

graph <- igraph::set_edge_attr(
    graph, "relatedness",
    value = dplyr::case_when(
        (igraph::E(graph)$weight / 2) > .9 ~ "High Relatedness",
        (igraph::E(graph)$weight / 2) > .5 ~ "Medium Relatedness",
        (igraph::E(graph)$weight / 2) > .2 ~ "Low Relatedness",
        TRUE ~ "Unrelated"
    )
)

graph <- igraph::set_edge_attr(
    graph, "relatedness_order",
    value = dplyr::case_when(
        (igraph::E(graph)$weight / 2) > .9 ~ 1,
        (igraph::E(graph)$weight / 2) > .5 ~ 2,
        (igraph::E(graph)$weight / 2) > .2 ~ 3,
        TRUE ~ 4
    )
)

graph <- igraph::set_edge_attr(
    graph, "relatedness",
    value = as.factor(igraph::E(graph)$relatedness)
)


clusters <- igraph::cluster_fast_greedy(graph)
graph_tidy <- tidygraph::as_tbl_graph(graph) %N>%
    dplyr::mutate(cluster = clusters$membership)



clustered_data <- data.frame(
    name = tidygraph::pull(graph_tidy, "name"),
    cluster = tidygraph::pull(graph_tidy, "cluster"),
    cluster_status = tidygraph::pull(graph_tidy, "cluster_status"),
    cluster_status_color = tidygraph::pull(graph_tidy, "cluster_status_color"),
    cluster_status_order = tidygraph::pull(graph_tidy, "cluster_status_order")
) |>
    dplyr::left_join(dd_epi, by = c("name" = "sample_id")) |>
    dplyr::rename(sample_id = name) |>
    dplyr::filter(
        !is.na(casecat.y),
        !is.na(test_date),
    ) |>
    dplyr::left_join(
        epi_clustered_df |> dplyr::select(
            sample_id,
            max_out_weight,
            cluster_status,
            cluster_status_color,
        ),
    )



min_date <- min(clustered_data$test_date, na.rm = TRUE) - 5
max_date <- max(clustered_data$test_date, na.rm = TRUE) + 5
cluster_colors <- RColorBrewer::brewer.pal(4, "Dark2")

cluster_annotations <- data.frame(
    cluster = c(1, 2, 3, 4),
    cluster_color = cluster_colors,
    min_date = min_date,
    max_date = max_date
)


clustered_data <- clustered_data |>
    dplyr::left_join(
        cluster_annotations,
        by = "cluster"
    ) |>
    dplyr::mutate(
        site = dplyr::case_when(
            site.y == "DD city" ~ "City",
            site.y == "DD University" ~ "University"
        )
    )
