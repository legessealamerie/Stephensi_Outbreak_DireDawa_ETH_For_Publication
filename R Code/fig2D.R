################################################################################
# Code to generate fig 2D
################################################################################

source("cluster_analysis.R")

library(patchwork)



center_offset <- 0.25
to_plot <- clustered_data |>
    dplyr::group_by(cluster, test_date) |>
    dplyr::mutate(
        y.idx = seq(
            unique(cluster) - (n() - 1) * (center_offset / 8),
            unique(cluster) + (n() - 1) * (center_offset / 8),
            length.out = n()
        )
    ) |>
    dplyr::filter(!is.na(y.idx)) |>
    dplyr::ungroup()

# Generate cluster plots
relatedness_colors <- rev(RColorBrewer::brewer.pal(4, "BuPu"))
cluster_plots <- list()
for (cluster_num in unique(graph_tidy %N>% dplyr::pull(cluster))) {
    cluster_subset <- graph_tidy %N>%
        dplyr::filter(cluster == cluster_num)
    cluster_plots[[cluster_num]] <- ggraph::ggraph(
        cluster_subset,
        layout = "fr"
    ) +
        ggraph::geom_edge_link(
            aes(color = forcats::fct_reorder(relatedness, relatedness_order)),
            edge_width = 2
        ) +
        scale_edge_alpha(guide = "none") +
        scale_edge_width(guide = "none") +
        scale_edge_color_manual(
            values = relatedness_colors,
            breaks = c(
                "High Relatedness", "Medium Relatedness",
                "Low Relatedness", "Unrelated"
            ),
            labels = c(
                "High", "Medium",
                "Low", "Unrelated"
            ),
            drop = FALSE
        ) +
        labs(edge_color = "Pairwise Relatedness") +
        ggraph::geom_node_point(
            size = 3,
            aes(shape = site.y, color = as.factor(cluster))
        ) +
        scale_color_manual(
            breaks = c(1, 2, 3, 4),
            values = cluster_colors,
            guide = "none"
        ) +
        scale_shape(guide = "none") +
        theme_void()
}


# Generate timeline plot for case incidence
timeline_plt_site <- ggplot(
    to_plot,
) +
    geom_rect(
        aes(
            xmin = min_date, xmax = max_date,
            ymin = cluster - center_offset, ymax = cluster + center_offset,
        ),
        fill = cluster_annotations$cluster_color,
        alpha = .2,
        data = cluster_annotations
    ) +
    geom_point(
        aes(x = test_date, y = y.idx, shape = site),
        size = 5,
        data = to_plot |> dplyr::filter(cluster_status == "High Relatedness"),
    ) +
    geom_point(
        aes(x = test_date, y = y.idx, shape = site),
        color = to_plot$cluster_color,
        size = 3
    ) +
    labs(
        x = "Test Date", y = element_blank(),
        shape = "Collection Site",
    ) +
    theme_minimal() +
    theme(
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 25),
        axis.title.x = element_text(size = 27),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(linetype = "dashed", color = "grey")
    )


# Generate marginal incidence timeline
date_marginal_plt <- ggplot(
    to_plot, aes(x = test_date)
) +
    geom_histogram(bins = 100) +
    xlim(min_date, max_date) +
    theme_void()


layout <- "
##EEEEEE
GDFFFFFF
GCFFFFFF
GBFFFFFF
GAFFFFFF
"

lineage_annotations <- ggplot(cluster_annotations, aes(y = cluster, x = 0)) +
    geom_text(
        aes(label = paste("Lineage", cluster, sep = " ")),
        size = 10,
        color = "black",
    ) +
    lims(y = c(.75, 4.25)) +
    theme_void()


# Construct fig 2D
cluster_timeline_plt <- (
    cluster_plots[[1]] +
        cluster_plots[[2]] +
        cluster_plots[[3]] +
        cluster_plots[[4]] +
        patchwork::plot_layout(tag_level = "new")) +
    date_marginal_plt +
    timeline_plt_site +
    lineage_annotations +
    patchwork::plot_layout(design = layout, tag_level = "new", guides = "collect") &
    theme(
        text = element_text(size = 30, family = "Calibri"),
    )

ggsave("figures/fig2D.png", cluster_timeline_plt, width = 21, height = 10)
