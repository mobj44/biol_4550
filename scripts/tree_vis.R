library(ape)
library(phytools)
library(ggtree)
library(ggplot2)

# load tree file
tree <- read.tree("./trees/ml_trees/iqtree_supermatrix.treefile")

# set outgroup
outgroup_tips <- grep("ursus", tree$tip.label, value = TRUE)

# root tree
rooted_tree <- root(tree,
                    outgroup = outgroup_tips,
                    resolve.root = TRUE)

root_node <- Ntip(rooted_tree) + 1

# time tree
cal <- makeChronosCalib(
    rooted_tree,
    node    = root_node,
    age.min = 45,
    age.max = 45
)


tree_time <- chronos(
    rooted_tree,
    cal   = cal,
    model = "correlated"
)


root_age <- max(node.depth.edgelength(tree_time))


breaks <- seq(0, round(root_age), by = 2)

# plot tree
p <- ggtree(tree_time, layout = "rectangular") +
    geom_tree(size = 0.5) +

    geom_tiplab(size = 3, hjust = 0) +

    geom_label2(
        aes(
            subset = !isTip & as.numeric(label) >= 70,
            label  = label
        ),
        size       = 2.4,
        fill       = "lightgreen",
        label.size = 0,
        nudge_x    = 0.1,
        vjust      = -0.5
    ) +

    scale_x_continuous(
        name   = "Time before present (Mya)",
        limits = c(0, root_age * 1.05),
        breaks = breaks,
        labels = round(root_age - breaks, 1)
    ) +

    coord_cartesian(clip = "off") +
    theme_tree2() +
    theme(
        plot.margin  = margin(10, 150, 10, 10),
        axis.text.x  = element_text(size = 10),
        axis.title.x = element_text(size = 12)
    )

p <- p +
    geom_label2(
        aes(
            subset = !isTip,
            label  = label,
            fill   = cut(
                as.numeric(label),
                breaks = c(-Inf, 70, 80, Inf),
                labels = c("low", "medium", "high")
            )
        ),
        size       = 2.4,
        alpha      = 0.7,
        label.size = 0,
        nudge_x    = 0.1,
        vjust      = -0.5,
        color      = "black"
    ) +
    scale_fill_manual(
        values = c(
            low    = "coral",     # < 70
            medium = "cornflowerblue",  # 70–79
            high   = "lightgreen"    # 80+
        ),
        guide = "none"
    )


# save tree
ggsave("final_tree.png", p,
       width = 12, height = 18, dpi = 300)
