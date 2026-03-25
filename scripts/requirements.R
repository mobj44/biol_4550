required_packages <- c("ape", "phangorn", "phytools", "ggtree", "ggplot2")

to_install <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
to_install
if(length(to_install) > 0) {
    install.packages(to_install)
}
