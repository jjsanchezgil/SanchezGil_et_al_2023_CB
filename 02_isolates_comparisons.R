require(ggVennDiagram, quietly = TRUE)
require(magrittr, quietly = TRUE)
require(ggtext, quietly = TRUE)
require(patchwork, quietly = TRUE)
require(ggplot2)
require(ggpubr)
require(extrafont)
require(pheatmap)
require(ggsci)
require(scales)

minmax <- function(x) {
  return((x - min(x) / (max(x) - min(x))))
}

jaccardIndex <- function(a, b) {
  #Computes Jaccard index (JI) based on vector composition
  #JI = |A AND B|/|A OR B| = |A AND B|/(|A| + |B| - |A AND B|)
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return (intersection/union)
}

jaccardIndexMatrix <- function(mat) {
  jacmat <- matrix(data = numeric(0), nrow = ncol(mat), ncol = ncol(mat))
  for (ro in 2:nrow(jacmat)) {
    for (co in 1:(ro - 1)) {
      roVec <- names(mat[, ro][mat[, ro] != 0])
      coVec <- names(mat[, co][mat[, co] != 0])
      jacmat[ro, co] <- jaccardIndex(roVec, coVec)
    }
  }
  rownames(jacmat) <- colnames(mat)
  colnames(jacmat) <- colnames(mat)
  return(jacmat)
}

jaccardDistanceMatrix <- function(mat) {
  jacmat <- jaccardIndexMatrix(mat)
  jacmat <- as.dist(1 - jacmat)
  return(jacmat)
}

weightedJaccard <- function(A) {
  A <- t(A)
  sim.jac <- matrix(0, nrow=nrow(A), ncol=nrow(A))
  rownames(sim.jac) <- rownames(A)
  colnames(sim.jac) <- rownames(A)

  #weighted jaccard
  pairs <- t(combn(1:nrow(A), 2))
  for (i in 1:nrow(pairs)){
    num <- sum(sapply(1:ncol(A), function(x)(min(A[pairs[i,1],x],A[pairs[i,2],x]))))
    den <- sum(sapply(1:ncol(A), function(x)(max(A[pairs[i,1],x],A[pairs[i,2],x]))))
    # sim.jac[pairs[i,1],pairs[i,2]] <- num/den
    sim.jac[pairs[i,2],pairs[i,1]] <- num/den
  }
  sim.jac[which(is.na(sim.jac))] <- 0
  # diag(sim.jac) <- 1
  return(as.dist(1 - sim.jac))
}

extrafont::font_import("/Windows/Fonts/", pattern = "Montserrat.*", prompt = F)
extrafont::loadfonts(device = "win")

#Read orthogroups table (OrthoFinder output)
og.iso <- read.delim("orthogroups_isolates.tsv", row.names = 1)
ogs <- rownames(og.iso)

##PA table
og.iso.dummy <- og.iso != ""

#List of OGs present in each isolate
ogs.list <- apply(og.iso.dummy, 2, \(x) ogs[x])


#Build Venn diagram
og.venn <- ogs.list %>%
  ggVennDiagram::Venn(.) %>%
  ggVennDiagram::process_data(.)

perc.digits <- 0
og.venn@region %<>%
  dplyr::filter(.data$component == "region") %>%
  dplyr::mutate(percent = paste(round(.data$count*100/sum(.data$count),
                                      digits = perc.digits),"%", sep="")) %>%
  dplyr::mutate(both = paste(.data$count,paste0("(",.data$percent,")"),sep = " "))

hex_codes <- c(pal_material("cyan")(4)[4],
               "grey80",
               pal_material("amber")(9)[4],
               pal_material("deep-orange")(4)[4],
               "grey70",
               "grey80",
               "grey90")

vd.all <- ggplot() +
  geom_sf(aes(fill = name), data = venn_region(og.venn), color = NA) +
  geom_sf(data = venn_setedge(og.venn), show.legend = FALSE, color = "white", linewidth = 2) +
  geom_sf_text(aes(label = name), check_overlap = T, size = 8, parse = T, family = "Montserrat", data = venn_setlabel(og.venn)) +
  geom_sf_text(aes(label = count), position = position_dodge2(2,padding = .2), lineheight = 0,
               size = 4.5, data = venn_region(og.venn)
  ) +
  theme_void() +
  scale_x_continuous(expand = expansion(add = 200)) +
  theme(legend.position = "none",
        strip.text = element_markdown())
vd.all

ragg::agg_jpeg("7isolates_all.jpg", width = 2000, height = 2000, res = 300)
vd.all
invisible(dev.off())

#Heatmap of non-core OGs in isolates ====
og.iso.noncore <- og.iso.dummy[!apply(og.iso.dummy, 1, all),] %>%
  .[apply(., 1, any),]
class(og.iso.noncore) <- "numeric"
og.ph1 <- pheatmap::pheatmap(mat = og.iso.noncore,
                          scale = "none")
or <- og.ph1$tree_row$order
oc <- og.ph1$tree_col$order
og.df <- og.iso.noncore[or, oc] %>% as.data.frame
og.df$og <- rownames(og.df)
og.df %<>% reshape2::melt(., id.var = "og")
og.df$og %<>% factor(., levels = rownames(og.iso.noncore)[rev(or)])
og.df$variable %<>% factor(., levels = colnames(og.iso.noncore)[oc])
og.ph2 <- ggplot(og.df,
       aes(x = variable, y = og, fill = as.factor(value))) +
  geom_tile() +
  scale_fill_manual(values = c("white", "#00BCD4")) +
  coord_flip() +
  theme(legend.position = "none",
        text = element_text(family = "Montserrat", size = 15),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

ragg::agg_jpeg("7isolates_heatmap_noncore_og.jpg", width = 3800, height = 2000, res = 300)
og.ph2
invisible(dev.off())


#Heatmap of all OGs in isolates ====
og.iso.all <- og.iso.dummy[apply(og.iso.dummy, 1, any),]
class(og.iso.all) <- "numeric"
og.phf1 <- pheatmap::pheatmap(mat = og.iso.all,
                             scale = "none")
f.or <- og.phf1$tree_row$order
f.oc <- og.phf1$tree_col$order
ogf.df <- og.iso.all[f.or, f.oc] %>% as.data.frame
ogf.df$og <- rownames(ogf.df)
ogf.df %<>% reshape2::melt(., id.var = "og")
ogf.df$og %<>% factor(., levels = rownames(og.iso.all)[rev(f.or)])
ogf.df$variable %<>% factor(., levels = colnames(og.iso.all)[f.oc])
ogf.ph2 <- ggplot(ogf.df,
                 aes(x = variable, y = og, fill = as.factor(value))) +
  geom_tile() +
  scale_fill_manual(values = c("white", "#00BCD4")) +
  coord_flip() +
  theme(legend.position = "none",
        text = element_text(family = "Montserrat", size = 15),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

ragg::agg_jpeg("7isolates_heatmap_all_og.jpg", width = 3800, height = 2000, res = 300)
ogf.ph2
invisible(dev.off())



#Heatmap of non-core KOs in isolates ====
#Read the KEGG key
kegg.key <- read.delim("KEGG_Pathway_key.txt")
kegg.key$pathway %<>% sub(x = ., "map", "ko")

kkey <- kegg.key$description
names(kkey) <- kegg.key$pathway

kkey[["ko01130"]] <- "Biosynthesis of antibiotics"
kkey[["ko00473"]] <- "D-Amino acid metabolism"
kkey[["ko01053"]] <- "Biosynthesis of siderophore group NRPs"
kkey[["ko02025"]] <- "Biofilm formation"

#Read KOs in isolates
m.paths <- readRDS("paths_matrix.RDS")
ko.iso.dummy <- t(m.paths > 0)
ko.iso.noncore <- ko.iso.dummy[!apply(ko.iso.dummy, 1, all),] %>%
  .[apply(., 1, any),]
class(ko.iso.noncore) <- "numeric"

ko.ph1 <- pheatmap::pheatmap(mat = ko.iso.noncore,
                             scale = "none")
ko.or <- ko.ph1$tree_row$order
ko.oc <- ko.ph1$tree_col$order

ko.df <- ko.iso.noncore[ko.or, ko.oc] %>% as.data.frame
ko.df$ko <- rownames(ko.df)
ko.df %<>% reshape2::melt(., id.var = "ko")
ko.df$ko %<>% factor(., levels = rownames(ko.iso.noncore)[ko.or])
ko.df$variable %<>% factor(., levels = colnames(ko.iso.noncore)[ko.oc])
ko.ph2 <- ggplot(ko.df,
                 aes(x = variable, y = ko, fill = as.factor(value))) +
  geom_tile() +
  scale_fill_manual(values = c("white", "#00BCD4")) +
  coord_flip() +
  scale_y_discrete(labels = kkey) +
  theme(legend.position = "none",
        text = element_text(family = "Montserrat", size = 10),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 5,
                                   angle = 90,
                                   hjust = 0,
                                   vjust = .2))

ragg::agg_jpeg("7isolates_heatmap_noncore_ko.jpg", width = 2000, height = 500, res = 300)
ko.ph2
invisible(dev.off())

#Matrices of genomic distances ----

#Based on ANI ====
ani <- read.delim("isolates.ani.out", header = F)
ani <- ani[ani$V1 != ani$V2, c(1, 2, 3)]
ani$V1 %<>% sub(x = ., ".fasta", "")
ani$V2 %<>% sub(x = ., ".fasta", "")
ani <- tidyr::spread(data = ani, key = V1, value = V3)
rownames(ani) <- ani$V2
ani[is.na(ani)] <- 100
ani <- ani[, -1]
for (i in 2:6) for (j in 1:(i-1)) ani[i, j] <- median(c(ani[i, j], ani[j, i]))
ani[upper.tri(ani, diag = T)] <- NA

ani.df <- reshape2::melt(ani)
ani.df$var2 <- rownames(ani)
ani.df <- ani.df[complete.cases(ani.df),]

g.ani1 <- ggplot(ani.df,
                 aes(x = var2, y = variable, fill = value)) +
  geom_raster() +
  scale_fill_gradientn(colours = pal_material("blue")(9)) +
  geom_text(aes(label = format(ani.df$value, digits = 4))) +
  theme(legend.position = "none",
        text = element_text(family = "Montserrat", size = 12),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = .75))

g.ani1

#Based on OG P/A ====
og.iso.dummy2 <- og.iso.dummy[, colnames(ani)]
og.jac <- jaccardDistanceMatrix(og.iso.dummy2) %>% as.matrix
og.jac[upper.tri(og.jac, diag = T)] <- NA
og.jac <- 100*(1 - og.jac)
og.jac.df <- reshape2::melt(og.jac)
og.jac.df <- og.jac.df[complete.cases(og.jac.df),]

g.og.jac <- ggplot(og.jac.df,
                 aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  scale_fill_gradientn(colours = pal_material("blue")(9)) +
  geom_text(aes(label = format(value, digits = 4))) +
  theme(legend.position = "none",
        text = element_text(family = "Montserrat", size = 12),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = .75))

g.og.jac

#Based on KO P/A ====
ko.iso.dummy

ko.jac <- jaccardDistanceMatrix(ko.iso.dummy) %>% as.matrix
ko.jac[upper.tri(ko.jac, diag = T)] <- NA
ko.jac <- 100*(1 - ko.jac)
ko.jac.df <- reshape2::melt(ko.jac)
ko.jac.df <- ko.jac.df[complete.cases(ko.jac.df),]

g.ko.jac <- ggplot(ko.jac.df,
                   aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  scale_fill_gradientn(colours = pal_material("blue")(9)) +
  geom_text(aes(label = format(value, digits = 4))) +
  theme(legend.position = "none",
        text = element_text(family = "Montserrat", size = 12),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = .75))

g.ko.jac

#Based on KO weighted ====
ko.iso.abundance <- t(m.paths)
ko.wjac <- weightedJaccard(ko.iso.abundance) %>% as.matrix
ko.wjac[upper.tri(ko.wjac, diag = T)] <- NA
ko.wjac <- 100*(1 - ko.wjac)
ko.wjac.df <- reshape2::melt(ko.wjac)
ko.wjac.df <- ko.wjac.df[complete.cases(ko.wjac.df),]

g.ko.wjac <- ggplot(ko.wjac.df,
                   aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  scale_fill_gradientn(colours = pal_material("blue")(9)) +
  geom_text(aes(label = format(value, digits = 4))) +
  theme(legend.position = "none",
        text = element_text(family = "Montserrat", size = 12),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = .75))

g.ko.wjac

#weighted Jaccard on PA is the same as the JI
jaccardDistanceMatrix(t(m.paths > 0)) == weightedJaccard(t(m.paths > 0))

all.df <- data.frame(
  Var1 = og.jac.df$Var1,
  Var2 = og.jac.df$Var2,
  "OG Jaccard" = og.jac.df$value,
  "ANI" = ani.df$value,
  "KO Jaccard" = ko.jac.df$value,
  "KO weighted Jaccard" = ko.wjac.df$value
)

g.comps <- PerformanceAnalytics::chart.Correlation(all.df[, -c(1, 2)],
                                                   histogram = TRUE,
                                                   method = "pearson")

ragg::agg_jpeg("7isolates_heatmap_comps.jpg", width = 2200, height = 2000, res = 300)
((g.ani1 + g.og.jac) / (g.ko.jac + g.ko.wjac))
invisible(dev.off())

ragg::agg_jpeg("7isolates_heatmap_comps_stats.jpg", width = 2000, height = 2000, res = 300)
PerformanceAnalytics::chart.Correlation(all.df[, -c(1, 2)],
                                        histogram = TRUE,
                                        method = "pearson")
invisible(dev.off())


##Figure 2 Venn ----
re <- which(!colnames(og.iso.dummy) %in% c("WCS417", "CHA0"))

#Orthogroups present in any of the rest (dummy)
og.rest.dummy <- rowSums(og.iso.dummy[, re]) %>% as.logical

##Unique CHA0 orthogroups
og.cha.dummy <- og.iso[, "CHA0"] != ""
og.417.dummy <- og.iso[, "WCS417"] != ""

og.shared <- og.cha.dummy & og.417.dummy
og.shared.uniq <- !og.rest.dummy & og.shared

#Total number shared
sum(og.shared.uniq)
#Shared orthogroups
ortho.shared <- ogs[og.shared.uniq]


hex_codes <- c(pal_material("cyan")(4)[4],
               "grey80", 
               pal_material("amber")(9)[4], 
               pal_material("deep-orange")(4)[4], 
               "grey70", 
               "grey80", 
               "grey90")
               
og.venn2 <- list(
  "italic('P. simiae')~WCS417" = ogs[og.417.dummy],
  "italic('P. protegens')~CHA0" = ogs[og.cha.dummy],
  "Rest~of~isolates" = ogs[og.rest.dummy]
) %>%
  ggVennDiagram::Venn(.) %>%
  ggVennDiagram::process_data(., shape_id == "301f")

perc.digits <- 0
og.venn2@region %<>%
  dplyr::filter(.data$component == "region") %>%
  dplyr::mutate(percent = paste(round(.data$count*100/sum(.data$count),
                                      digits = perc.digits),"%", sep="")) %>%
  dplyr::mutate(both = paste(.data$count,paste0("(",.data$percent,")"),sep = " "))

vd1 <- ggplot() + geom_sf(aes(fill = name), data = venn_region(og.venn2), color = NA) +
  scale_fill_manual(values = hex_codes) +
  scale_color_manual(values = hex_codes) +
  geom_sf(data = venn_setedge(og.venn2), show.legend = FALSE, color = "white", linewidth = 2) +
  geom_sf_text(aes(label = name), check_overlap = T, size = 5, parse = T, family = "montserrat", data = venn_setlabel(og.venn2)) +
  geom_sf_text(aes(label = both), lineheight = 0, nudge_y = c(.2, .2, .2, -.2, -.4, -.4, .2),
               # fill = c("white", "white", "white", "white", "white", "white", "white"),
               color = c("grey20", "grey20", "grey20", "black", "grey20", "grey20", "grey20"),
               size = 5, data = venn_region(og.venn2) 
               # label.size = NA, label.padding = unit(4, "mm"),
               # label.r = unit(.2, "lines")) +
  ) +
  theme_void() + 
  scale_x_continuous(expand = expansion(add = 3)) +
  theme(legend.position = "none",
        strip.text = element_markdown())

ragg::agg_jpeg("7isolates_sharedgenes_venn3.jpg", width = 1500, height = 1500, res = 300)
vd1
invisible(dev.off())

##Figure 2 barplot ----
all.annot <- read.delim("all_annotations.tsv")

isolates <- c("AGN" = "CHA0",
              "CHP" = "WCS417",
              "FBB" = "RS158",
              "KMP" = "WCS134",
              "ODA" = "WCS317",
              "JLC" = "WCS358",
              "EEH" = "WCS365")

all.annot$isoCode <- sapply(all.annot$query_name, 
                            stringr::str_trunc, 
                            width = 3, 
                            ellipsis = "") %>% unname
all.annot$Isolates <- isolates[all.annot$isoCode]

all.genes <- rep(rownames(og.iso), lapply(rownames(og.iso), \(x) og.iso[x, ] %>% c %>% unlist %>% unname %>% strsplit(., ",") %>% unlist) %>% sapply(length))
names(all.genes) <- lapply(rownames(og.iso), \(x) og.iso[x, ] %>% c %>% unlist %>% unname %>% strsplit(., ",") %>% unlist) %>% unlist %>% stringr::str_replace_all(.," ", "")
all.annot$OG <- all.genes[all.annot$query_name]
all.annot <- all.annot[, c("Isolates", "OG", "isoCode", "KEGG_Pathway", "KEGG_ko", "query_name")]
all.annot <- all.annot[!is.na(all.annot$OG),]

all.ko <- all.annot[all.annot$KEGG_Pathway != "", "KEGG_Pathway"] %>% 
  sapply(\(x) strsplit(x, ",")) %>% 
  sapply(\(x) x[startsWith(x, "ko")]) %>% 
  unname %>% 
  unlist


iso.names <- unique(all.annot$Isolates)
ortho <- vector(mode = "list", length = length(iso.names))
names(ortho) <- iso.names

all.ogs <- all.annot$OG %>% unique %>% sort
m.ogs <- matrix(0, nrow = length(iso.names), ncol = length(all.ogs), 
                dimnames = list(iso.names, all.ogs))

for (iso in iso.names) {
  m.ogs[iso, unique(all.annot[all.annot$Isolates == iso, "OG"])] <- 1
}

ogs.417 <- m.ogs["WCS417", ] == 1
ogs.cha <- m.ogs["CHA0", ] == 1
ogs.rest <- (m.ogs[!rownames(m.ogs) %in% c("WCS417", "CHA0"), ] == 1) %>% apply(., 2, any) 

sum(!ogs.rest & (ogs.417 & ogs.cha))
ogs.shared <- colnames(m.ogs)[!ogs.rest & (ogs.417 & ogs.cha)]

all.annot$KEGG_Pathway %<>% sapply(\(x) sub(x = x, "ko02026", "ko02025")) %>% sapply(\(x) sub(x = x, "ko05111", "ko02025"))

paths <- vector(mode = "list", length = length(iso.names))
names(paths) <- iso.names
for (iso in unique(all.annot$Isolates)) {
  paths[[iso]] <- all.annot[all.annot$KEGG_Pathway != "" &
                              all.annot$Isolates == iso, "KEGG_Pathway"] %>%
    sapply(\(x) strsplit(x, ",")) %>% 
    sapply(\(x) x[startsWith(x, "ko")]) %>% 
    unname %>% 
    unlist
}
all.paths <- paths %>% unlist %>% unname %>% sort %>% unique
m.paths <- matrix(0, nrow = length(iso.names), ncol = length(all.paths), 
                  dimnames = list(iso.names, all.paths))
for (iso in iso.names) {
  iso.table <- table(paths[[iso]])
  m.paths[iso, names(iso.table)] <- unname(iso.table)
}

annot.shared <- all.annot[all.annot$OG %in% ogs.shared,]
shared.paths <- vector(mode = "list", length = length(iso.names))
names(shared.paths) <- iso.names
for (iso in unique(all.annot$Isolates)) {
  shared.paths[[iso]] <- annot.shared[annot.shared$KEGG_Pathway != "" &
                                        annot.shared$Isolates == iso, "KEGG_Pathway"] %>%
    sapply(\(x) strsplit(x, ",")) %>% 
    sapply(\(x) x[startsWith(x, "ko")]) %>% 
    unname %>% 
    unlist
}
all.shared.paths <- shared.paths %>% unlist %>% unname %>% sort %>% unique
m.shared.paths <- matrix(0, nrow = length(iso.names), ncol = length(all.shared.paths), 
                         dimnames = list(iso.names, all.shared.paths))
for (iso in iso.names) {
  iso.table <- table(shared.paths[[iso]])
  m.shared.paths[iso, names(iso.table)] <- unname(iso.table)
}
m.shared.paths <- m.shared.paths[c("WCS417", "CHA0"),]
m.paths <- m.paths[!rownames(m.paths) %in% c("WCS417", "CHA0"), colnames(m.shared.paths)]

fisher.417 <- vector(mode = "list", length = ncol(m.shared.paths))
names(fisher.417) <- colnames(m.shared.paths)
for (i in seq_along(colnames(m.shared.paths))) {
  path <- colnames(m.shared.paths)[i]
  fisher.417[[path]] <- data.frame(path = c(m.shared.paths["WCS417", i], 
                                            sum(m.paths[, i])), 
                                   other = c(sum(m.shared.paths["WCS417", -i]), 
                                             sum(m.paths[, -i])), 
                                   row.names = c("in", "out")) %>% 
    fisher.test(alternative = "greater") %>% 
    '['(c("p.value", "estimate"))
}

fisher.cha <- vector(mode = "list", length = ncol(m.shared.paths))
names(fisher.cha) <- colnames(m.shared.paths)
for (i in seq_along(colnames(m.shared.paths))) {
  path <- colnames(m.shared.paths)[i]
  fisher.cha[[path]] <- data.frame(path = c(m.shared.paths["CHA0", i], 
                                            sum(m.paths[, i])), 
                                   other = c(sum(m.shared.paths["CHA0", -i]), 
                                             sum(m.paths[, -i])), 
                                   row.names = c("in", "out")) %>% 
    fisher.test(alternative = "greater") %>% 
    '['(c("p.value", "estimate"))
}

kegg.key <- read.delim("KEGG_Pathway_key.txt")
kegg.key$pathway %<>% sub(x = ., "map", "ko")

kkey <- kegg.key$description 
names(kkey) <- kegg.key$pathway 

kkey[["ko01130"]] <- "Biosynthesis of antibiotics"
kkey[["ko00473"]] <- "D-Amino acid metabolism"
kkey[["ko01053"]] <- "Biosynthesis of siderophore group NRPs"
kkey[["ko02025"]] <- "Biofilm formation"


enrich.df <- data.frame(rbind(do.call(rbind, fisher.417),
                              do.call(rbind, fisher.cha)),
                        isolate = c(rep("WCS417", length(fisher.417)),
                                    rep("CHA0", length(fisher.cha))),
                        pathway = c(names(fisher.417), names(fisher.cha)),
                        row.names = NULL)
enrich.df$p.value %<>% as.numeric()
enrich.df$estimate %<>% as.numeric()
enrich.df$adjusted <- p.adjust(enrich.df$p.value, "fdr")
enrich.df$pathway <- factor(enrich.df$pathway, levels = enrich.df[order(enrich.df$p.value, decreasing = T), "pathway"] %>% unique)
enrich.df$isolate <- factor(enrich.df$isolate, levels = c("WCS417", "CHA0"))

hex_codes.bar <- c(pal_material("amber")(4)[4],
                   pal_material("cyan")(9)[4])

content <- data.frame(counts = c(m.shared.paths[rownames(m.shared.paths) == "WCS417", ], 
                                 m.shared.paths[rownames(m.shared.paths) == "CHA0", ]), 
                      isolate = rep(c("WCS417", "CHA0"), each = ncol(m.shared.paths)),
                      path = colnames(m.shared.paths))
content$path <- factor(content$path, levels = content$path[order(enrich.df$p.value, decreasing = T)] %>% unique)
content$isolate <- factor(content$isolate, levels = c("WCS417", "CHA0"))
content$path <- factor(content$path, levels = content[order(content$counts, decreasing = F), "path"] %>% unique)

paths.more.than.one <- content %>%
  dplyr::group_by(path) %>%
  dplyr::mutate(maxc = max(counts)) %>%
  dplyr::filter(maxc > 1) %>%
  dplyr::filter(!path %in% c("ko01100", "ko01120")) %>%
  as.data.frame %>%
  .[order(.$maxc, decreasing = T),] %>%
  '[['("path") %>% unique



cn1 <- ggplot(content[content$path %in% paths.more.than.one, ], 
              aes(x = counts, y = path, fill = isolate)) + 
  geom_col(position = "dodge", width = .8) +
  scale_y_discrete(labels = kkey) +
  scale_x_continuous(position = "top", breaks = seq(0, 12, 2)) +
  scale_fill_manual(values = hex_codes.bar, labels = c("*P. simiae* WCS417", "*P. protegens* CHA0")) +
  xlab("No. annotations") +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_blank(),
        text = element_text(family = "Montserrat", size = 12),
        axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray80")
  ) 
cn1


enrich.df$pathway <- factor(enrich.df$pathway, levels = enrich.df[order(content$counts, decreasing = T), "pathway"] %>% unique)
en1 <- ggplot(enrich.df[!enrich.df$pathway %in% c("ko01100", "ko01120") & 
                          enrich.df$pathway %in% paths.more.than.one, ],
              aes(y = pathway, x = -log10(p.value), fill = isolate)) + 
  geom_col(position = "dodge", width = .8) +
  geom_vline(xintercept = -log10(.05), lty = "dashed", color = "gray60") +
  geom_vline(xintercept = -log10(.001), lty = "dashed", color = "gray60") +
  scale_y_discrete(labels = kkey) +
  scale_x_continuous(position = "top", breaks = seq(0, 10, 2)) +
  scale_fill_manual(values = hex_codes.bar,
                    labels = c("*P. simiae* WCS417", "*P. protegens* CHA0")) +
  xlab("-log<sub>10</sub>*p*-value") +
  theme(axis.title.x.top = element_markdown(),
        axis.title.y = element_blank(),
        text = element_text(family = "Montserrat", size = 12),
        axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray80")
  ) 
en1


#HORIZONTAL ==== 
content <- data.frame(counts = c(m.shared.paths[rownames(m.shared.paths) == "WCS417", ], 
                                 m.shared.paths[rownames(m.shared.paths) == "CHA0", ]), 
                      isolate = rep(c("WCS417", "CHA0"), each = ncol(m.shared.paths)),
                      path = colnames(m.shared.paths)) %>%
  dplyr::filter(!path %in% c("ko01100", "ko01120"))
content$isolate <- factor(content$isolate, levels = c("WCS417", "CHA0"))
content$path <- factor(content$path, levels = content[order(content$counts, decreasing = T), "path"] %>% unique)
paths.more.than.one <- content %>%
  dplyr::group_by(path) %>%
  dplyr::mutate(maxc = max(counts)) %>%
  dplyr::filter(maxc > 1) %>%
  dplyr::mutate(path = as.character(path)) %>%
  dplyr::filter(!path %in% c("ko01100", "ko01120")) %>%
  as.data.frame %>%
  .[order(.$maxc, decreasing = T),] %>%
  '[['("path") %>% unique
content <- content[content$path %in% paths.more.than.one, ]
cn2 <- ggplot(content, 
              aes(x = path, y = counts, fill = isolate)) + 
  geom_col(position = "dodge", width = .8) +
  scale_x_discrete(labels = kkey) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  scale_fill_manual(values = hex_codes.bar, labels = c("*P. simiae* WCS417", "*P. protegens* CHA0")) +
  ylab("No. annotations") +
  theme(axis.title.y = element_markdown(),
        axis.title.x = element_blank(),
        text = element_text(family = "Montserrat", size = 12),
        axis.text.y = element_markdown(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(),
        legend.title = element_blank(),
        legend.position = c(.8, .5),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray80")
  ) 
cn2

paths.more.than.one.e <- content %>%
  dplyr::group_by(path) %>%
  dplyr::mutate(maxc = max(counts)) %>%
  dplyr::filter(maxc > 1) %>%
  '[['("path")
# ce <- as.character(enrich.df$pathway)
# cc <- as.character(content$path)
# cm <- ce[match(ce[ce %in% cc], cc)] %>% unique
# enrich.df <- enrich.df[enrich.df$pathway %in% cm, ]
# enrich.df$pathway %<>% factor(., levels = cm)
enrich.df <- data.frame(rbind(do.call(rbind, fisher.417),
                              do.call(rbind, fisher.cha)),
                        isolate = c(rep("WCS417", length(fisher.417)),
                                    rep("CHA0", length(fisher.cha))),
                        pathway = c(names(fisher.417), names(fisher.cha)),
                        row.names = NULL)
enrich.df <- enrich.df[enrich.df$pathway %in% as.character(content$path),]
enrich.df$p.value %<>% as.numeric()
enrich.df$estimate %<>% as.numeric()
enrich.df$adjusted <- p.adjust(enrich.df$p.value, "fdr")
enrich.df$pathway <- factor(enrich.df$pathway, levels = levels(content$path) %>% unique)
enrich.df$isolate <- factor(enrich.df$isolate, levels = c("WCS417", "CHA0"))

en2 <- ggplot(enrich.df,
              aes(x = pathway, y = -log10(p.value), fill = isolate)) + 
  geom_col(position = "dodge", width = .8) +
  geom_hline(yintercept = -log10(.05), lty = "dashed", color = "gray60") +
  geom_text(aes(x = 25, y = -log10(.05) + .5, label = "italic(p)==.05"),
            family = "Montserrat",
            parse = TRUE,
            fontface = "italic") +
  geom_hline(yintercept = -log10(.001), lty = "dashed", color = "gray60") +
  geom_text(aes(x = 25, y = -log10(.001) + .5, label = "italic(p)==.001"),
            family = "Montserrat",
            parse = TRUE,
            fontface = "italic") +
  scale_x_discrete(labels = kkey) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  scale_fill_manual(values = hex_codes.bar,
                    labels = c("*P. simiae* WCS417", "*P. protegens* CHA0")) +
  ylab("-log<sub>10</sub>*p*-value") +
  theme(axis.title.x.top = element_markdown(),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(),
        text = element_text(family = "Montserrat", size = 12),
        axis.text.x = element_text(hjust = 1, angle = 30, size = 6),
        legend.text = element_markdown(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray80")
  ) 

en2
layout.summary2 <- c("
AAAAAA
BBBBBB
")
min.margin <- kkey[levels(content$path)[1]] %>% nchar
fig2_b <- cn2 / en2 + theme(plot.margin = margin(l = min.margin)) + plot_layout(design = layout.summary2)
fig2_b











