library(ggsci)
library(readxl)
library(magrittr)
library(pheatmap)
library(ggplot2)
library(network)
library(ggnet)
library(scales)
library(ggtext)

minmax <- function(vec, nmin = 0, nmax = 1) {
  vmin <- min(vec)
  vmax <- max(vec)
  scaled_vec <- nmin + (vec - vmin)*(nmax - nmin)/(vmax - vmin)
  return(scaled_vec)
}

scaleWeight <- function(vec, alpha = 1, nmin = 0) {
  nvec <- alpha*abs(scale(abs(vec)))
  nvec[nvec < nmin] <- nmin
  return(nvec)
}

# scales::show_col(ggsci::pal_material(palett = "teal")(15))
pal <- c(
  ggsci::pal_material(palette = "red")(10)[10],
  ggsci::pal_material(palette = "red")(10)[5],
  ggsci::pal_material(palette = "green")(10)[5],
  ggsci::pal_material(palette = "blue-grey")(10)[2],
  ggsci::pal_material(palette = "purple")(10)[2],
  ggsci::pal_material(palette = "teal")(10)[2],
  ggsci::pal_material(palette = "orange")(10)[2]
)

iols_genes <- c(
  "PS417_11845", "PS417_11850", "PS417_11855", "PS417_11860",
  "PS417_11865", "PS417_11870", "PS417_11875", "PS417_11880",
  "PS417_11885", "PS417_11890", "PS417_11895", "PS417_11900"
)
iols_genenames <- c(
  "iolR", "iolC", "iolE", "iolB",
  "iolI", "iolD", "iolG", "mocA",
  "iatA", "iatC", "iatB", "dksA"
)

t_score_min <- 1.2
fsRaw <- readxl::read_excel("WCS417_fitness_scores.xlsx")
fsRaw[fsRaw$sysName %in% iols_genes, "geneName"] <- iols_genenames
#Check with fsRaw[fsRaw$sysName %in% iols, 1:5]

tsRaw <- readxl::read_excel("WCS417_t_scores.xlsx")
tsRaw[tsRaw$sysName %in% iols_genes, "geneName"] <- iols_genenames
#Check with tsRaw[tsRaw$sysName %in% iols, 1:5]

rownames(tsRaw) <- rownames(fsRaw)

annot <- read.delim("WCS417_meta.txt")
annot <-  annot[, c("expName", "expDesc", "expGroup", "condition_2", "units_1", "concentration_1")]
colnames(annot) <- c("Experiment", "Description", "Group", "Condition", "Units", "Concentration")
rownames(annot) <- annot$Experiment

ino_exp <- annot[annot$expDesc == "m-Inositol (C)", "expName"]

fs <- fsRaw[, 6:ncol(fsRaw)] %>% as.matrix()
rownames(fs) <- fsRaw$sysName
rownames(fs)[which(rownames(fs) %in% iols_genes)] <- iols_genenames
colnames(fs) <- strsplit(colnames(fs), " ") %>% sapply(., '[', 1)

ts <- tsRaw[, 6:ncol(tsRaw)] %>% as.matrix()
rownames(ts) <- tsRaw$sysName
rownames(ts)[which(rownames(ts) %in% iols_genes)] <- iols_genenames
colnames(ts) <- strsplit(colnames(fs), " ") %>% sapply(., '[', 1)

#Remove cluster of genes which fitness spans all conditions
cls_fs <- pheatmap(fs, silent = T)
cls_fs <- cutree(cls_fs$tree_row, k = 2)
fs <- fs[!as.logical(cls_fs - 1), ]
ts <- ts[!as.logical(cls_fs - 1), ]

fs <- t(fs) %>% as.data.frame()
fs <- cbind(annot[rownames(fs), c("Group", "Condition", "Units", "Concentration")], fs)
fs <- fs %>% dplyr::group_by(Group, Condition, Units, Concentration) %>%
  dplyr::summarise(dplyr::across(.fns = median))
fs$Concentration <- sub("/.0+$", "", fs$Concentration) %>% as.numeric()
fs$Name <- fs$Condition
fs <- as.data.frame(fs)
fs$NameGroup <- paste(fs$Name, fs$Group)
fs[fs$NameGroup %in% names(table(fs$NameGroup))[table(fs$NameGroup) != 1], "Name"] <- with(fs[fs$NameGroup %in% names(table(fs$NameGroup))[table(fs$NameGroup) != 1], ], paste(Condition, sub("/.0+$", "", Concentration), sub("vol%", "%", Units)))
fs.meta <- fs[, c("Group", "Condition", "Units", "Concentration", "Name", "NameGroup")]
fs <- fs[, !colnames(fs) %in% c("Group", "Condition", "Units", "Concentration", "Name", "NameGroup")] %>% as.matrix()

ts <- t(ts) %>% as.data.frame()
ts <- cbind(annot[rownames(ts), c("Group", "Condition", "Units", "Concentration")], ts)
ts <- ts %>% dplyr::group_by(Group, Condition, Units, Concentration) %>%
  dplyr::summarise(dplyr::across(.fns = median))
ts$Concentration <- sub("/.0+$", "", ts$Concentration) %>% as.numeric()
ts$Name <- ts$Condition
ts <- as.data.frame(ts)
ts$NameGroup <- paste(ts$Name, ts$Group)
ts[ts$NameGroup %in% names(table(ts$NameGroup))[table(ts$NameGroup) != 1], "Name"] <- with(ts[ts$NameGroup %in% names(table(ts$NameGroup))[table(ts$NameGroup) != 1], ], paste(Condition, sub("/.0+$", "", Concentration), sub("vol%", "%", Units)))
ts <- ts[, !colnames(ts) %in% c("Group", "Condition", "Units", "Concentration", "Name", "NameGroup")] %>% as.matrix()

fs[abs(ts) < 0.5] <- 0
fs[abs(fs) < 0.5] <- 0

fsn.conds <- rowSums(abs(fs[, iols_genenames])) != 0
fsn <- fs[fsn.conds, iols_genenames]
dim(fsn)
fsn <- t(fsn)
fs.net <- network(fsn,
                  bipartite = TRUE,
                  ignore.eval = FALSE,
                  names.eval = "weights"
                  )
fs.net %v% "vertex.names" <- c(rownames(fsn), fs.meta$Name[fsn.conds])
fs.net %v% "vertex.names" %<>% sub("Ethylmethylimidazolium", x = ., replacement = "EMI") %>% sub("Chloramphenicol", x = ., replacement = "Chl") %>% sub("Phosphomycin", x = ., replacement = "Pho") %>% sub("Spectinomycin", x = ., replacement = "Spe") %>% sub("Nalidixate", x = ., replacement = "Nal") %>% sub("Neomycin", x = ., replacement = "Neo") %>% sub("Furfuraldehyde", x = ., replacement = "Fur") %>% sub("Bacitracin", x = ., replacement = "Bac") %>% sub("Gentamicin", x = ., replacement = "Gen") %>% sub("Cisplatin", x = ., replacement = "Cis") %>% sub("Doxycycline", x = ., replacement = "Dox") %>% sub(x = ., "Methylglyoxal", "MeGly") %>% sub(x = ., "Thallium acetate", "Tl acetate") %>% sub(x = ., "Cobalt chloride", "CoCl2") %>% sub(x = ., "Aluminium chloride", "AlCl3")
col = c("actor" = pal[1], "event" = "grey")
node.size <- abs(fs.net %e% "weights") %>% minmax(., 1, 9)
set.edge.attribute(fs.net, "color", ifelse(fs.net %e% "weights" > 0, pal[2], pal[3]))
set.vertex.attribute(fs.net, "node", ifelse(fs.net %v% "vertex.names" %in% iols_genenames, 60, 17))
set.vertex.attribute(fs.net, "text", ifelse(fs.net %v% "vertex.names" %in% iols_genenames, 16, 10))
set.vertex.attribute(fs.net, "shape", c(rep(16, length(iols_genenames)), rep(c(15, 19, 17, 18), times = as.factor(fs.meta[fsn.conds, "Group"]) %>% as.numeric() %>% table())))
#dev.new(width=100000, height=100000)
#set.seed(2)
ggn <- ggnet2(fs.net, seed = 1,
       color = "mode", layout.par = list("niter" = 1500,
                                         "area" = 30000000,
                                         "repulse.rad" = 30000000),
       palette = col,
       shape = "shape",
       size = "node",
       label.size = "text",
       edge.size = .9*(abs(fs.net %e% "weights")),
       label = FALSE,
       legend.position = "none",
       edge.color = "color",
       edge.alpha = scaleWeight(fs.net %e% "weights", alpha = .8, nmin = .35),
       layout.exp = 0.0001) +
  geom_point(aes(color = color, shape = as.character(shape)), size = 18, color = "white") +
  geom_point(aes(color = color, shape = as.character(shape)), size = 14, alpha = 0.5, color = rep(c(pal[1], pal[-c(1, 2, 3)]), times = c(length(iols_genenames), table(fs.meta[fsn.conds, "Group"])))) +
  geom_point(aes(color = color, shape = as.character(shape)), size = 16, color = rep(c(pal[1], pal[-c(1, 2, 3)]), times = c(length(iols_genenames), table(fs.meta[fsn.conds, "Group"])))) +
  geom_text(aes(label = label), color = c(rep("white", 12), rep("black", 66)), fontface = "bold", size = 4.5) 

ragg::agg_jpeg("price_figure_7e.jpg", width = 4700, height = 3200, res = 300)
ggn
invisible(dev.off())

