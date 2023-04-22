require(magrittr)
require(reshape2)
require(stringr)
require(dplyr)
require(ggplot2)
require(ggtext)
require(extrafont)

extrafont::font_import("/Windows/Fonts/", pattern = "Montserrat.*", prompt = F)
extrafont::loadfonts(device = "win")

scale.vec <- function(vec, use = "mean") {
  if (use == "mean") {
    scaled <- (vec - mean(vec, na.rm = TRUE)) / sd(vec, na.rm = TRUE)
  } else if (use == "median") {
    scaled <- (vec - median(vec, na.rm = TRUE)) / sd(vec, na.rm = TRUE)
  } else stop("Parameter 'use' should be 'mean' or 'median'")
  return(scaled)
}

geom.mean <- function(vec) {
  return(exp(mean(log(vec + 1))))
}

tpm_raw <- readRDS("vesga_2021_tpm.RDS")

iolgenes <- 2627:2638 %>% as.character() %>% paste0("PPRCHA0_", .)
iolnames <- c("iolR", "iolC", "iolE", "iolB", "iolI", "iolD", "iolG", "hypothetical", "iatA", "iatC", "iatB", "dksA")
tpm_iol <- tpm_raw[iolgenes,]
tpm_iol <- reshape2::melt(tpm_iol)
tpm_iol$Var1 %<>% iolnames[.] %>% as.factor() %>% ordered(levels = iolnames)
names(tpm_iol) <- c("Gene", "Medium", "TPM")
tpm_iol$Medium %<>% as.character() %>%
  sapply(function(x) ifelse(startsWith(x, "PL"),
                            yes = strsplit(x, "_")[[1]],
                            no = stringr::str_remove(x, pattern = "[1-9]")))
qsgenes <- c("PPRCHA0_0923", "PPRCHA0_1214", "PPRCHA0_2929", "PPRCHA0_4134",
             "PPRCHA0_4135", "PPRCHA0_4136", "PPRCHA0_4137", "PPRCHA0_4216", 
             "PPRCHA0_4229", "PPRCHA0_4230", "PPRCHA0_4231", "PPRCHA0_2823", 
             "PPRCHA0_2713", "PPRCHA0_1053", "PPRCHA0_1266", "PPRCHA0_5330", 
             "PPRCHA0_5704", "PPRCHA0_2619", "PPRCHA0_2620", "PPRCHA0_2621", 
             "PPRCHA0_3216", "PPRCHA0_5218", "PPRCHA0_5219", "PPRCHA0_5220", 
             "PPRCHA0_1470", "PPRCHA0_1473", "PPRCHA0_1634", "PPRCHA0_1635",
             "PPRCHA0_1689", "PPRCHA0_1651", "PPRCHA0_5883", "PPRCHA0_1690",
             "PPRCHA0_1691", "PPRCHA0_1692", "PPRCHA0_1254", "PPRCHA0_4363",
             "PPRCHA0_3621", "PPRCHA0_3622", "PPRCHA0_3623", "PPRCHA0_3624",
             "PPRCHA0_5858", "PPRCHA0_5860", "PPRCHA0_5859", "PPRCHA0_5861",
             "PPRCHA0_5862", "PPRCHA0_5857", "PPRCHA0_5856", "PPRCHA0_5855",
             "PPRCHA0_4811", "PPRCHA0_2326", "PPRCHA0_3580", "PPRCHA0_4497",
             "PPRCHA0_2825", "PPRCHA0_2826", "PPRCHA0_2827", "PPRCHA0_2828",
             "PPRCHA0_2829", "PPRCHA0_2830", "PPRCHA0_2831", "PPRCHA0_2832",
             "PPRCHA0_3736", "PPRCHA0_3505", "PPRCHA0_3506", "PPRCHA0_3507", 
             "PPRCHA0_3513", "PPRCHA0_3510", "PPRCHA0_3509", "PPRCHA0_3514",
             "PPRCHA0_3515", "PPRCHA0_3257", "PPRCHA0_3078", "PPRCHA0_3079",
             "PPRCHA0_3080", "PPRCHA0_3081", "PPRCHA0_3082", "PPRCHA0_3083",
             "PPRCHA0_3084", "PPRCHA0_6002", "PPRCHA0_3016", "PPRCHA0_6004",
             "PPRCHA0_3012", "PPRCHA0_5994", "PPRCHA0_5995", "PPRCHA0_5996",
             "PPRCHA0_5998", "PPRCHA0_5999", "PPRCHA0_5960", "PPRCHA0_5989",
             "PPRCHA0_5990", "PPRCHA0_5991", "PPRCHA0_5992",
             iolgenes)
qsnames <- c("rpoN", "rpoS", "pvdQ", "fpvA", 
             "pvdD", "pvdJ", "pvdI", "pvdH", 
             "pvdL", "pvdS", "pvdY", "pltR", 
             "phzF_1", "rhlB", "rhlG", "rhlE_1", 
             "rhlE_2", "hcnA", "hcnB", "hcnC", 
             "aprA", "pilA", "pilC", "pilD", 
             "algU", "mucD", "flgF", "flgG",
             "fliA", "fliC", "fliL", "cheY",
             "cheZ", "cheA", "recA", "gyrA", 
             "prnA", "prnB", "prnC", "prnD",
             "phlA", "phlB", "phlC", "phlD",
             "phlE", "phlF", "phlG", "phlH",
             "phzF_2", "phzF_3", "gacA", "gacS",
             "pltA", "pltB", "pltC", "pltD",
             "pltE", "pltF", "pltG", "pltZ",
             "phy", "pchA", "pchB", "pchC",
             "pchD", "pchE", "pchF", "pchR",
             "fetA", "mtlR", "mltE", "mltF",
             "mltG", "mltK", "mltD", "mltY",
             "mltZ", "vgrG1a", "vgrG1b", "rhsB/A",
             "ghh1", "tssA", "tssB", "tssC",
             "tssE", "tssF", "tssG", "tssM",
             "tssL", "tssK", "tssJ",
             iolnames)
qsgroups <- c("General regulators", "General regulators", "Pyoverdin", "Pyoverdin",
              "Pyoverdin", "Pyoverdin", "Pyoverdin", "Pyoverdin",
              "Pyoverdin", "Pyoverdin", "Pyoverdin", "Pyoluteorin",
              "Phenazine", "Rhamnolipids", "Rhamnolipids", "Rhamnolipids",
              "Rhamnolipids", "HCN", "HCN", "HCN",
              "Protease", "Pilus", "Pilus", "Pilus",
              "Biofilm", "Biofilm", "Flagella", "Flagella",
              "Flagella", "Flagella", "Flagella", "Chemotaxis",
              "Chemotaxis", "Chemotaxis", "Housekeeping", "Housekeeping",
              "Pyrrolnitrin", "Pyrrolnitrin", "Pyrrolnitrin", "Pyrrolnitrin",
              "DAPG", "DAPG", "DAPG", "DAPG",
              "DAPG", "DAPG", "DAPG", "DAPG",
              "Phenazine", "Phenazine", "General regulators", "General regulators",
              "Pyoluteorin", "Pyoluteorin", "Pyoluteorin", "Pyoluteorin",
              "Pyoluteorin", "Pyoluteorin", "Pyoluteorin", "Pyoluteorin",
              "Inositol metabolism", "Pyochelin", "Pyochelin", "Pyochelin",
              "Pyochelin", "Pyochelin", "Pyochelin", "Pyochelin",
              "Pyochelin", "Mannitol metabolism", "Mannitol metabolism", "Mannitol metabolism",
              "Mannitol metabolism", "Mannitol metabolism", "Mannitol metabolism", "Mannitol metabolism",
              "Mannitol metabolism", "T6SS", "T6SS", "T6SS",
              "T6SS", "T6SS", "T6SS", "T6SS",
              "T6SS", "T6SS", "T6SS", "T6SS",
              "T6SS", "T6SS", "T6SS",
              rep("Inositol metabolism", length(iolnames)))
names(qsnames) <- qsgenes
names(qsgroups) <- qsgenes
tpm <- tpm_raw
tpm_qs <- tpm[, colnames(tpm) %in% qsgenes] %>% t()
tpm_qs <- tpm_raw[qsgenes,]
tpm_qs <- reshape2::melt(tpm_qs)
tpm_qs$Var1 %<>% qsnames[.] %>% as.factor() %>% ordered(levels = qsnames)
tpm_qs$Groups <- qsgroups[tpm_qs$Var1] 
names(tpm_qs) <- c("Gene", "Medium", "TPM", "Groups")
tpm_qs$Medium %<>% as.character() %>% 
  sapply(function(x) ifelse(startsWith(x, "PL"), 
                            yes = strsplit(x, "_")[[1]],
                            no = stringr::str_remove(x, pattern = "[1-9]")))
tpm_qs <- tpm_qs[!tpm_qs$Medium %in% c("GM", "LB"),]
tpm_qs$Medium %<>% 
  sub("LB", "LB<br>medium", .) %>%
  sub("GM", "GM<br>medium", .) %>%
  sub("PL24", "*P. xylostella* gut<br>24h", .) %>%
  sub("PL36", "*P. xylostella* gut<br>36h", .) %>%
  sub("Galleria", "*G. mellonella*<br>hemolymph", .) %>%
  sub("Wheat", "Wheat<br>rhizosphere", .)

tpm_qs_scaled <- tpm_qs[tpm_qs$Gene != "phy",] %>% group_by(Medium, Gene) %>% mutate(TPMs = mean(TPM))
tpm_qs_scaled %<>% group_by(Gene) %>% mutate(mTPMs = scale.vec(log2(TPMs + 1)))
tpm_qs_scaled2 <- tpm_qs[tpm_qs$Gene != "phy",] %>% group_by(Medium, Gene) %>% mutate(TPMs = geom.mean(TPM)) 
tpm_qs_scaled2 %<>% group_by(Gene) %>% mutate(mTPMs = scale(log2(TPMs + 1)))
tpm_qs_scaled2$Groups %<>% sub(x = ., "Inositol metabolism", "Inositol<br>metabolism")
gtpm <- ggplot(tpm_qs_scaled2[!tpm_qs_scaled2$Groups %in% c("Protease",
                                                            "Mannitol metabolism",
                                                            "General regulators",
                                                            "Phenazine",
                                                            "Pyrrolnitrin",
                                                            "Housekeeping",
                                                            "T6SS",
                                                            "Rhamnolipids",
                                                            "Pyoluteorin", 
                                                            "Pyoverdin"), ],
               aes(x = Medium, fill = mTPMs, y = Gene, group = Groups)) + 
  geom_tile() + 
  facet_grid(rows = vars(Groups), scales = "free", space = "free") + 
  labs(fill = "Gene-scaled log<sub>2</sub>(TPM)") +
  theme_minimal() +
  scale_x_discrete(expand = c(0, 0)) +
  theme(text = element_text(size = 18, family = "Montserrat"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        strip.text.y = element_markdown(angle = 0, hjust = 0),
        axis.text.x = element_markdown(),
        axis.text.y = element_markdown(face = "italic"),
        legend.title = element_markdown(),
        legend.position = "top",
        panel.grid.major = element_line(colour = "gray95"),
        panel.background = element_rect(fill = "white")) + 
  scale_fill_gradient2(low = "white", high = "#FA2E00", midpoint = 0)
ragg::agg_jpeg("vesga_figure_7f.jpg", width = 3500, height = 4000, res = 300)
gtpm
invisible(dev.off())


###REDUCED----
gtpm2 <- ggplot(tpm_qs_scaled2[tpm_qs_scaled2$Groups == c("Inositol<br>metabolism"), ],
               aes(x = Medium, fill = mTPMs, y = Gene, group = Groups)) + 
  geom_tile() + 
  facet_grid(rows = vars(Groups), scales = "free", space = "free") + 
  labs(fill = "Gene-scaled log<sub>2</sub>(TPM)") +
  theme_minimal() +
  scale_x_discrete(expand = c(0, 0)) +
  theme(text = element_text(size = 12, family = "Montserrat"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        strip.text.y = element_blank(),
        axis.text.x = element_markdown(size = 8),
        axis.text.y = element_markdown(face = "italic"),
        legend.title = element_markdown(),
        legend.position = "top",
        panel.grid.major = element_line(colour = "gray95"),
        panel.background = element_rect(fill = "white")) + 
  scale_fill_gradient2(low = "white", high = "#FA2E00", midpoint = 0)
ragg::agg_jpeg("vesga_figure_7f_simplified.jpg", width = 1500, height = 1800, res = 300)
gtpm2
invisible(dev.off())
