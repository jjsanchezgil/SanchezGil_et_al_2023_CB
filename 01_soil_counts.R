require(magrittr, quietly = TRUE)
require(dplyr, quietly = TRUE)
require(multcomp, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(ggtext, quietly = TRUE)
require(ggpubr, quietly = TRUE)
require(patchwork, quietly = TRUE)
require(extrafont, quietly = TRUE)
require(ggsci, quietly = TRUE)
require(scales, quietly = TRUE)
require(emmeans, quietly = TRUE)

extrafont::font_import("/Windows/Fonts/", pattern = "Montserrat.*", prompt = F)
extrafont::loadfonts(device = "win")

#2.46 * 10^5 copies/ng Salinibacter DNA
#4 ng were added
salinibacter.conversion <- 2.46e5 *4

##Define the sequences in the study
sq_SAL <- "TGGGGAATCTTGCACAATGGGGTCACCCCTGATGCAGCCATGCCGCGTGGAGGAAGACACCCCTATGGGGCGTAAACTCCTTTTCTGAATGAAGAAACCCCTGTAGCTTCAGGGCGCGACGGTAGTTCAGGAATAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGTGCAGGCGGGGCAGCAAGTCGGATGTGAAACCCCATGGCTTAACCATGGAGGTGCATTCGAAACTGTTGCTCTTGAGTCCCGGAGAGGCTGTCGGAATTCGTGGTGTAGCGGTGAAATGCGTAGATATCACGAGGAACACCAGAGGCGAAAGCGGACAGCTGGACGGGTACTGACGCTCAGGCACGAAAGCGTGGGGAGCAAACA"
sq_CLP <- "TGGGGAATTTTCCGCAATGGGCGAAAGCCTGACGGAGCAATGCCGCGTGGAGGTAGAAGGCCTACGGGTCCTGAACTTCTTTTCCCAGAGAAGAAGCAATGACGGTATCTGGGGAATAAGCATCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACCAAGCTTGAGTACGGTAGGGGCAGAGGGAATTTCCGGTGGAGCGGTGAAATGCGTAGAGATCGGAAAGAACACCAACGGCGAAAGCACTCTGCTGGGCCGACACTGACACTGAGAGACGAAAGCTAGGGGAGCGAATG"
sq_MIT <- "TGGGGAATCTTGGACAATGGGCGAAAGCCCGATCCAGCAATATCGCGTGAGTGAAGAAAGGCAATGCCGCTTGTAAAGCTCTTTCGTCGAGTGCGCGATCATGACAGGACTCGAGGAAGAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAAAACGGGGGGGGCAAGTGTTCTTCGGAATGACTGGGCGTAAAGGGCACGTAGGCGGTGAATCGGGTTGAAAGTGAAAGTCGCCAAAAAGTGGCGGAATGCTTTCGAAACCAATTCACTTGAGTGAGACAGAGGAGAGTGGAATTTCGTGTGGAGGGGTGAAATCTACAGATCTACGAAGGAACGCCAAAAGCGAAGGCAGCTCTCTGGGTCCCTACCGACGCTGGGGTGCGAAAGCATGGGGAGCGAACG"
sq_CHA <- "TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCACTTTAAGTTGGGAGGAAGGGCAGTTACCTAATACGTGATTGTTTTGACGTTACCGACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTGGCAAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACA"
sq_417 <- "TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCACTTTAAGTTGGGAGGAAGGGCAGTTACCTAATACGTGATTGTTTTGACGTTACCGACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACTGACTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTAATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACA"
sq_365 <- "TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCACTTTAAGTTGGGAGGAAGGGCAGTAAGCGAATACCTTGCTGTTTTGACGTTACCGACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGAATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTGGCAAGCTAGAGTACGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACA"
sq_358 <- "TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCACTTTAAGTTGGGAGGAAGGGCAGTAAGCGAATACCTTGCTGTTTTGACGTTACCGACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGAATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTGGCAAGCTAGAGTACGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACA"
sq_317 <- "TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCACTTTAAGTTGGGAGGAAGGGTTGTAGATTAATACTCTGCAATTTTGACGTTACCGACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACA"
sq_158 <- "TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCACTTTAAGTTGGGAGGAAGGGCAGTAAGTTAATACCTTGCTGTTTTGACGTTACCGACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTCGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTGGCGAGCTAGAGTACGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACA"
sq_134 <- "TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCACTTTAAGTTGGGAGGAAGGGTTGTAGATTAATACTCTGCAATTTTGACGTTACCGACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTCGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGTCGAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACA"
sq_134.2 <- "TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCACTTTAAGTTGGGAGGAAGGGTTGTAGATTAATACTCTGCAATTTTGACGTTACCGACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTCGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGTCGAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACA"
sq_V3V4 <- c("SAL" = sq_SAL, "CLP" = sq_CLP, "MIT" = sq_MIT,
             "358" = sq_358, "417" = sq_417, "CHA" = sq_CHA, "134.2" = sq_134.2,
             "158" = sq_158, "365" = sq_365, "134" = sq_134, "317" = sq_317)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

#working directory where sample data is
setwd("./")

weights <- read.csv("root_weights.csv", sep = ";", dec = ",")
st.full <- readRDS("sequence.table.clean.RDS") %>% as.data.frame()
#REMOVE 365!
st.full <- st.full[(rownames(st.full) %>% strsplit("-") %>% sapply('[', 1) %>% sub(x = ., "CTR_", "")) != "365", ]
st <- st.full[, unique(sq_V3V4)]
dim(st) #Should be 144 x 9

st$Sample <- rownames(st) %>% sub(x = ., "-R1.fastq.gz", "")
st$Treatments <- rownames(st) %>% strsplit("-") %>% sapply('[', 1) %>% sub(x = ., "CTR_", "")
st$Compartment <- rownames(st) %>% strsplit("-") %>% sapply('[', 3) %>% sub(x = ., "CTR_", "")
st$Inoculation <- ifelse(startsWith(rownames(st), "CTR"), "Uninoculated", "Inoculated")
st$sCounts <- st[[sq_V3V4[["SAL"]]]]
st$Mito <- st[, sq_V3V4[["MIT"]]]
st$Chlp <- st[, sq_V3V4[["CLP"]]]
st$Depth <- rowSums(st.full)
st$mDepth <- with(st, Depth - Mito - Chlp)

st$sConversion <- salinibacter.conversion / st$sCounts
st$pCounts <- with(st, Mito + Chlp)
st$pCopies <- with(st, pCounts * sConversion)
st$Counts <- sapply(rownames(st), \(x) st[x, sq_V3V4[[st[x, "Treatments"]]]])
st$Rabn <- with(st, 100 * Counts / mDepth)
st$Copies <- with(st, Counts * sConversion)
st$RCopies <- with(st, Copies / pCopies)


weights.plane <- weights
weights.vec <- weights.plane[["Weight"]]
names(weights.vec) <- weights.plane[["SampleID"]]
st$Weight <- weights.vec[st$Sample %>% sub(x = ., "CTR_...-", "CTR-")]
st$Density <- with(st, Copies/Weight)
st$SalRatio <- with(st, sCounts/Depth)
st$mLoad <- with(st, mDepth * sConversion / Weight)

st <- st[, 10:ncol(st)]
st$Treatments <- as.character(st$Treatments)
st$Compartment %<>% sub(x = ., "P", "Root") %>%
  sub(x = ., "S", "Rhizosphere") %>%
  sub(x = ., "B", "Soil") %>%
  factor(., levels = c("Soil", "Rhizosphere", "Root"))

mocks <- st %>% filter(Inoculation == "Uninoculated") %>%
  dplyr::group_by(Treatments) %>%
  summarise(meanMock = mean(Density)) %>%
  as.data.frame()
rownames(mocks) <- mocks$Treatments
st$cDensity <- st$Density - mocks[st$Treatments, "meanMock"]
st$cDensity1 <- st$cDensity + 1

st.density <- st[st$Inoculation == "Inoculated", ]
st.density <- st.density[st.density$Density > 0 &
                           st.density$sCounts > 10, ]

#Remove based on the density of Salinibacter
#to avoid errors in the inoculation density
#Remove them based on the 99th and 1st quantile
lower.bound <- with(st.density, sCounts/Depth) > quantile(with(st.density, sCounts/Depth), .01, na.rm = T)
upper.bound <- with(st.density, sCounts/Depth) < quantile(with(st.density, sCounts/Depth), .99, na.rm = T)
st.density <- st.density[lower.bound & upper.bound, ]

outliers.hi <- st.density %>% dplyr::group_by(Compartment, Treatments) %>% dplyr::mutate(higher = Density/median(Density) < 10) %>% '[['("higher")
outliers.lo <- st.density %>% dplyr::group_by(Compartment, Treatments) %>% dplyr::mutate(higher = median(Density)/Density < 10) %>% '[['("higher")
st.density <- st.density[outliers.hi & outliers.lo, ]

xlabs <- c("WCS134", "RS158", "WCS317", "WCS358", "WCS417", "CHA0")
names(xlabs) <- c("134", "158", "317", "358", "417", "CHA")

st.density$TreatmentsFull <- xlabs[st.density$Treatments]
st.density$TreatmentsFull <- factor(st.density$TreatmentsFull, levels = xlabs)
st.density$Treatments %<>% factor
####################
#Across compartments
st.treatments <- c()
st.compartments <- c()
st.labels <- c()

for (cmp in unique(st.density$Compartment)) {
  m <- MASS::glm.nb(cDensity1 ~ Treatments, data = st.density[st.density$Compartment == cmp,]) %>%
    emmeans(object = ., specs = "Treatments") %>%
    cld(., adjust = "sidak", alpha = .05) %>%
    as.data.frame %>%
    '['(c("Treatments", ".group"))
  mx.nm <- m$.group %>% as.numeric %>% as.character %>%
    strsplit(., "") %>% sapply(., max) %>% max
  ngroups <- m$.group %>% as.numeric %>% as.character %>%
    strsplit(., "") %>% sapply(\(x) as.numeric(mx.nm) - as.numeric(x) + 1) %>%
    sapply(\(x) letters[sort(x)]) %>% sapply(paste0, collapse = "")
  m$ngroups <- ngroups

  st.treatments <- c(st.treatments, m$Treatments)
  st.compartments <- c(st.compartments, rep(cmp, nrow(m)))
  st.labels <- c(st.labels, m$ngroups)
}
st.labels.full <- data.frame(
  Treatments = st.treatments,
  Compartment = st.compartments,
  Label = st.labels
)
st.labels.full$Compartment %<>% factor(., levels = c("Soil", "Rhizosphere", "Root"))
st.labels.full$TreatmentsFull <- xlabs[st.labels.full$Treatments]


hex_codes.isolates <- c(
  "grey40",
          "grey40",
          "grey40",
          "grey40",
          pal_material("amber")(4)[4],
          pal_material("cyan")(9)[4])


##Figure 1a ====
#st.density$Compartment %<>% factor(., levels = c("Soil", "Rhizosphere", "Root"))
cc4 <- ggplot(st.density,
              aes(x = TreatmentsFull, y = cDensity, fill = TreatmentsFull)) +
  # geom_boxplot(outlier.shape = NA, alpha = .5, lwd = .2, width = .6) +
  geom_jitter(aes(color = TreatmentsFull), size = 1.5, width = .2) +
  facet_wrap(~Compartment) +
  geom_errorbar(data = st.density %>%
                  dplyr::group_by(TreatmentsFull, Compartment) %>%
                  dplyr::transmute(cDensity = median(cDensity),
                                   std = sd(cDensity)),
                aes(ymin = cDensity - std,
                    ymax = cDensity + std,
                    x = TreatmentsFull,
                    color = TreatmentsFull), linewidth = .7) +
  geom_text(mapping = aes(x = TreatmentsFull, y = 2e8, label = Label),
            data = st.labels.full,
            family = "Montserrat", size = 5) +
  ylab("*16S rRNA* gene copies per gram sample") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.y = element_markdown(),
        strip.text.x = element_markdown(color = "black",
                                        face = "bold"),
        legend.position = "none",
        text = element_text(family = "Montserrat", size = 18),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray80")) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l",
                      long = unit(4, "mm"),
                      mid = unit(3, "mm"),
                      short = unit(1, "mm"),
                      size = .2) +
  scale_color_manual(values = hex_codes.isolates)

ragg::agg_jpeg("7isolates_figure_1a.jpg", width = 5000, height = 2000, res = 300)
cc4
invisible(dev.off())

st.density$Behaviour <- "Increasing or equal"
st.density[st.density$Treatments %in% c("358", "158", "317"), "Behaviour"] <- "Decreasing"
cc5 <- ggplot(st.density,
       aes(x = Compartment, y = cDensity, fill = TreatmentsFull)) +
  # geom_boxplot(outlier.shape = NA, alpha = .5, lwd = .2, width = .6) +
  geom_jitter(aes(color = TreatmentsFull), size = 1.5, width = .1) +
  geom_errorbar(data = st.density %>%
                  dplyr::group_by(TreatmentsFull, Compartment, Behaviour) %>%
                  dplyr::transmute(cDensity = median(cDensity),
                                   std = sd(cDensity)),
                aes(ymin = cDensity - std,
                    ymax = cDensity + std,
                    x = Compartment,
                    color = TreatmentsFull), linewidth = .7, width = .25) +
  facet_wrap(~Behaviour) +
  geom_text(data = st.density %>%
              dplyr::filter(Compartment == "Root") %>%
              dplyr::group_by(TreatmentsFull, Compartment, Behaviour) %>%
              dplyr::summarise(cDensity = median(cDensity)),
            mapping = aes(x = 3.2, y = cDensity, label = TreatmentsFull),
            family = "Montserrat", size = 3, hjust = 0) +
  geom_line(data = st.density %>%
              dplyr::group_by(TreatmentsFull, Compartment, Behaviour) %>%
              dplyr::transmute(cDensity = median(cDensity),
                               std = sd(cDensity)),
            aes(y = cDensity,
                x = Compartment,
                group = TreatmentsFull,
                color = TreatmentsFull), linewidth = .7) +
  ylab("*16S rRNA* gene copies per gram sample") +
  expand_limits(x = c(1, 4)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.y = element_markdown(),
        strip.text.x = element_markdown(color = "black",
                                        face = "bold"),
        legend.position = "none",
        text = element_text(family = "Montserrat", size = 18),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray80")) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l",
                      long = unit(4, "mm"),
                      mid = unit(3, "mm"),
                      short = unit(1, "mm"),
                      size = .2) +
  scale_color_manual(values = hex_codes.isolates)
ragg::agg_jpeg("7isolates_trajectories.jpg", width = 2500, height = 2000, res = 300)
cc5
invisible(dev.off())


#Detection limit ----
s8 <- ggplot(st.density, aes(x = Compartment, y = sConversion/Weight)) +
  geom_boxplot() +
  geom_jitter() +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  ylab("Detection limit (copies / g)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.y = element_markdown(),
        strip.text.x = element_markdown(color = "black",
                                        face = "bold"),
        legend.position = "none",
        text = element_text(family = "Montserrat", size = 25),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray80"))
ragg::agg_jpeg("figure_S8.jpg", width = 2500, height = 2000, res = 300)
s8
invisible(dev.off())


