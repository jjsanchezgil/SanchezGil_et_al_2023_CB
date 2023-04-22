require(magrittr)
require(ggplot2)
require(ggpubr)
require(ggsci) 
require(ggtext) 
library(patchwork)
library(scales)
library(ggsignif)
library(reshape2)
library(rstatix)
library(MASS)
library(emmeans)
library(multcomp)
extrafont::font_import("/Windows/Fonts/", pattern = "Montserrat.*", prompt = F)
extrafont::loadfonts(device = "win")

nice.log.breaker = function(x) {
  brks <- extended_breaks(Q = c(1, 5))(log10(x))
  10^(brks[brks %% 1 == 0])
}
nice.breaker = function(x) {
  brks <- extended_breaks(Q = c(1, 5))(x)
  brks[brks %% 1 == 0]
}

i.cfu <- readxl::read_excel("First exp_sample_sheet_cfu_individual.xlsx") %>% as.data.frame()
colnames(i.cfu)
m.cfu <- readxl::read_excel("First_exp_sample_sheet_cfu_mixed.xlsx") %>% as.data.frame()
colnames(m.cfu)

i.cfu$Genotype %<>% factor(., levels = c("Control", "WT", "iol"))
i.cfu$Label <- paste(i.cfu$Source, i.cfu$Time) %>% sub(x = ., "Bulk", "Soil")
i.cfu$Label %<>% factor(., levels = c("Soil 0 dpi", "Soil 21 dpi", "Root 21 dpi"))
m.cfu$Label %<>% factor(., levels = c("Bulk 0 dpi", "Bulk 21 dpi", "Root 21 dpi"))

#Labels for genotypes 
gen.labels <- c("&Delta;*iol*", "WT")
names(gen.labels) <- c("iol", "WT")

####################################
##Individual survival at time point 21 ----
median.t0 <- i.cfu[i.cfu$Time == "0 dpi", ] %>% 
  dplyr::group_by(Genotype) %>% 
  dplyr::summarise(med.cfu = median(CFU_g)) %>% 
  as.data.frame()
i.cfu[i.cfu$Time == "21 dpi" & 
        i.cfu$Genotype == "iol", 
      "RelSurvival"] <- i.cfu[i.cfu$Time == "21 dpi" & i.cfu$Genotype == "iol", 
                              "CFU_g"] / median.t0[median.t0$Genotype == "iol", "med.cfu"]
i.cfu[i.cfu$Time == "21 dpi" & 
        i.cfu$Genotype == "WT", 
      "RelSurvival"] <- i.cfu[i.cfu$Time == "21 dpi" & i.cfu$Genotype == "WT", 
                              "CFU_g"] / median.t0[median.t0$Genotype == "WT", "med.cfu"]

ggplot(i.cfu[i.cfu$Time == "21 dpi",], 
       aes(x = Genotype, y = RelSurvival, fill = Genotype)) +
  ylab("CFU per gram") +
  facet_wrap(~Source) +
  geom_jitter(aes(color = Genotype)) +
  geom_errorbar(data = i.cfu %>%
                  dplyr::filter(Time == "21 dpi") %>%
                  dplyr::group_by(Genotype, Source) %>%
                  dplyr::mutate(dens = median(RelSurvival),
                                std = sd(dens)),
                aes(ymin = dens - std,
                    ymax = dens + std,
                    x = Genotype,
                    color = Genotype), linewidth = 1) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_markdown(angle = 30, hjust = 1),
        axis.title.y = element_text(),
        legend.position = "none",
        text = element_text(family = "Montserrat", size = 27),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray80")) +
  scale_x_discrete(labels = gen.labels) 

####################################
# Ratio in competition at time point 21 ----
m.cfu$Ratio <- m.cfu$CFU_g_WT / m.cfu$CFU_g_KO
ggplot(m.cfu, 
       aes(x = Label, y = log2(Ratio), fill = Label)) +
  ylab("log<sub>2</sub> WT:??*iol* ratio") + 
  geom_jitter(aes(color = Label)) +
  geom_errorbar(data = m.cfu %>%
                  dplyr::group_by(Label) %>%
                  dplyr::mutate(dens = median(log2(Ratio)),
                                std = sd(log2(dens)),
                                meandens = mean(log2(Ratio))),
                aes(ymin = dens - std,
                    ymax = dens + std,
                    x = Label,
                    color = Label), linewidth = 1) +
  stat_summary(aes(x = Label, y = log2(Ratio)),
               fun = mean, geom = "point", shape = "-", size = 15, inherit.aes = F) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_markdown(angle = 30, hjust = 1),
        axis.title.y = element_markdown(),
        legend.position = "none",
        text = element_text(family = "Montserrat", size = 27),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray80")) +
  scale_x_discrete(labels = gen.labels) +
  # scale_y_continuous(breaks = seq(0, 10, 1), limits = c(0, 10)) +
  stat_compare_means(ref.group = "Bulk 0 dpi") +
  stat_compare_means(ref.group = "Bulk 21 dpi", vjust = 4)



####################################
# Ratio in competition for EXPERIMENT 1 ----
m.cfu1 <- read.csv("Second_exp_sample_sheet_cfu_mixed.csv", sep = ";") 
m.cfu1$Ratio <- m.cfu1$CountsWT / m.cfu1$CountsKO
m.cfu1$Time <- paste(m.cfu1$TP, "dpi") %>% sub(x = ., "7", "21")
m.cfu1$Label <- paste(m.cfu1$Comp, m.cfu1$Time)
m.cfu1 <- m.cfu1[m.cfu1$Label != "Root 0 dpi", ]
m.cfu1$Label %<>% factor(., levels = c("Soil 0 dpi", "Soil 21 dpi", "Root 21 dpi"))

ggplot(m.cfu1, 
       aes(x = Label, y = log2(Ratio), fill = Label)) +
  ylab("log<sub>2</sub> WT:??*iol* ratio") + 
  geom_jitter(aes(color = Label)) +
  geom_errorbar(data = m.cfu1 %>%
                  dplyr::group_by(Label) %>%
                  dplyr::mutate(dens = median(log2(Ratio), na.rm = T),
                                std = sd(dens)),
                aes(ymin = dens - std,
                    ymax = dens + std,
                    x = Label,
                    color = Label), linewidth = 1) +
  stat_summary(aes(x = Label, y = log2(Ratio)),
               fun = mean, geom = "point", shape = "-", size = 15, inherit.aes = F) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_markdown(angle = 30, hjust = 1),
        axis.title.y = element_markdown(),
        legend.position = "none",
        text = element_text(family = "Montserrat", size = 27),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray80")) +
  scale_x_discrete(labels = gen.labels) +
  stat_compare_means(ref.group = "Soil 0 dpi") +
  stat_compare_means(ref.group = "Soil 21 dpi", vjust = 2)



# MIX FULL dataset ----
m.cfu.all <- rbind(data.frame(m.cfu1[, c("Label", "Ratio")], 
                              Experiment = 1, 
                              WT = m.cfu1[, "CountsWT"],
                              KO = m.cfu1[, "CountsKO"],
                              Prop = 100*with(m.cfu1, CountsWT/(CountsWT + CountsKO))), 
                   data.frame(m.cfu[, c("Label", "Ratio")], 
                              Experiment = 2, 
                              WT = m.cfu[, "CFU_g_WT"],
                              KO = m.cfu[, "CFU_g_KO"],
                              Prop = 100*with(m.cfu, CFU_g_WT/(CFU_g_WT + CFU_g_KO))))
m.cfu.all$Label %<>% sub(x = ., "Bulk", "Soil")
m.cfu.all <- m.cfu.all[!is.na(m.cfu.all$Ratio) &
                         !is.infinite(m.cfu.all$Ratio), ]
m.cfu.all$logRatio <- log2(m.cfu.all$Ratio)
m.cfu.all$Label %<>% factor(., levels = c("Soil 0 dpi", "Soil 21 dpi", "Root 21 dpi"))
m.cfu.all$Experiment %<>% factor(.)
#summary vs Experiment
#Check if Experiment affects the effect between Label's (doesn't affect it)
aov(logRatio ~ Label * Experiment, data = m.cfu.all) %>% summary()
aov(Prop ~ Label * Experiment, data = m.cfu.all) %>% summary()

#Label:Experiment is non-significant, I remove interaction
#Check if Experiment affects significantly in log2Ratio
aov1 <- aov(logRatio ~ Label + Experiment, data = m.cfu.all) 
aov2 <- update(aov1, . ~ . - Experiment, data = m.cfu.all) 
anova(aov1, aov2) #p = .009073
#Check if Experiment affects significantly in proportion
aov3 <- aov(Prop ~ Label + Experiment, data = m.cfu.all) 
aov4 <- update(aov3, . ~ . - Experiment, data = m.cfu.all) 
anova(aov3, aov4) #p = .009857

####################################
# Ratio in competition for BOTH EXPERIMENTS ====
m.cfu.all$log2Ratio <- m.cfu.all$logRatio
m.cfu.all$Proportion <- m.cfu.all$Prop
#Use these to switch the graph between proportions and ratios
m.cfu.all$Prop <- m.cfu.all$log2Ratio
y.lab <- "log<sub>2</sub> WT:&Delta;*iol* ratio"

pal.ratio <- c(pal_material(palette = "brown")(10)[4],
               pal_material(palette = "brown")(10)[4],
               pal_material(palette = "green")(10)[6])
hline.exp1 <- mean(m.cfu.all$Prop[m.cfu.all$Experiment == "1" & m.cfu.all$Label == "Soil 0 dpi"])
hline.exp2 <- mean(m.cfu.all$Prop[m.cfu.all$Experiment == "2" & m.cfu.all$Label == "Soil 0 dpi"])

aov.plot <- aov(Prop ~ Label * Experiment, data = m.cfu.all) 
tuk <- TukeyHSD(aov.plot, which = "Label")

comps.ratio <- sapply(rownames(tuk$Label), strsplit, "-")
p.ratio <- tuk$Label[,"p adj"]
p.ratio <- ifelse(p.ratio < 1e-4, format(p.ratio, scientific = F, digits = 1), formatC(p.ratio, digits = 1))
p.format <- paste0(p.ratio,
                   rstatix::p_mark_significant(p.ratio) %>% gsub(x = ., "[^*]", ""))
pos.ratio <- sapply(comps.ratio, \(x) max(m.cfu.all$Prop[m.cfu.all$Label %in% x]))

mx <- ggplot(m.cfu.all, 
             aes(x = Label, y = Prop, fill = Label)) +
  ylab(y.lab) + 
  geom_jitter(aes(color = Label, shape = Experiment), width = .2) +
  geom_errorbar(data = m.cfu.all %>%
                  dplyr::group_by(Label) %>%
                  dplyr::mutate(dens = mean(Prop, na.rm = T),
                                std = sd(Prop)),
                aes(ymin = dens - std,
                    ymax = dens + std,
                    x = Label,
                    color = Label), 
                linewidth = 1, width = .1) +
  geom_hline(yintercept = hline.exp1, lty = "solid", color = "grey60") +
  geom_hline(yintercept = hline.exp2, lty = "dashed", color = "grey60") +
  stat_summary(aes(x = Label, 
                   y = Prop, 
                   color = Label,
                   shape = Experiment),
               fun = mean, geom = "point", size = 4, inherit.aes = F) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_markdown(),
        axis.title.y = element_markdown(),
        legend.position = "none",
        text = element_text(family = "Montserrat", size = 27),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray80")) +
  scale_x_discrete(labels = gen.labels) +
  scale_color_manual(values = pal.ratio) +
  geom_signif(comparisons = comps.ratio,
              annotations = p.format,
              y_position = pos.ratio * c(1, 1, .85),
              family = "Montserrat") +
  annotation_logticks(sides = "l", 
                      long = unit(4, "mm"), 
                      mid = unit(3, "mm"), 
                      short = unit(1, "mm"),
                      size = .2, scaled = T, base = 2)
mx

####################################
# CFU in competition for BOTH EXPERIMENTS ====

m.cfu.ind <- reshape2::melt(m.cfu.all[m.cfu.all$KO != 0, ], 
                            id = colnames(m.cfu.all)[!colnames(m.cfu.all) %in% c("WT", "KO")],
                            value.name = "cfu", variable.name = "Genotype")
m.cfu.ind$logCFU <- log10(m.cfu.ind$cfu)
m.cfu.ind <- m.cfu.ind[!is.infinite(m.cfu.ind$logCFU), ]
m.cfu.ind$Genotype %<>% sub(x = ., "KO", "iol")
m.cfu.ind$Genotype %<>% factor(., levels = c("WT", "iol"))

i.cfu1 <- read.csv("Second_exp_sample_sheet_cfu_individual.csv", sep = ";") 
i.cfu1$Time <- paste(i.cfu1$TP, "dpi") %>% sub(x = ., "7", "21")
i.cfu1$Label <- paste(i.cfu1$Comp, i.cfu1$Time)
i.cfu1 <- i.cfu1[i.cfu1$Label != "Root 0 dpi",]
i.cfu1$Genotype %<>% sub(x = ., "KO", "iol")
i.cfu1$Genotype %<>% factor(., levels = c("WT", "iol"))
i.cfu1$Label %<>% factor(., levels = c("Soil 0 dpi", "Soil 21 dpi", "Root 21 dpi"))
colnames(i.cfu1)[5] <- "CFU_g"
i.cfu.all <- rbind(data.frame(i.cfu1[, c("Label", "Genotype", "CFU_g")], 
                              Experiment = 1), 
                   data.frame(i.cfu[, c("Label", "Genotype", "CFU_g")], 
                              Experiment = 2)) 
colnames(i.cfu.all)[3] <- "cfu"
i.cfu.all$logCFU <- log10(i.cfu.all$cfu)
i.cfu.all$Experiment %<>% as.factor(.)

label.label <- sapply(strsplit(levels(i.cfu.all$Label) %>% as.character, " "), '[', 1)
names(label.label) <- levels(i.cfu.all$Label) %>% as.character

cfu.all <- rbind(data.frame(
  i.cfu.all[i.cfu.all$Genotype != "Control" &
              i.cfu.all$cfu != 0, ][, c("Label", "Experiment", "Genotype", "cfu")],
  Treatment = "Single"),
  data.frame(m.cfu.ind[, c("Label", "Experiment", "Genotype", "cfu")],
             Treatment = "Mix"
  ))
cfu.all$logCFU <- log10(cfu.all$cfu)
cfu.all$Genotype %<>% factor(., levels = c("Control", "WT", "iol"))
cfu.all$Treatment %<>% factor(., levels = c("Single",
                                            "Mix"))
label.label <- sapply(strsplit(levels(cfu.all$Label) %>% as.character, " "), '[', 1)
names(label.label) <- levels(cfu.all$Label) %>% as.character

# p-values for a and b plots ====
pvals <- vector(mode = "numeric", length = 0)
treats <- vector(mode = "character", length = 0)
labels <- vector(mode = "character", length = 0)
for (trt in unique(cfu.all$Treatment)) {
  for (lbl in c("Soil 0 dpi", "Soil 21 dpi", "Root 21 dpi")) {
    if (trt == "Single") {
      pvals <- c(pvals, aov(logCFU ~ Genotype + Experiment, 
                               data = cfu.all[cfu.all$Treatment == trt & 
                                                cfu.all$Label == lbl,]) %>% 
                   summary %>% unlist %>% '['("Pr(>F)1"))
    } else {
      pvals <- c(pvals, compare_means(data = m.cfu.ind[m.cfu.ind$Label == lbl,],
                                      paired = T, formula = logCFU ~ Genotype,
                                      method = "anova") %>% '[['("p.adj"))
    }
    treats <- c(treats, trt)
    labels <- c(labels, lbl)
  }
}

pvals <- paste0(ifelse(pvals < 1e-4, format(pvals, scientific = F, digits = 1), formatC(pvals, digits = 1)),
                rstatix::p_mark_significant(pvals) %>% gsub(x = ., "[^*]", "")) %>% sub(x = ., " ", "")
df.annot <- data.frame(Treatment = treats, 
                       Label = labels,
                       pval = pvals,
                       Genotype = "WT")
df.annot$Label %<>% factor(., c(levels = "Soil 0 dpi", "Soil 21 dpi", "Root 21 dpi"))

cfu.all$interaction <- interaction(cfu.all$Genotype, cfu.all$Treatment, cfu.all$Label)
cfu.groups <- glm.nb(cfu ~ interaction + Experiment, data = cfu.all[cfu.all$Label != "Soil 0 dpi",]) %>%
  glht(., linfct = mcp(interaction = "Tukey")) %>%
  cld(type = "response", decreasing = T) %>% 
  '[['("mcletters") %>% 
  '[['("Letters")
df.groups <- names(cfu.groups) %>% 
  strsplit(., ".", fixed = T) %>% 
  do.call(rbind, .) %>% 
  as.data.frame %>%
  cbind(., cfu.groups)
colnames(df.groups) <- c("Genotype", "Treatment", "Label", "Group")
df.groups$Label %<>% factor(., c(levels = "Soil 21 dpi", "Root 21 dpi"))

pal.counts <- c(pal_material("cyan")(10)[c(8, 4)]) 
# Plot b (mixed paired) ====
mi1 <- ggplot(m.cfu.ind[m.cfu.ind$Label != "Soil 0 dpi", ], 
              aes(x = Genotype, y = logCFU, fill = Genotype)) + 
  labs(y = "CFU per gram sample") + 
  geom_jitter(aes(color = Genotype, shape = Experiment), width = .2) +
  geom_line(aes(group = log2Ratio),
            alpha = 0.8,
            size = .2) +
  geom_errorbar(data = m.cfu.ind %>%
                  dplyr::filter(Label != "Soil 0 dpi") %>%
                  dplyr::group_by(Label, Genotype) %>%
                  dplyr::mutate(dens = mean(logCFU, na.rm = T),
                                std = sd(logCFU, na.rm = T)),
                aes(ymin = dens - std,
                    ymax = dens + std,
                    x = Genotype,
                    color = Genotype), 
                linewidth = 1, width = .1) +
  facet_wrap(~Label, labeller = labeller(Label = label.label)) +
  theme(
    axis.title = element_blank(), 
    axis.text.x = element_markdown(), 
    axis.title.y = element_markdown(),
    text = element_text(family = "Montserrat", size = 25),
    legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "gray95"),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "gray50"),
    strip.text = element_text(color = "white")
  ) +
  scale_x_discrete(labels = gen.labels) +
  scale_color_manual(values = pal.counts) +
  scale_y_continuous(labels = trans_format("c", math_format(10^.x)),
                     breaks = nice.breaker, limits = c(5, 8)) +
  geom_text(data = df.annot[df.annot$Treatment == "Mix" &
                              df.annot$Label != "Soil 0 dpi", ], 
            aes(x = "WT", y = 8, group = Treatment, label = pval),
            family = "Montserrat", size = 4, nudge_x = .5) +
  geom_text(data = df.groups[df.groups$Treatment == "Mix", ], 
            aes(x = Genotype, y = 7.7, group = Treatment, label = Group),
            family = "Montserrat", size = 4) +
  annotation_logticks(sides = "l", 
                      long = unit(4, "mm"), 
                      mid = unit(3, "mm"), 
                      short = unit(1, "mm"),
                      size = .2, scaled = T) +
  geom_segment(aes(x = "WT", xend = "iol", y = 7.9, yend = 7.9))
mi1



# Plot A (individual) ====
sn1 <- ggplot(i.cfu.all[i.cfu.all$Label != "Soil 0 dpi" & 
                          i.cfu.all$Genotype != "Control" &
                          i.cfu.all$cfu != 0, ], 
              aes(x = Genotype, y = logCFU, fill = Genotype)) + 
  labs(y = "CFU per gram sample") + 
  geom_jitter(aes(color = Genotype, shape = Experiment), width = .2) +
  geom_errorbar(data = i.cfu.all %>%
                  dplyr::filter(Label != "Soil 0 dpi") %>%
                  dplyr::filter(cfu != 0) %>%
                  dplyr::filter(Genotype != "Control") %>%
                  dplyr::group_by(Label, Genotype) %>%
                  dplyr::mutate(dens = mean(logCFU, na.rm = T),
                                std = sd(logCFU, na.rm = T)),
                aes(ymin = dens - std,
                    ymax = dens + std,
                    x = Genotype,
                    color = Genotype), 
                linewidth = 1, width = .1) +
  facet_wrap(~Label, labeller = labeller(Label = label.label)) +
  theme(
    axis.title = element_blank(), 
    axis.text.x = element_markdown(), 
    axis.title.y = element_markdown(),
    text = element_text(family = "Montserrat", size = 25),
    legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "gray95"),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "grey50"),
    strip.text = element_text(color = "white")
  ) +
  scale_x_discrete(labels = gen.labels) +
  scale_color_manual(values = pal.counts) +
  scale_y_continuous(labels = trans_format("c", math_format(10^.x)),
                     breaks = nice.breaker, limits = c(5, 8)) +
  geom_text(data = df.annot[df.annot$Treatment == "Single" &
                              df.annot$Label != "Soil 0 dpi", ], 
            aes(x = "WT", y = 8, group = Treatment, label = pval),
            family = "Montserrat", size = 4, nudge_x = .5) +
  geom_text(data = df.groups[df.groups$Treatment == "Single", ], 
            aes(x = Genotype, y = 7.7, group = Treatment, label = Group),
            family = "Montserrat", size = 4) +
  annotation_logticks(sides = "l", 
                      long = unit(4, "mm"), 
                      mid = unit(3, "mm"), 
                      short = unit(1, "mm"),
                      size = .2, scaled = T) +
  geom_segment(aes(x = "WT", xend = "iol", y = 7.9, yend = 7.9))

sn1

lay.cfu <- c("AAABBBCCC
             AAABBBCCC
             AAABBBCCC
             AAABBBCCC
             ")
cfu.full <- ((sn1 + theme(text = element_text(size = 20),
                          strip.text.x = element_text(color = "black", face = "bold"),
                          strip.background = element_rect(fill = "gray80"),
                          axis.title.y = element_markdown(margin = margin(r = 10))))) + 
  (mi1 + theme(legend.position = "right", 
               legend.key = element_rect(fill = "white"),
               legend.title = element_blank(),
               legend.text = element_text(size = 12),
               strip.text.x = element_text(color = "black", face = "bold"),
               strip.background = element_rect(fill = "gray80"),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(), 
               axis.line.y = element_blank(),
               text = element_text(size = 20)) +
     scale_shape_discrete(labels = c("1" = "Experiment 1", "2" = "Experiment 2")) +
     guides(shape = guide_legend(override.aes = list(size = 4)),
            color = "none", group = "none", fill = "none")) + 
  (mx + 
     scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 1)) +
     theme(text = element_text(size = 20),
           axis.text.x = element_markdown())) + 
  plot_layout(design = lay.cfu, guides = "collect") 
ragg::agg_jpeg("final_figure.jpg", width = 5000, height = 2000, res = 300)
cfu.full
invisible(dev.off())

ragg::agg_jpeg("final_figure_reduced.jpg", width = 4500, height = 1500, res = 300)
cfu.full
invisible(dev.off())





