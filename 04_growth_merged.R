library(reshape2)
library(magrittr)
library(ggplot2)
library(ggpubr)
require(showtext)
require(patchwork)
require(ggsci)
require(ggtext)
require(dplyr)
require(emmeans)
require(multcomp)


extrafont::font_import("/Windows/Fonts/", pattern = "Montserrat.*", prompt = F)
extrafont::loadfonts(device = "win")

pal.growth <- c(pal_material("blue-grey")(5)[c(4)],
                pal_material("blue")(5)[c(4)],
                pal_material("orange")(5)[c(4)]
)

gr1 <- readRDS("growth_28_7.RDS")
gr2 <- readRDS("growth_6_2.RDS")

cols <- c("Time", "cell", "OD", "genotype", 
          "treatment", "pvd", "pch", 
          "ginc", "pinc", "cinc")

gr.full <- rbind(data.frame(gr1[gr1$nut == "1X", cols],
                            Experiment = 1),
                 data.frame(gr2[, cols],
                            Experiment = 2))
gr.full$treatment %<>% 
  sub(x = ., "Glucose 10 mM", "Glc") %>%
  sub(x = ., "Inositol 10 mM", "Ino") %>%
  factor(., levels = c("Control", "Glc", "Ino"))
gr.full$Experiment %<>% factor(.)
gr.full <- merge(gr.full, gr.full %>%
        dplyr::filter(Time == "0") %>%
        dplyr::group_by(Experiment, treatment,genotype) %>%
                   dplyr::summarise(mOD = mean(OD)),
                 by = c("Experiment", "treatment", "genotype"))
gr.full <- merge(gr.full, gr.full %>%
                   dplyr::filter(Time == "0") %>%
                   dplyr::group_by(Experiment, treatment, genotype) %>%
                   dplyr::summarise(med.pvd = mean(pvd)),
                 by = c("Experiment", "treatment", "genotype"))

gr.full$nOD <- gr.full$OD - gr.full$mOD
gr.full$npvd <- gr.full$pvd - gr.full$med.pvd


gr.full$cOD <- NA
blks <- gr.full %>% dplyr::group_by(Time, Experiment, genotype) %>% dplyr::filter(treatment == "Control") %>% dplyr::summarise(bOD = median(nOD))
for (i in 1:nrow(gr.full)) {
  gr.full[i, "cOD"] <- gr.full[i, "nOD"] - blks[blks$genotype == gr.full[i, "genotype"] & blks$Time == gr.full[i, "Time"] & blks$Experiment == gr.full[i, "Experiment"], "bOD", drop = TRUE]
}

gr.full %<>%
  dplyr::group_by(Experiment, genotype, treatment, Time) %>%
  dplyr::mutate(q90 = quantile(npvd, .9)) %>%
  dplyr::filter(npvd < q90) %>%
  as.data.frame


#Growth in 1X CC ====
gcc1 <- ggplot(gr.full[gr.full$Time %in% seq(0, 48, .5),] %>%
         dplyr::group_by(cell, Experiment) %>%
         dplyr::arrange(Time, .by_group = T) %>%
         dplyr::mutate(nOD = zoo::rollapply(nOD, mean, width = 6, align = "left", fill = NA)) %>%
         dplyr::filter(Time %in% seq(0, 48, 1)),# / time.inc),
       aes(x = as.factor(Time),
           y = nOD,
           fill = treatment,
           color = treatment)) +
  # ylim(c(-.03,.06)) + 
  # xlim(c(20, 24)) +
  geom_boxplot(alpha = .2, 
               outlier.shape = NA,
               linewidth = .3,
               show.legend = F,
               position = "identity") +
  geom_smooth(aes(group = treatment),
              method = "loess",
              span = .1, show.legend = F) +
  facet_wrap(~genotype) +
  ylab("Growth in CC (OD<sub>600</sub>)") +
  xlab("Time (h)") +
  geom_hline(yintercept = 0) +
  geom_point(size = NA, shape = 19) +
  scale_x_discrete(breaks = seq(0, 48, 12)[-5] %>% as.character()) +
  scale_fill_manual(values = pal.growth) +
  scale_color_manual(values = pal.growth,
                     labels = c("Control" = "Control",
                                "Glc" = "Glucose",
                                "Ino" = "Inositol")) +
  guides(fill = "none",
         color = guide_legend(override.aes = list(fill = NA,
                                                  shape = 19,
                                                  size = 4.5,
                                                  alpha = 1))) +
  theme(axis.title.y = element_markdown(),
        axis.text.x = element_markdown(),
        axis.text.y = element_markdown(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        text = element_text(family = "Montserrat", size = 18),
        panel.background = element_rect(fill = "white"),
        strip.text = element_markdown(colour = "black", face = "bold"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray80"))

ragg::agg_jpeg("growth_cc_merged.jpg", width = 2200, height = 2000, res = 300)
gcc1
invisible(dev.off())

#Growth curve analysis CC ====
library(growthcurver)

gcv <- gr.full[gr.full$Time %in% seq(0, 48, .5),] %>% 
  .[complete.cases(.),] %>% 
  dplyr::group_by(Experiment, treatment, genotype, cell) %>% 
  dplyr::summarise(r = growthcurver::SummarizeGrowth(Time, nOD)$vals$r, 
                   k = growthcurver::SummarizeGrowth(Time, nOD)$vals$k, 
                   s= growthcurver::SummarizeGrowth(Time, nOD)$vals$sigma)

gcv.control <- gcv %>% 
  dplyr::filter(treatment == "Control") %>% 
  dplyr::group_by(Experiment) %>% 
  dplyr::summarise(mr = median(r), 
                   mk = median(k), 
                   ms = median(s)) %>%
  as.data.frame
gcv$nr <- sapply(unique(gcv$Experiment), \(x) gcv$r[gcv$Experiment == x] - gcv.control[gcv.control$Experiment == x, "mr"]) %>% c
gcv$nk <- sapply(unique(gcv$Experiment), \(x) gcv$k[gcv$Experiment == x] - gcv.control[gcv.control$Experiment == x, "mk"]) %>% c
gcv$ns <- sapply(unique(gcv$Experiment), \(x) gcv$s[gcv$Experiment == x] - gcv.control[gcv.control$Experiment == x, "ms"]) %>% c
gcv <- gcv[gcv$k < 1, ]

m1 <- glm(nk ~ treatment * genotype * Experiment, data = gcv)
summary(m1) #Experiment doesn't seem significant per se
m2 <- glm(nk ~ treatment * genotype, data = gcv)
summary(m2)

gcv$int <- interaction(gcv$treatment, gcv$genotype)
groups <- glm(nk ~ int, data = gcv) %>% 
  emmeans::emmeans(., specs = "int") %>%
  cld(., decreasing = T, alpha = .05, adjust = "sidak") %>% 
  as.data.frame %>%
  '['(c("int", ".group", "upper.CL"))
#
mx.nm <- groups$.group %>% as.numeric %>% as.character %>% 
  strsplit(., "") %>% sapply(., max) %>% max
ngroups <- groups$.group %>% as.numeric %>% as.character %>% 
  strsplit(., "") %>% sapply(\(x) as.numeric(mx.nm) - as.numeric(x) + 1) %>% 
  sapply(\(x) letters[sort(x)]) %>% sapply(paste0, collapse = "")
groups$ngroups <- ngroups
df.labels <- data.frame(
  genotype = groups$int %>% as.character %>% strsplit("\\.") %>% sapply('[', 2),
  treatment = groups$int %>% as.character %>% strsplit("\\.") %>% sapply('[', 1),
  Group = ngroups,
  nk = groups$upper.CL
)
df.labels$genotype %<>% factor(., levels = c("WT", "&Delta;*iol*"))

gcv1 <- ggplot(gcv, aes(x = treatment, y = nk, color = treatment)) +
  geom_errorbar(data = gcv %>%
                  dplyr::group_by(treatment, genotype) %>%
                  dplyr::mutate(dens = mean(nk, na.rm = T),
                                std = sd(nk)),
                aes(ymin = dens - std,
                    ymax = dens + std,
                    x = treatment,
                    color = treatment), 
                linewidth = 1, width = .1,
                position = position_dodge(.25)) +
  geom_jitter(aes(color = treatment, shape = Experiment),
              width = .2) +
  facet_wrap(~genotype) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_text(data = df.labels, 
            aes(x = treatment, y = nk + .05, 
                group = genotype, label = Group),
            color = "black",
            position = position_dodge(.75),
            family = "Montserrat",
            size  = 5) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = .9),
        axis.title.y = element_markdown(),
        legend.position = "none",
        text = element_text(family = "Montserrat", size = 18),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.text.x = element_markdown(color = "black", face = "bold"),
        strip.background = element_rect(fill = "gray80")) +
  scale_color_manual(values = pal.growth) + 
  scale_x_discrete(labels = c("Control" = "Control",
                              "Glc" = "Glucose",
                              "Ino" = "Inositol")) +
  ylab("*K* vs control")

gcv1

#Growth in M9 ----

m9 <- readxl::read_excel("growth_M9_20_9_22.xlsx")
m9 <- melt(m9, id.vars = c(1, 2))

m9$col <- sapply(m9$variable, stringr::str_remove, "[A-Z]") %>% as.numeric()
m9$row <- sapply(m9$variable, stringr::str_remove, "[0-9]+")
m9 <- m9[m9$col <= 6,]

m9$genotype <- ifelse(m9$col %in% c(1, 3, 5), yes = "WT", no = "KO")

treat <- c(rep("Control", 2), rep("Glc 100 mM", 2), rep("Ino 100 mM", 2))
names(treat) <- 1:6
m9$treatment <- treat[as.character(m9$col)]
m9$Time <- rep(rep(seq(0, 80, 0.5), 8), 6)

m9$nOD <- NA
blks <- m9 %>% dplyr::group_by(Time) %>% dplyr::filter(row == "H") %>% dplyr::summarise(bOD = median(value))
for (i in 1:nrow(m9)) {
  m9[i, "nOD"] <- m9[i, "value"] - blks[blks$Time == m9[i, "Time"], "bOD", drop = TRUE]
}

m9$genotype <- ifelse(m9$genotype == "KO", yes = "&Delta;*iol*", no = "WT") %>%
  factor(., levels = c("WT", "&Delta;*iol*"))

m9.plot <- m9[m9$row != "H" &
                m9$Time %in% seq(0, 72, 1),] %>%
  dplyr::group_by(treatment, Time, genotype) %>%
  dplyr::mutate(mOD = mean(nOD), sOD = sd(nOD)) %>%
  dplyr::filter((nOD < (mOD + 2*sOD)),(nOD > (mOD - 2*sOD)))

desic.lm <- lm(nOD ~ I(.01 + Time), m9.plot[m9.plot$treatment == "Control",])
m9.plot$cOD <- m9.plot$nOD - predict(desic.lm, m9.plot) + .01

g9 <- ggplot(m9.plot %>%
               dplyr::filter(Time %in% seq(0, 72, 6)),
             aes(x = Time,
                 y = cOD,
                 fill = treatment,
                 color = treatment)) +
  geom_jitter(aes(x = Time, color = treatment), alpha = .3) +
  geom_errorbar(data = m9.plot %>%
                  dplyr::group_by(Time, treatment, genotype) %>%
                  dplyr::mutate(mOD = mean(cOD),
                                std = sd(cOD)),
                aes(x = Time, color = treatment,
                    ymin = mOD - std,
                    ymax = mOD + std),
                linewidth = .5, alpha = .7, width = 3) +
  geom_smooth(aes(group = treatment),
              method = "gam", level = .99,
              alpha = .2) +
  facet_wrap(~genotype) +
  ylab("Growth in M9 medium (OD<sub>600</sub>)") +
  xlab("Time (h)") +
  scale_x_continuous(breaks = seq(0, 72, 24)) +
  scale_y_continuous(limits = c(0.005, .03)) +
  scale_fill_manual(values = pal.growth,
                    labels = c("Control" = "Control",
                               "Glc 100 mM" = "Glucose",
                               "Ino 100 mM" = "Inositol")) +
  scale_color_manual(values = pal.growth,
                     labels = c("Control" = "Control",
                                "Glc 100 mM" = "Glucose",
                                "Ino 100 mM" = "Inositol")) +
  guides(fill = "none",
         color = guide_legend(override.aes = list(fill = NA,
                                                  shape = 19,
                                                  size = 4.5,
                                                  alpha = 1))) +
  theme(axis.title.x = element_markdown(vjust = -2),
        axis.text.x = element_markdown(),
        axis.title.y = element_markdown(),
        legend.key = element_blank(),
        legend.key.height = unit(0, "mm"),
        legend.key.width = unit(0, "mm"),
        legend.key.size = unit(0, "mm"),
        legend.text = element_markdown(),
        legend.title = element_blank(),
        text = element_text(family = "Montserrat", size = 18),
        panel.background = element_rect(fill = "white"),
        strip.text = element_markdown(colour = "black", face = "bold"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray80")) 
g9

#Growth curve analysis M9 ====
library(growthcurver)

m9.plot <- m9[m9$row != "H" &
                m9$Time %in% seq(0, 72, .5),] %>%
  dplyr::group_by(treatment, Time, genotype) %>%
  dplyr::mutate(mOD = mean(nOD), sOD = sd(nOD)) %>%
  dplyr::filter((nOD < (mOD + 2*sOD)),(nOD > (mOD - 2*sOD)))

desic.lm <- lm(nOD ~ I(.01 + Time), m9.plot[m9.plot$treatment == "Control",])
m9.plot$cOD <- m9.plot$nOD - predict(desic.lm, m9.plot) + .01

m9.plot$cell <- paste0(m9.plot$row, m9.plot$col) %>% factor(.)
#m9o <- m9.plot
#m9.plot <- m9o
m9.plot %<>% 
  dplyr::group_by(treatment, genotype, cell) %>% 
  dplyr::mutate(n = dplyr::n()) %>%
  dplyr::filter(n > 10) %>%
  dplyr::arrange(Time, .by_group = T) %>%
  dplyr::summarise(r = growthcurver::SummarizeGrowth(Time, cOD)$vals$r, 
                   k = growthcurver::SummarizeGrowth(Time, cOD)$vals$k, 
                   kp = growthcurver::SummarizeGrowth(Time, cOD)$vals$k_p, 
                   rp = growthcurver::SummarizeGrowth(Time, cOD)$vals$r_p, 
                   s= growthcurver::SummarizeGrowth(Time, cOD)$vals$sigma) %>%
  dplyr::filter(k < 1)# & rp < .05)

m9.control <- m9.plot %>% 
  dplyr::filter(treatment == "Control") %>% 
  dplyr::group_by(genotype) %>%
  dplyr::filter(kp < .05) %>%
  dplyr::summarise(mr = mean(r), 
                   mk = mean(k), 
                   ms = mean(s)) %>%
  as.data.frame
m9.plot$nr <- NA
m9.plot$nk <- NA
m9.plot$ns <- NA
for(r in 1:nrow(m9.plot)) {g <- m9.plot[r, "genotype", drop = T] %>% as.character; m9.plot[r, "nr"] <- m9.plot[r, "r", drop = T] - m9.control[m9.control$genotype == g, "mr"]}
for(r in 1:nrow(m9.plot)) {g <- m9.plot[r, "genotype", drop = T] %>% as.character; m9.plot[r, "nk"] <- m9.plot[r, "k", drop = T] - m9.control[m9.control$genotype == g, "mk"]}
for(r in 1:nrow(m9.plot)) {g <- m9.plot[r, "genotype", drop = T] %>% as.character; m9.plot[r, "ns"] <- m9.plot[r, "s", drop = T] - m9.control[m9.control$genotype == g, "ms"]}

m2 <- glm(nk ~ treatment * genotype, data = m9.plot)
summary(m2)

m9.plot$int <- interaction(m9.plot$treatment, m9.plot$genotype)
groups <- glm(nk ~ int, data = m9.plot) %>% 
  emmeans::emmeans(., specs = "int") %>%
  cld(., decreasing = T, alpha = .05, adjust = "sidak") %>% 
  as.data.frame %>%
  '['(c("int", ".group", "upper.CL"))
#
mx.nm <- groups$.group %>% as.numeric %>% as.character %>% 
  strsplit(., "") %>% sapply(., max) %>% max
ngroups <- groups$.group %>% as.numeric %>% as.character %>% 
  strsplit(., "") %>% sapply(\(x) as.numeric(mx.nm) - as.numeric(x) + 1) %>% 
  sapply(\(x) letters[sort(x)]) %>% sapply(paste0, collapse = "")
groups$ngroups <- ngroups
df.labels <- data.frame(
  genotype = groups$int %>% as.character %>% strsplit("\\.") %>% sapply('[', 2),
  treatment = groups$int %>% as.character %>% strsplit("\\.") %>% sapply('[', 1),
  Group = ngroups,
  nk = groups$upper.CL
)
df.labels$genotype %<>% factor(., levels = c("WT", "&Delta;*iol*"))

gcm9 <- ggplot(m9.plot, aes(x = treatment, y = nk, color = treatment)) +
  geom_errorbar(data = m9.plot %>%
                  dplyr::group_by(treatment, genotype) %>%
                  dplyr::mutate(dens = mean(nk, na.rm = T),
                                std = sd(nk)),
                aes(ymin = dens - std,
                    ymax = dens + std,
                    x = treatment,
                    color = treatment), 
                linewidth = 1, width = .1,
                position = position_dodge(.25)) +
  geom_jitter(aes(color = treatment),
              width = .2) +
  facet_wrap(~genotype) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_text(data = df.labels, 
            aes(x = treatment, y = nk + .001, 
                group = genotype, label = Group),
            color = "black",
            position = position_dodge(.75),
            family = "Montserrat",
            size  = 5) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = .9),
        axis.title.y = element_markdown(),
        legend.position = "none",
        text = element_text(family = "Montserrat", size = 18),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.text.x = element_markdown(color = "black", face = "bold"),
        strip.background = element_rect(fill = "gray80")) +
  scale_color_manual(values = pal.growth) + 
  scale_x_discrete(labels = c("Control" = "Control",
                              "Glc 100 mM" = "Glucose",
                              "Ino 100 mM" = "Inositol")) +
  ylab("*K* vs control")

gcm9



#Get g9 from growth_M9.R
ragg::agg_jpeg("growth_all_merged.jpg", width = 5000, height = 2000, res = 300)
(g9 + theme(legend.position = "none") + gcm9) + gcc1 + gcv1 +
  plot_layout(guides = "collect", nrow = 1) + plot_annotation(tag_levels = "a")
invisible(dev.off())


