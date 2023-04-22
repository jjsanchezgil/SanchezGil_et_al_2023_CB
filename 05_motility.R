library(reshape2)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggpubr)
require(extrafont)
require(patchwork)
require(ggsci)
require(ggtext)
require(magick)
extrafont::font_import("/Windows/Fonts/", pattern = "Montserrat.*", prompt = F)
extrafont::loadfonts(device = "win")

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      q1 = unname(quantile(x[[col]], .25)),
      q2 = unname(quantile(x[[col]], .75)))
  }
  data_sum <- ddply(data, groupnames, .fun=summary_func,
                    varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

mixed_pal <- c(pal_futurama()(10)[c(9, 7, 3)],
               pal_material("orange")(10)[c(5, 10)]
               #pal_futurama()(10)[1]
)


plate_diam_mm <- 88
plate_diam_px <- 520

conversion <- plate_diam_mm / plate_diam_px

swim <- read.csv("swimming_data_1.csv", sep = ";")
swim$radius_px <- swim$radius
#Radius in mm
swim$radius <- swim$radius * conversion
swim$area_mm2 <- (swim$radius ^ 2) * pi
swim$plate <- paste0(swim$treatment, swim$batch, swim$repno)

swim$time <- ifelse(swim$timepoint == "48h", yes = 48, no = ifelse(swim$timepoint == "65h", yes = 65 - 48, no = 70 - 65))

swim.wi <- swim %>% filter(genotype %in% c("WT", "iol"))
swim.wi <- swim.wi[swim.wi$timepoint == "48h", ] %>% mutate(r48 = radius) %>% '['(-c(1, 2, 4, 9, 10, 11)) %>%
  left_join(swim.wi[swim.wi$timepoint == "65h", ] %>% mutate(r65 = radius) %>% '['(c("plate", "genotype", "r65")), by = c("plate", "genotype")) %>%
  left_join(swim.wi[swim.wi$timepoint == "70h", ] %>% mutate(r70 = radius) %>% '['(c("plate", "genotype", "r70")), by = c("plate", "genotype"))

swim.wi$t48 <- swim.wi$r48 / 48
swim.wi$t65 <- (swim.wi$r65 - swim.wi$r48) / (65 - 48)
swim.wi$t70 <- (swim.wi$r70 - swim.wi$r65) / (70 - 65)

#OUTPUT: there is a batch effect
# Df Sum Sq Mean Sq F value   Pr(>F)
# treatment                  4  43.96   10.99   4.618 0.001674 **
# batch                      3  42.13   14.04   5.900 0.000863 ***
# genotype                   1  85.96   85.96  36.115 2.05e-08 ***
# treatment:batch           12  37.27    3.11   1.305 0.224402
# treatment:genotype         4 131.48   32.87  13.809 2.71e-09 ***
# batch:genotype             3  20.47    6.82   2.867 0.039468 *
# treatment:batch:genotype  12  32.42    2.70   1.135 0.338667
# Residuals                120 285.62    2.38
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

aov(data = swim.wi[swim.wi$batch != "batch5", ], r48 ~ treatment * batch * genotype) %>% summary()
#OUTPUT: there is no batch effect if batch 5 is removed
# Df Sum Sq Mean Sq F value   Pr(>F)
# treatment                 4  49.28   12.32   5.857  0.00031 ***
# batch                     2   1.13    0.56   0.268  0.76579
# genotype                  1  84.21   84.21  40.031 9.51e-09 ***
# treatment:batch           8  18.08    2.26   1.075  0.38806
# treatment:genotype        4 138.31   34.58  16.438 3.82e-10 ***
# batch:genotype            2  15.21    7.60   3.615  0.03089 *
# treatment:batch:genotype  8  16.06    2.01   0.954  0.47680
# Residuals                90 189.32    2.10
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#
# #Swimming speed (mm/h) for the first 48h excluding batch 5
# ##splitting per batch
# ##WT
# ggplot(swim.wi[swim.wi$genotype == "WT", ],
#        aes(x = treatment,
#            y = t48)) +
#   geom_boxplot() + geom_jitter() + facet_wrap(~batch) +
#   stat_compare_means(ref.group = "Control", method = "t.test", label.y.npc = 0.95)
# ##iol
# ggplot(swim.wi[swim.wi$genotype == "iol", ],
#        aes(x = treatment,
#            y = t48)) +
#   geom_boxplot() + geom_jitter() + facet_wrap(~batch) +
#   stat_compare_means(ref.group = "Control", method = "t.test", label.y.npc = 0.95)
#
# #Swimming speed (mm/h) for the first 48h excluding batch 5
# ##all batches merged
# ##WT
# ggplot(swim.wi[swim.wi$genotype == "WT" & swim.wi$batch != "batch5", ],
#        aes(x = treatment,
#            y = t48)) +
#   geom_boxplot() + geom_jitter() +
#   stat_compare_means(ref.group = "Control", method = "t.test")
#
# ##iol
# ggplot(swim.wi[swim.wi$genotype == "iol" & swim.wi$batch != "batch5", ],
#        aes(x = treatment,
#            y = t48)) +
#   geom_boxplot() + geom_jitter() +
#   stat_compare_means(ref.group = "Control", method = "t.test")
#
# #Swimming speed (mm/h) for 48-65h excluding batch 5
# ##all batches merged
# ##WT
# ggplot(swim.wi[swim.wi$genotype == "WT" & swim.wi$batch != "batch5", ],
#        aes(x = treatment,
#            y = t65)) +
#   geom_boxplot() + geom_jitter() +
#   stat_compare_means(ref.group = "Control", method = "t.test")
#
# ##iol
# ggplot(swim.wi[swim.wi$genotype == "iol" & swim.wi$batch != "batch5", ],
#        aes(x = treatment,
#            y = t65)) +
#   geom_boxplot() + geom_jitter() +
#   stat_compare_means(ref.group = "Control", method = "t.test")
#
# #Swimming speed (mm/h) for 65-70h excluding batch 5
# ##all batches merged
# ##WT
# ggplot(swim.wi[swim.wi$genotype == "WT" & swim.wi$batch != "batch5", ],
#        aes(x = treatment,
#            y = t70)) +
#   geom_boxplot() + geom_jitter() +
#   stat_compare_means(ref.group = "Control", method = "t.test")
#
# ##iol
# ggplot(swim.wi[swim.wi$genotype == "iol" & swim.wi$batch != "batch5", ],
#        aes(x = treatment,
#            y = t70)) +
#   geom_boxplot() + geom_jitter() +
#   stat_compare_means(ref.group = "Control", method = "t.test")


#swim.sp <- reshape2::melt(swim.wi, id.vars = 1:7, value.name = "speed", measure.vars = c("t48", "t65", "t70"), variable.name = "stage")

#summ.sp <- data_summary(swim.sp[swim.sp$batch != "batch5", ], varname = "speed", groupnames = c("treatment", "genotype", "stage"))
# ggplot(summ.sp,
#        aes(x = stage, y = speed, color = treatment, group = treatment, shape = genotype)) +
#   geom_pointrange(aes(ymin = q1, ymax = q2), fatten = 8) +
#   geom_line() +
#   facet_wrap(~genotype)
#

swim.wi <- swim.wi[swim.wi$batch != "batch5", ]

flabs <- ifelse(swim.wi$genotype == "iol", yes = "&Delta;*iol*", no = "WT")
names(flabs) <- swim.wi$genotype
xlabs <- paste0(sapply(swim.wi$treatment, \(x) ifelse(startsWith(x, "Glc"), "Glucose", ifelse(startsWith(x, "Ino"), "Inositol", "Control"))),
                "<br>",
                sapply(swim.wi$treatment, \(x) ifelse(endsWith(x, "1"), "1 mM", ifelse(endsWith(x, "10"), "10 mM", "Control")))
) %>% sub(x = ., pattern = "Control<br>", "")
names(xlabs) <- swim.wi$treatment

swim.wi$genotype <- as.factor(swim.wi$genotype) %>% relevel(., "WT")
#Boxplots showing differences per treatment in each genotype
sw1 <- ggplot(swim.wi[swim.wi$batch != "batch5" &
                        swim.wi$plate != "Ino10batch33", ],
              aes(x = treatment,
                  y = r48,
                  fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = .5, lwd = .2) +
  geom_jitter(position = position_jitterdodge(),
              size = 1,
              aes(color = treatment)) +
  stat_compare_means(ref.group = "Control", method = "t.test", label = "p.signif", hide.ns = TRUE, size = 7, label.y.npc = .95) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_markdown(lineheight = 0.4),
        axis.title.y = element_markdown(lineheight = 0.5),
        strip.text.x = element_markdown(color = "black",
                                        lineheight = .25,
                                        face = "bold"),
        legend.position = "none",
        text = element_text(family = "Montserrat", size = 18),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray80")) +
  facet_wrap(~genotype, labeller = labeller(genotype = flabs)) +
  scale_x_discrete(labels = xlabs) +
  scale_color_manual(values = mixed_pal) + scale_fill_manual(values = mixed_pal) +
  ylab("Swimming radius (mm)") +
  ylim(c(4, 14))


sw2 <- ggplot(swim.wi[swim.wi$batch != "batch5", ],
              aes(x = genotype, y = r48, fill = treatment)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, lwd = .2) +
  geom_point(aes(color = treatment), size = 1) +
  geom_line(aes(group = plate, color = treatment),
            alpha = 0.8,
            size = .2) +
  stat_compare_means(paired = TRUE, aes(group = genotype),
                     method = "t.test",
                     label = "p.signif",
                     hide.ns = TRUE,
                     family = "Montserrat",
                     size = 7,
                     label.y.npc = .8,
                     label.x.npc = .4) +
  facet_wrap(~treatment, labeller = labeller(treatment = xlabs), nrow = 1) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_markdown(lineheight = 0.5),
        axis.title.y = element_markdown(lineheight = 0.5),
        strip.text.x = element_markdown(color = "black",
                                        lineheight = 0.4,
                                        face = "bold"),
        legend.position = "none",
        text = element_text(family = "Montserrat", size = 18),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray80")) +
  scale_x_discrete(labels = flabs) +
  scale_color_manual(values = mixed_pal) + scale_fill_manual(values = mixed_pal) +
  ylab("Swimming<br>radius (mm)")

l1 <- c("
        AAAA
        AAAA
        AAAA
        AAAA
        AAAA
        BBBB
        BBBB
        ")
sw.plot <- sw1 / sw2 + plot_layout(design = l1) +
  plot_annotation(tag_levels = "a")
sw.plot

sw.pictures <- magick::image_ggplot(image = magick::image_read("full_legend_swimming.jpg", density = 300)) +
  theme(text = element_text(family = "Montserrat", size = 18))

l2 <- c("
        AAAAAAACC
        AAAAAAACC
        AAAAAAACC
        AAAAAAACC
        AAAAAAACC
        AAAAAAACC
        BBBBBBBBB
        BBBBBBBBB
        ")

sw.full <- sw1 + sw2 + sw.pictures + 
  plot_layout(design = l2) + 
  plot_annotation(tag_levels = "a")

ragg::agg_jpeg("swimming.jpg", width = 5000, height = 2300, res = 300)
sw.full
invisible(dev.off())

#Exp 2 swimming ====
ds <- read.csv("swimming_data_2.csv")
comps <- as.data.frame(compare_means(Area_cm2 ~ Treatment, data = ds[ds$Genotype == "WT",]))
comps <- comps[comps$p.format < 0.05 , c("group1", "group2", "p.format")]
trans <- sapply(unique(ds$Treatment), as.character)
comps$group1 <- sapply(comps$group1, function(x) which(trans == x))
comps$group2 <- sapply(comps$group2, function(x) which(trans == x))
comps$ord <- comps$group2 - comps$group1
comps <- comps[order(comps$group2, -comps$group1),]
comps <- lapply(1:nrow(comps), function(x) c(comps[x, "group1"], comps[x, "group2"]))
npg_pal <- pal_futurama()(10)[c(9, 3, 7, 1, 9, 9, 9)]

ggplot(data = ds, aes(x = Treatment, y = Area_cm2, fill = Treatment)) + 
  geom_boxplot(alpha = .5,outlier.shape = NA) + #facet_grid(.~Genotype) + 
  facet_wrap(~Genotype) + #Facets by Wash
  stat_compare_means(aes(x = Treatment, y = Area_cm2, group = Treatment),
                     method = "anova", size = 7, hide.ns = FALSE, 
                     label = "p.format", family = "montserrat",
                     label.y = 3.5, vjust = .1) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_markdown(), 
        axis.title.y = element_markdown(),
        strip.text.x = element_markdown(color = "white"),
        legend.position = "none",
        text = element_text(family = "montserrat", size = 20),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray30")) +
  geom_jitter(size = 1, aes(color = Treatment)) + 
  scale_color_manual(values = npg_pal) + scale_fill_manual(values = npg_pal) +
  ylab("Colony surface (cm<sup>2</sup>)") 

#Swarming ====

mixed_pal <- c(pal_futurama()(10)[c(9, 7, 3)],
               pal_material("orange")(10)[c(5, 10)]
)

sw <- read.csv("swarming.csv", sep = ";") 
sw$Area_mm2 %<>% sub(",", ".", .) %>% as.numeric()
sw$Area_cm2 %<>% sub(",", ".", .) %>% as.numeric()
sw <- sw[sw$Genotype %in% c("WT", "iol"),]
sw$Radius <- sqrt(sw$Area_mm2 / pi)
sw$Genotype <- factor(sw$Genotype, levels = c("WT", "iol"))

flabs <- ifelse(sw$Genotype == "iol", yes = "Δ*iol*", no = "WT")
names(flabs) <- sw$Genotype
xlabs <- paste0(sapply(sw$Treatment, \(x) ifelse(startsWith(x, "Glc"), "Glucose", ifelse(startsWith(x, "Ino"), "Inositol", "Control"))),
                "<br>",
                sapply(sw$Treatment, \(x) ifelse(endsWith(x, "1 mM"), "1 mM", ifelse(endsWith(x, "10 mM"), "10 mM", "Control")))
) %>% sub(x = ., pattern = "Control<br>", "")
names(xlabs) <- sw$Treatment


#X axis is treatment
ggplot(sw[sw$Genotype %in% c("WT", "iol"), ],
       aes(x = Treatment, y = Area_cm2, fill = Genotype)) + 
  geom_boxplot() +
  geom_jitter() + 
  facet_wrap(~Genotype, nrow = 1) +
  stat_compare_means(aes(group = Treatment), ref.group = "Control", label = "p.signif", size = 7) +
  ylab("Swarming area (mm<sup>2</sup>)") +
  theme(axis.title.y = element_markdown(),
        axis.title.x = element_blank(), 
        axis.text.x = element_markdown(), 
        axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_markdown(),
        text = element_text(family = "montserrat", size = 20),
        panel.background = element_rect(fill = "white"),
        strip.text = element_text(colour = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray30"))

swrm1 <- ggplot(sw[sw$Genotype %in% c("WT", "iol"), ],
                aes(x = Treatment, y = Area_cm2, fill = Treatment)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, lwd = .2) +
  geom_jitter(position = position_jitterdodge(),
              size = .3,
              aes(color = Treatment)) +
  stat_compare_means(ref.group = "Control", method = "t.test", label = "p.signif", hide.ns = TRUE, size = 10, label.y.npc = .95) +
  ylab("Swarm surface (mm<sup>2</sup>)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_markdown(lineheight = 0.4),
        axis.title.y = element_markdown(lineheight = 0.5),
        strip.text.x = element_markdown(color = "black",
                                        lineheight = .25,
                                        face = "bold"),
        legend.position = "none",
        text = element_text(family = "montserrat", size = 27),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray80")) +
  facet_wrap(~Genotype, labeller = labeller(Genotype = flabs), nrow = 1) +
  scale_x_discrete(labels = xlabs) +
  scale_color_manual(values = mixed_pal) + scale_fill_manual(values = mixed_pal)

#X axis is genotype
swrm2 <- ggplot(sw[sw$Genotype %in% c("WT", "iol"), ],
                aes(x = Genotype, y = Area_cm2, fill = Treatment)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, lwd = .2) +
  geom_point(aes(color = Treatment), size = .3) + 
  geom_line(aes(group = Label, color = Treatment),
            alpha = 0.8,
            size = .2) +
  stat_compare_means(paired = TRUE, aes(group = Genotype),
                     method = "t.test",
                     label = "p.signif",
                     hide.ns = TRUE,
                     family = "montserrat",
                     size = 12,
                     label.y.npc = .8,
                     label.x.npc = .4) +
  ylab("Swarm<br>surface (mm<sup>2</sup>)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_markdown(lineheight = 0.4),
        axis.title.y = element_markdown(lineheight = 0.5),
        strip.text.x = element_markdown(color = "black",
                                        lineheight = .25,
                                        face = "bold"),
        legend.position = "none",
        text = element_text(family = "montserrat", size = 27),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "gray80")) +
  facet_wrap(~Treatment, labeller = labeller(Treatment = xlabs), nrow = 1) +
  scale_x_discrete(labels = flabs) +
  scale_color_manual(values = mixed_pal) + scale_fill_manual(values = mixed_pal)

swrm.l1 <- c("
        AAAA
        AAAA
        AAAA
        AAAA
        AAAA
        BBBB
        BBBB
        ")
swrm.plot <- swrm1 / swrm2 + plot_layout(design = swrm.l1) +
  plot_annotation(tag_levels = "a")
swrm.plot

