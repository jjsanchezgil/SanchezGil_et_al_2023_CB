require(ggplot2, quietly = TRUE)
require(ggtext, quietly = TRUE)
require(ggsci, quietly = TRUE)
require(ggrepel, quietly = TRUE)
require(extrafont, quietly = TRUE)
require(readxl, quietly = TRUE)

extrafont::font_import("/Windows/Fonts/", pattern = "Montserrat.*", prompt = F)
extrafont::loadfonts(device = "win")

p.pal <- c(pal_material("blue")(6)[6],
           pal_material("red")(6)[6])

ft <- read_excel("cole_2017_data_iol.xlsx")
ft <- ft[ft$Gene.Symbol!= "mcp", ]
ft$PlantFit <- -ft$RPL.t
ft$PlantSig <- -log10(ft$RPL.p)
ft$NmeshFit <- -ft$NRF.mean
ft$NmeshSig <- -log10(ft$NRF.p)
ft$PMSig <- ft$RPL.NRF.p
gf1 <- ggplot(ft, aes(x = PlantFit, y = PlantSig)) +
  ylab("-log<sub>10</sub> *p*-value (FDR-corrected)") +
  ylim(c(0, 6)) +
  xlab("Contribution to root colonisation") +
  xlim(c(-1, 30)) +
  labs(color = "-log<sub>10</sub>*p*") +
  geom_point(aes(color = PlantSig)) +
  geom_text_repel(aes(label = Gene.Symbol),
                  family = "Montserrat",
                  fontface = "italic") +
  geom_hline(yintercept = -log10(.05),
             lty = "dashed", color = "gray60") +
  geom_text(aes(x = 28, y = -log10(.05) + .2, label = "italic(p)==.05"),
            family = "Montserrat",
            parse = TRUE,
            fontface = "italic") +
  geom_hline(yintercept = -log10(.001),
             lty = "dashed", color = "gray60") +
  geom_text(aes(x = 28, y = -log10(.001) + .2, label = "italic(p)==.001"),
            family = "Montserrat",
            parse = TRUE,
            fontface = "italic") +
  scale_color_gradient(low = p.pal[1], high = p.pal[2]) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_markdown(),
        strip.text.x = element_markdown(color = "black",
                                        face = "bold"),
        legend.title = element_markdown(),
        legend.position = "bottom",
        text = element_text(family = "Montserrat", size = 18),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray95"),
        axis.line = element_line(color = "black"))

ragg::agg_jpeg("cole_figure_7d.jpg", width = 1500, height = 2000, res = 300)
gf1
invisible(dev.off())

