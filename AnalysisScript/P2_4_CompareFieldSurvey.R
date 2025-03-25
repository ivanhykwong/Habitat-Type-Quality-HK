library(terra)

# load habitat quality map and field survey points
hq_2010 <- rast("D:/HabitatQuality/workspace/quality_c_y2010.tif")
pt_2010 <- read.csv("D:/HabitatQuality/validation/pt_2010_hq.csv")
pt_2010 <- vect(pt_2010, geom=c("Rect_X","Rect_Y"), crs=crs(hq_2010))

# extract pixel value
df <- as.data.frame(extract(hq_2010, pt_2010, method="bilinear", bind=TRUE, na.rm=TRUE))
df$Value3 <- as.numeric(factor(df$Value2, levels = c("Low","Medium","High")))
df <- df[!is.na(df$quality_c_test1),]

# Correlation and plot
cor(df$quality_c_test1, df$Value3, method = "kendall")
cor.test(df$quality_c_test1, df$Value3, method = "kendall")
plot(df$quality_c_test1, as.numeric(df$Value3))
boxplot(quality_c_test1 ~ Value2, data=df)

# Statistical test
kruskal.test(quality_c_test1 ~ Value2, data=df)
pairwise.wilcox.test(df$quality_c_test1, df$Value2, p.adjust.method = "BH")

# Boxplot

library(ggplot2)
library(ggtext)
library(dplyr)

setwd("D:/HabitatQuality/validation")
df$value_fact <- factor(df$Value2, levels=c("Low","Medium","High"))
p <- ggplot(df, aes(x=value_fact, y=quality_c_test1)) + 
  stat_boxplot(geom="errorbar", width = 0.3) +
  geom_boxplot(color="black", outlier.color="black") +
  scale_x_discrete(labels=c("Low<br>*(32)*","Medium<br>*(241)*","High<br>*(177)*"))+
  xlab("Ecological value assessed by ecologists in field") +
  ylab("Habitat quality from InVEST model") +
  theme_classic() +
  theme(axis.title=element_text(size=10),
        axis.text.x = element_markdown(lineheight=1.2, colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.title.y = element_text(margin=margin(r=5)),
        axis.title.x = element_text(margin=margin(t=5)))
p
ggsave("boxplot_validate2010.pdf",width=10,height=10,unit="cm")

df_stat <- df %>% group_by(value_fact) %>%
  summarise(min=min(quality_c_test1),
            q1=quantile(quality_c_test1,0.25),
            median=median(quality_c_test1),
            q3=quantile(quality_c_test1,0.75),
            max=max(quality_c_test1),
            mean=mean(quality_c_test1),
            sd=sd(quality_c_test1))
df_stat
