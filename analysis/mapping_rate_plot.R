# plot the mapping rate of CAPS
# Zhiyuan Hu
# 4 July 2020

library(ggplot2)
df_plot <- data.frame(method = c("CAPS","ACE-seq"),
                      mapping_rate = c(90.73, 71.78))

df_plot$method <- factor(df_plot$method, levels = c("CAPS","ACE-seq"))

ggplot(df_plot, aes(x = method, y = mapping_rate)) + geom_bar(stat = "identity", fill = "grey", col = "grey20") + 
  geom_text(aes(label=round(mapping_rate, 1)), vjust=-0.4) +
  theme_light() + xlab("") + ylab("Mapping rate (%)") + theme(axis.text = element_text(color = "black")) + ylim(0,100)

ggsave("plots/basic_statistics20200726/CAPS_mapping_rate.pdf",width = 2.3, height = 3.5)
