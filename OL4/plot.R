library(ggplot2)

# Create violin plot with Species on x-axis and Sepal.Length on y-axis
p <- ggplot(iris, aes(x = Species, y = Sepal.Length, fill = Species)) +
  geom_violin() +

  # Set y-axis range
  scale_y_continuous(limits = c(0.5, 7)) +

  # Set fill colors for each Species
  scale_fill_manual(values = c("#C44E52", "#55A868", "#4C72B0")) +

  # Set plot title and formatting
  labs(title = "Sepal Length Distribution", 
       x = "Species", 
       y = "Sepal Length") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# save the picture
ggsave("violin_plot.pdf", plot=p, height = 3, width = 5)