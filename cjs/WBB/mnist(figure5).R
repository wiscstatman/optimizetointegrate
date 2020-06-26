# results procuded by "mnist.py"
acc = read.table("test_1.txt")
acc_df = data.frame(acc)
colnames(acc_df) = 'accuracy'
p <- ggplot(acc_df, aes(x=accuracy)) + 
  xlim(0,1) +
  geom_density(color='red',fill='red',alpha=0.2)+
  theme_bw() + 
  ##ggtitle("Weighted Bayesian Bootstrap")+
  theme(text = element_text(size=22),
        plot.title = element_text(hjust = 0.5))
pdf('acc.pdf',height = 6,width = 8)
print(p)
dev.off()
