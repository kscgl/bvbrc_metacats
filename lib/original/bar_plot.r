exponent <- seq(from=0, to=-10, by=-1)
yaxis <- 10^exponent
yrange <- c(1, 1e-10)
par(bg="lightgray")
barplot2(sigpvals,ylim=yrange,yaxt="n",log="y")
axis(2,at=yaxis,labels=yaxis,cex.axis=0.8,las=2)
abline(h=1)
abline(h=pvalcutoff,col="red")
cutofflabel <- paste("cutoff:",pvalcutoff)
mtext(cutofflabel, at=pvalcutoff, side=4, col="red")
barplot2(sigpvals,main="Chi-square Test P-values Bar Plot",xlab="Positions",ylab="P-Values",names.arg=positions,col="#2B60DE",ylim=yrange,yaxt="n",log="y",font.lab=2, plot.grid=TRUE, add=TRUE)
dev.off()
}