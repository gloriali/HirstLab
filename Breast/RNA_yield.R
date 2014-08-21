# version 1 
summary<-data.frame(tissue=c("lum","myo"),
                    yield=c(mean(c(4.91, 1.83, 4.25)), mean(c(2.14, 0.31, 1.16))), 
                    col=c(rgb(200,50,0, maxColorValue = 255),rgb(50,200,50, maxColorValue = 255)))
error.bar <- function(x, y, upper, lower, length=0.1,...){
  arrows(x,upper,x,lower,angle=90,code=3,length=length,lwd=1.8,...)
}
ee_up<-c(max(c(4.91, 1.83, 4.25)), max(c(2.14, 0.31, 1.16)))
ee_low<-c(min(c(4.91, 1.83, 4.25)), min(c(2.14, 0.31, 1.16)))



pdf("~/快盘/REMC/figures/RNA_yield.pdf", height = 5)
par(mar = c(5,5,2,8), xpd = TRUE, cex.axis = 1.2, cex.lab = 1.5)
barx <- barplot(summary$yield, col = as.character(summary$col), ylim = c(0, 5), width = 0.5, cex.names = 1.2, names.arg = summary$tissue, axis.lty = 1, xlab = "Cell type", ylab="RNA yield")  
error.bar(barx, summary$yield, ee_up, ee_low)
legend("topright", as.character(summary$tissue), inset=c(-0.2, 0), col=as.character(summary$col), lty = 1, lwd = 8)
dev.off()

# version 2
y <- matrix(c(4.91, 1.83, 4.25, 2.14, 0.31, 1.16), ncol = 3, byrow = T)
cols <- rep(c(rgb(200,50,0, maxColorValue = 255),rgb(50,200,50, maxColorValue = 255)), each = 3)
pdf("~/快盘/REMC/figures/RNA_yield.pdf", height = 5)
par(mar = c(5,5,2,8), xpd = TRUE, cex.axis = 1.2, cex.lab = 1.5)
barplot(y, beside=TRUE, col=c(rgb(200,50,0, maxColorValue = 255),rgb(50,200,50, maxColorValue = 255)), ylim=c(0,5),cex.names=1.2,names.arg=c("RM084", "RM080", "RM035"), axis.lty=1, ylab="RNA yield")  
legend("topright", as.character(summary$tissue), inset=c(-0.2, 0), col=as.character(summary$col), lty = 1, lwd = 8)
dev.off()
