#
# Generate random data simulating genomic
# and RH markers and how they relate
#

# 26 chromosomes and an Un
chrs <- c(1:26, "Un")

# build the random data frame
chr <- vector()
pos <- vector()
rhchr <- vector()
rhpos <- vector()

for (c in chrs) {

	chr   <- c(chr, rep(c, 1000))
	pos   <- c(pos, round(seq(1,1000000, length=1000) + runif(1000, max=350)))

	rhchr <- c(rhchr, rep(c, 1000))
	rhpos <- c(rhpos, seq(1,15000, length=1000) + runif(1000, max=35))

}

d <- data.frame(chr=chr, pos=pos, rhchr=rhchr, rhpos=rhpos)



#
# Actually plot the comparison
#

# we need to find the maximum position
# for each chromosome in genomic and RH 
# co-ordinates and also the cumulative sum
maxp <- aggregate(d$pos, by=list(d$chr), max)
maxp$cs <- cumsum(maxp$x)
max1 <- aggregate(d$rhpos, by=list(d$rhchr), max)
max1$cs <- cumsum(max1$x)

# find the maxima
cs <- sum(maxp$x)
c1 <- sum(max1$x)

# empty plot capable of showing all
# chromosomes end to end
par(mar=c(3.5,3.5,1,1))
plot(1:10, 1:10, pch="", xlim=c(0,cs), ylim=c(0,c1), xlab="", ylab="", xaxt="n", yaxt="n", bty="n")

# grid lines
for (i in 1:(nrow(max1)-1)) {
	lines(c(1,cs), c(max1$cs[i],max1$cs[i]), lty=2, col="grey")
}
for (i in 1:(nrow(maxp)-1)) {
	lines(c(maxp$cs[i],maxp$cs[i]), c(1,c1), lty=2)
}

# colours
library(RColorBrewer)
cols <- colorRampPalette(brewer.pal(8, "Set1"))(length(chrs))

# index
idx <- 1

# plot each chromosome in turn
for (chr in chrs) {

	# get a data frame representing
	# just this chromosome
	sd <- d[d$chr==chr,]

	# order by position
	sd <- sd[order(sd$pos),]

	# set the starting co-ordinate on the x-axis
	# which is the maximum of the previous chromosome
	mp <- 0
	if (chr != 1) {
		mp <- maxp$cs[maxp[,1]==chrs[(1:length(chrs))[chrs==chr]-1]]
	}

	# literally plot each point
	for (i in 1:nrow(sd)) {

		# set the starting co-ordinate on the y-axis
		# which is the maximum of the previous chromosome
		ychr <- sd$rhchr[i]
		m1 <- 0
		if (ychr != 1) {
			m1 <- max1$cs[max1[,1]==chrs[(1:length(chrs))[chrs==chr]-1]]
		}

		# plot the points
		points(sd$pos[i]+mp, sd$rhpos[i]+m1, pch=16, col=cols[idx], cex=0.8) 
	
	}

	idx <- idx + 1

}

# plot X and Y axis
ypos <- 1:length(chrs)
ypos[1] <- max1$cs[1] / 2
for (i in 2:length(chrs)) {
	ypos[i] <- max1$cs[i-1] + (max1$cs[i] - max1$cs[i-1]) / 2
}

xpos <- 1:length(chrs)
xpos[1] <- maxp$cs[1] / 2
for (i in 2:length(chrs)) {
	xpos[i] <- maxp$cs[i-1] + (maxp$cs[i] - maxp$cs[i-1]) / 2
}

axis(side=2, at=ypos, labels=paste("chr",chrs, sep=""), las=2, tick=FALSE, line=-1.2)
axis(side=1, at=xpos, labels=paste("chr",chrs, sep=""), las=2, tick=FALSE, line=-1)

axis(side=1, at=cs/2, labels="Genomic pos", tick=FALSE, line=1)
axis(side=2, at=c1/2, labels="RH1 pos", tick=FALSE, line=1)
