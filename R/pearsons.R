require(qualV)

if ('%(filetype)s' == 'postscript') {
%(filetype)s('%(output_filename)s', horizontal=FALSE, paper="special", width=12, height=12) # otherwise postscript defaults to A4, rotated images
} else {
%(filetype)s('%(output_filename)s')
}

par(mar=c(5, 5, 1, 1))
a <- read.csv('%(inputfname)s', header=T)

aexp = a$%(experiment_field)s
apre = a$PredictedDDG
rvalue <- cor(aexp, apre)          # apply the cor function
cor(aexp, apre)

maevalue <- MAE(aexp, apre)

head(a)
reg1 <- lm(a$%(experiment_field)s~a$PredictedDDG)
plot(
	a$%(experiment_field)s, 
	a$PredictedDDG,
	main="%(title)s",
	pch=19,
	xlab=expression(paste(plain("Experimental ")*Delta*Delta*plain(G))),
	ylab=expression(paste(plain("Computational ")*Delta*Delta*plain(G))),
	xlim=c(min(a$%(experiment_field)s)-1,max(a$%(experiment_field)s)+1),
	ylim=c(min(a$PredictedDDG)-1,max(a$PredictedDDG)+1),
	cex.lab=1.5,
	)

abline(reg1)

pu <- par()$usr
#x <- pu[2] * 0.083 + pu[1] * 0.917 # when cex == 1
#y <- pu[3] * 0.05 + pu[4] * 0.95 # when cex == 1)
x <- pu[2] * 0.1 + pu[1] * 0.9 # when cex == 2
y <- pu[3] * 0.05 + pu[4] * 0.95 # when cex == 2 
text(x, y, sprintf("R = %%0.4f", round(cor(aexp, apre), digits = 4)), cex=2.0)

#x <- pu[2] * 0.1 + pu[1] * 0.9 # when cex == 1
#y <- pu[3] * 0.08 + pu[4] * 0.92 # when cex == 1
x <- pu[2] * 0.12 + pu[1] * 0.88 # when cex == 2
y <- pu[3] * 0.1 + pu[4] * 0.9 # when cex == 2
 
text(x, y, sprintf("MAE = %%0.4f", round(maevalue, digits = 4)), cex=2.0)
#text(min(a$%(experiment_field)s) + 2, max(a$PredictedDDG), sprintf("R = %%4f", round(cor(aexp, apre), digits = 4)))
#text(min(a$%(experiment_field)s) + 2, max(a$PredictedDDG) - 2, sprintf("MAEv = %%0.4f", round(maevalue, digits = 4)))

dev.off()