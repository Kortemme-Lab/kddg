library(ggplot2)
library(gridExtra) 
library(scales)
library(qualV)

if ('%(filetype)s' == 'pdf'){
	%(filetype)s('%(output_filename)s', paper="special", width=12, height=12) # otherwise postscript defaults to A4, rotated images
} else if ('%(filetype)s' == 'png'){
	%(filetype)s('%(output_filename)s', height=4096, width=4096, bg="white", res=600)
} else if ('%(filetype)s' == 'postscript'){
	%(filetype)s('%(output_filename)s', horizontal=FALSE, paper="special", width=12, height=12) # otherwise postscript defaults to A4, rotated images
}

#Use alpha 0.1 for PDF
#Use alpha 0.25 for PNG
if ('%(filetype)s' == 'pdf')
{
	txtalpha <- 0.3
	redtxtalpha <- 0.3
} else if ('%(filetype)s' == 'postscript') { # postscript does not handle transparency
	txtalpha <- 1.0
	redtxtalpha <- 1.0
} else if ('%(filetype)s' == 'png'){
	txtalpha <- 0.25
	redtxtalpha <- 0.5
}

par(mar=c(5, 5, 1, 1))
a <- read.csv('%(inputfname)s', header=T)
res <- resid(mod <- lm(PredictedDDG~%(experiment_field)s, data = a))
res.qt <- quantile(res, probs = c(0.05,0.95))
want <- which(res >= res.qt[1] & res <= res.qt[2])

coefs <- coef(lm(PredictedDDG~%(experiment_field)s, data = a[want,]))
fitcoefs = coef(lm(PredictedDDG~0 + %(experiment_field)s, data = a[want,]))
fitlmv_predicted <- as.numeric(fitcoefs[1])
rvalue <- cor(a[want,]$Experimental, a[want,]$Predicted)

paste('PYTHON_VALUE', 'float', 'correlation', rvalue)

lmv_intercept <- as.numeric(coefs[1])
lmv_PredictedDDG <- as.numeric(coefs[2])

xlabel <- expression(paste(plain("Experimental ")*Delta*Delta*plain("G (kcal/mol)")))
ylabel <- expression(paste(plain("Predicted ")*Delta*Delta*plain(G)))

# To change the font size of the axis labels (tick labels), use e.g.:
# 	p <- p + theme(axis.text.x=element_text(size=22))
# To change the font of the axis titles, use e.g.:
# 	p <- p + theme(axis.title.x = element_text(face="bold", colour="#990000", size=20),

# shape I(20) is a small dot, I(19) is a large dot, I(4) is a cross

p <- qplot(%(experiment_field)s, PredictedDDG, main="%(title)s", data=a[want,], xlab=xlabel, ylab=ylabel, shape = I(19), alpha = I(txtalpha)) + # label=ProThermID
		geom_abline(size = 0.25, intercept = lmv_intercept, slope = lmv_PredictedDDG) +
		geom_abline(color="blue",size = 0.25, intercept = 0, slope = fitlmv_PredictedDDG) +
		geom_point(data=a[-want,], color = "red", shape = I(19), alpha = I(redtxtalpha))

if ('%(filetype)s' == 'pdf'){
 	p <- p + theme(plot.title = element_text(size=25))
 	p <- p + theme(axis.title.x = element_text(size=45, vjust=-1.5)) # vjust for spacing
	p <- p + theme(axis.title.y = element_text(size=45))
	p <- p + theme(axis.text.x=element_text(size=25))
	p <- p + theme(axis.text.y=element_text(size=25))
}

#courier
#Times New Roman
#sans family

# Create labels for cor(y,x) and MAE 
# Using hjust=0 in geom_text sets text to be left-aligned
 
minx <- min(a$%(experiment_field)s)
maxx <- max(a$%(experiment_field)s)
miny <- min(a$PredictedDDG)
maxy <- max(a$PredictedDDG)

# fontface can be plain, bold, italic
# sans, serif, mono does not work for my PostScript driver
# palatino, bookman, helvetica, times works for PostScript but both look the same

if ('%(filetype)s' == 'postscript')
{
	fface <- "bookman"
} else {
	fface <- "sans"
}

xpos <- minx + ((maxx - minx) * 0.05)
ypos_cor <- maxy - ((maxy - miny) * 0.015)
ypos_mae <- maxy - ((maxy - miny) * 0.085) # different variable name as these seem to be evaluated later (if we use the same label, if affects the printing of cor(y,x) as well)
p <- p + geom_text(hjust=0, size=8, aes(xpos, ypos_cor, fontface="plain", family = fface, label=sprintf("cor(y,x) = %%f", round(rvalue, digits = 4)))) 

aexp = a[want,]$%(experiment_field)s
apre = a[want,]$PredictedDDG
maevalue <- MAE(aexp, apre)

paste('PYTHON_VALUE', 'float', 'MAE', maevalue)

p <- p + geom_text(hjust=0, size=8, aes(xpos, ypos_mae, fontface="plain", family = fface, label=sprintf("MAE = %%0.4f", round(maevalue, digits = 4))))


# Plot graph
p

dev.off()