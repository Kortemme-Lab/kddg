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
head(a)
coefs <- coef(lm(PredictedDDG~%(experiment_field)s, data = a))
# Sanity check
coefs
fitcoefs = coef(lm(PredictedDDG~0 + %(experiment_field)s, data = a))
fitlmv_PredictedDDG <- as.numeric(fitcoefs[1])

# coefs contains two values: (Intercept) and %(experiment_field)s
lmv_intercept <- as.numeric(coefs[1])
lmv_PredictedDDG <- as.numeric(coefs[2])

lm(a$PredictedDDG~a$%(experiment_field)s)
fitcoefs

xlabel <- expression(paste(plain("Experimental ")*Delta*Delta*plain("G (kcal/mol)")))
ylabel <- expression(paste(plain("Predicted ")*Delta*Delta*plain(G)))
rvalue <- cor(a$PredictedDDG, a$%(experiment_field)s)
rvalue

# To change the font size of the axis labels (tick labels), use e.g.:
# 	p <- p + theme(axis.text.x=element_text(size=22))
# To change the font of the axis titles, use e.g.:
# 	p <- p + theme(axis.title.x = element_text(face="bold", colour="#990000", size=20),

# shape I(20) is a small dot, I(19) is a large dot, I(4) is a cross

p <- qplot(%(experiment_field)s, PredictedDDG, main="%(title)s", data=a, xlab=xlabel, ylab=ylabel, shape = I(19), alpha = I(txtalpha)) + # label=ProThermID
		geom_abline(size = 0.25, intercept = lmv_intercept, slope = lmv_PredictedDDG) +
		geom_abline(color="blue",size = 0.25, intercept = 0, slope = fitlmv_PredictedDDG  )

if ('%(filetype)s' == 'pdf'){
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

aexp = a$%(experiment_field)s
apre = a$PredictedDDG
maevalue <- MAE(aexp, apre)
p <- p + geom_text(hjust=0, size=8, aes(xpos, ypos_mae, fontface="plain", family = fface, label=sprintf("MAE = %%0.4f", round(maevalue, digits = 4))))

#p <- p + geom_text(size=6, fontfamily = 'sans family', aes(min(a$%(experiment_field)s) + 3+1,max(a$PredictedDDG)+1, label=sprintf("cor(y,x) = %%f", round(rvalue, digits = 4)))) 


# Plot graph
p

#text(x, y, , cex=2.0)

#p + geom_text(size=4, fontfamily = 'courier', aes(min(a$%(experiment_field)s) + 3+1,max(a$PredictedDDG)+1, label=sprintf("cor(y,x) = %%f", round(rvalue, digits = 4)))) + 
#		geom_text(size=2, colour = 'red', fontfamily = 'courier', data=a[68, ], label=a$ProThermID[68], vjust=1) +
#		geom_text(size=2, colour = 'red', fontfamily = 'courier', data=a[76, ], label=a$ProThermID[76], vjust=1) + 
#		geom_text(size=2, colour = 'red', fontfamily = 'courier', data=a[307, ], label=a$ProThermID[307], vjust=1) + 
#		geom_text(size=2, colour = 'red', fontfamily = 'courier', data=a[650, ], label=a$ProThermID[650], vjust=1) + 
#		geom_text(size=2, colour = 'red', fontfamily = 'courier', data=a[727, ], label=a$ProThermID[727], vjust=1) + 
#		geom_text(size=2, colour = 'red', fontfamily = 'courier', data=a[842, ], label=a$ProThermID[842], vjust=1) + 
#		geom_text(size=2, colour = 'red', fontfamily = 'courier', data=a[918, ], label=a$ProThermID[918], vjust=1) + 
#		geom_text(size=2, colour = 'red', fontfamily = 'courier', data=a[268, ], label=a$ProThermID[268], vjust=1)  

dev.off()