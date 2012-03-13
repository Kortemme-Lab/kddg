import sys
import os
import re
import pickle
import subprocess
from string import join
import time
import common.colortext as colortext
import common.ddgproject as ddgproject
from common.rosettahelper import kJtokcal
from tempfile import mkstemp

ddGdb = ddgproject.ddGDatabase()
dbfields = ddgproject.FieldNames()

def _mean(points, numpoints):
	'''Points is expected to be a list (or iterable), numpoints should be an integer.'''
	return sum(points) / numpoints 

def _createMAEFile(results, outfname, average_fn = _mean):
	output = []
	output.append("X,Error")
	X = 1
	for r in results:
		predicted_ddG = pickle.loads(r["ddG"])["data"]["ddG"]
		
		expID = r["ExperimentID"]
		expScores = ddGdb.callproc('GetScores', parameters = (expID,))
		
		n = 0
		scores = []
		for e in expScores:
			c = e["NumberOfMeasurements"]
			n += c
			scores.append(e["ddG"] * c)
		 # Note the sign negation of average_fn(scores, n) as Rosetta convention is reverse to ProTherm
		point = abs(- kJtokcal(average_fn(scores, n)) - predicted_ddG)
		output.append("%s,%s"% (X,point))
		X += 1
		
	F, fname = mkstemp(dir = ".")
	F = os.fdopen(F, "w")
	F.write(join(output, "\n"))
	F.write("\n")
	F.close()
	return fname

def _createAveragedInputFile(results, outfname, average_fn = _mean):
	output = []
	output.append("Experimental,Predicted,ProThermID")
	for r in results:
		predicted_ddG = pickle.loads(r["ddG"])["data"]["ddG"] 
		expID = r["ExperimentID"]
		expScores = ddGdb.callproc('GetScores', parameters = (expID,))
		
		n = 0
		scores = []
		sources = []
		for e in expScores:
			sources.append(e["SourceID"])
			c = e["NumberOfMeasurements"]
			n += c
			scores.append(e["ddG"] * c)
		 # Note the sign negation as Rosetta convention is reverse to ProTherm.
		 # The scores stored in the database are in kJ/mol. Convert to kcal/mol
		eavg = - kJtokcal(average_fn(scores, n))
		output.append("%s,%s,%s"% (eavg, predicted_ddG, join(sorted(sources)," & ") or "."))

	F, fname = mkstemp(dir = ".")
	F = os.fdopen(F, "w")
	F.write(join(output, "\n"))
	F.write("\n")
	F.close()
	return fname

def _runRScript(RScript):
	F, rscriptname = mkstemp(dir = ".")
	F = os.fdopen(F, "w")
	F.write(RScript)
	F.close()
	p = subprocess.Popen(["R","CMD", "BATCH", rscriptname]) 
	while True:
		time.sleep(0.3)
		errcode = p.poll()
		if errcode != None:
			break
	rout = "%s.Rout" % rscriptname
	os.remove(rscriptname)
	if errcode != 0:
		if os.path.exists(rout):
			F = open(rout, "r")
			colortext.warning(F.read())
			F.close()
			os.remove(rout)
		raise colortext.Exception("The R script failed with error code %d." % errcode)
	os.remove(rout)
	
def _R_mean_unsigned_error(inputfname, outfname):
	RScript='''
pdf('%(outfname)s')
par(mar=c(5, 5, 1, 1))
a <- read.csv('%(inputfname)s', header=T)
head(a)
plot(a$Error, pch=19,xlab='Experiment #',ylab='MAE',cex.lab=1.5, xlim=c(1,max(a$X)+1), ylim=c(0,max(a$Error)+1), xaxt='n')
axis(side=1, at=0:23, cex.axis=.5)
?axis
dev.off()''' % vars()
	# todo: The graph size is not being set correctly i.e. the max function call is wrong above
	_runRScript(RScript)


def _R_correlation_coefficient(inputfname, outfname):
	RScript='''
pdf('%(outfname)s')
par(mar=c(5, 5, 1, 1))
a <- read.csv('%(inputfname)s', header=T)

aexp = a$Experimental   # the eruption durations 
apre = a$Predicted      # the waiting period 
rvalue <- cor(aexp, apre)          # apply the cor function
cor(aexp, apre)

head(a)
reg1 <- lm(a$Experimental~a$Predicted)
plot(
	a$Experimental, 
	a$Predicted,
	pch=19,
	xlab=expression(paste(plain("Experimental ")*Delta*Delta*plain(G))),
	ylab=expression(paste(plain("Computational ")*Delta*Delta*plain(G))),
	xlim=c(min(a$Experimental)-1,max(a$Experimental)+1),
	ylim=c(min(a$Predicted)-1,max(a$Predicted)+1),
	cex.lab=1.5,
	)
text(min(a$Experimental) + 2, max(a$Predicted), sprintf("R = %%f", round(cor(aexp, apre), digits = 4)))
abline(reg1)
dev.off()''' % vars()
	# todo: The graph size is not being set correctly i.e. the max function call is wrong above
	_runRScript(RScript)


def plot(RFunction, filecreator, results, outfname, average_fn = _mean):
	'''Results is expect to be a list of dicts each of which has the keys ExperimentID and ddG.'''
	if not results:
		raise Exception("The results set is empty.")
	elif len(results) == 1:
		raise Exception("The results set only has one data point. At least two points are required.")
	else:
		inputfname = filecreator(results, outfname, average_fn)
		try:
			colortext.printf("Running %s" % RFunction)
			print(inputfname)
			RFunction(inputfname, outfname)
		except Exception, e:
			import traceback
			colortext.error(traceback.format_exc())
			os.remove(inputfname)
			raise Exception(e)
		os.remove(inputfname)
			

# R notes for Shane
# par is used to specify graphical parameters. mar sets the margin sizes - c(bottom, left, top right)
# read.table is useful for reading delimited files. e.g. read.table("c:/mydata.csv", header=TRUE, sep=",", row.names="id")
# plot creates a plot
#	pch=19 means use a diamond character (web-search for 'R plot symbols' pch)
#	a$V1, a$V2 choose the first and second columns of data table a
#	xlab and ylab set the labels
#	cex	indicates the amount by which plotting text and symbols should be scaled relative to the default. 1=default, 1.5 is 50% larger, 0.5 is 50% smaller, etc.
#	cex.lab	defines the magnification of x and y labels relative to cex
#	abline(k,l) adds a straight line where k and l are the intercept and slope values respectively

