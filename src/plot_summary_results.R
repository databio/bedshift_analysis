# Script to produce plots in R/ggplot from summary result table output files.

library("data.table")
library("ggplot2")
library("pepr")

# project_name = "pep_main"
# pMain = Project(paste0(project_name, "/project_config.yaml"))
# stMain = sampleTable(pMain)
# stMain$universe

project_name = "pep_universe"
project_name = "pep_main"

p = Project(paste0(project_name, "/project_config.yaml"))
result_folder =paste0("results/", project_name)
st = sampleTable(p)
results = list(
	"jaccard" = paste0(result_folder, "/scores/jaccard_results.txt"),
	"coverage" = paste0(result_folder, "/scores/coverage_results.txt"),
	"cosine" = paste0(result_folder, "/scores/cosine_results.txt"),
	"euclidean" = paste0(result_folder, "/scores/euclidean_results.txt"))

# file="pep_demo/results/scores/jaccard_results.txt"
# name="jaccard"
transformResults = function(file, name="") {
	rawdata = fread(file)
	cm = colMeans(rawdata)
	tcm = as.data.frame(cbind(cm))
	colnames(tcm) = "value"
	sampleNames = rownames(tcm)
	tcm$sample_name = sampleNames
	tcm$metric = name
	return(tcm)
}

st[is.na(op2), scenario1:=op1]
st[is.na(op2), scenarioClass1:=toupper(op1)]
st[!is.na(op2), scenario1:=paste0(op1, op1deg, "+",toupper(op2))]
st[!is.na(op2), scenario2:=paste0(toupper(op1), "+", op2, op2deg)]
st[!is.na(op2), scenarioClass1:=paste0(op1, "+", toupper(op2))]
st[!is.na(op2), scenarioClass2:=paste0(toupper(op1), "+", op2)]
st$scenario1
st$scenario2
st$scenarioClass1
st$scenarioClass2

readin = lapply(results, transformResults)
rbl = rbindlist(readin)
rbl$metric = rep(names(readin), sapply(readin, NROW))
rbl
rbl2 = copy(rbl)

rbl$scenario=rep(st$scenario1,4)
rbl$scenarioClass=rep(st$scenarioClass1,4)
rbl$rep=rep(st$rep1,4)
rbl$univ=rep(gsub(".*GRCh38-ccREs.?(.*).bed", "\\1", unlist(st$universe)),4)
rbl$file=rep(substr(basename(unlist(st$file)), 1,5),4)


rbl2$scenario=rep(st$scenario2,4)
rbl2$scenarioClass=rep(st$scenarioClass2,4)
rbl2$rep=rep(st$rep2,4)
rbl2$univ=rep(gsub(".*GRCh38-ccREs.?(.*).bed", "\\1", unlist(st$universe)),4)
rbl2$file=rep(substr(basename(unlist(st$file)), 1,5),4)

rblc = rbind(rbl, rbl2[!is.na(scenario),])
rblc


# Results table
dcast(rblc, sample_name + file ~ metric, value.var="value")[file=="713f5",]
dcast(rblc, sample_name + file ~ metric, value.var="value")[file=="713f5",]



ymin = floor(min(rblc$value)*10)/10

p = ggplot(rblc, aes(y=value, x=rep, color=metric)) + geom_line(stat="identity", size=1.25, alpha=0.75) + 
	theme_classic() + ylim(ymin, 1) + 
	theme(panel.grid.major.y=element_line("#CCCCCC", size=0.2)) +
	theme(strip.text.x = element_text(size = 6)) +
	scale_x_continuous(breaks=1:3, labels=c("L", "M", "H")) + 
	geom_point(alpha=0.75) + theme(strip.background = element_blank()) + xlab(element_blank()) + ylab("Similarity score")

if (length(unique(rbl$univ)) >1) {
	p = p + facet_wrap(scenario + univ ~ ., scales="free_x", nrow=3)
} else if (length(unique(rbl$file) > 1)) {
	p = p + facet_wrap(scenario + file ~ ., scales="free_x", nrow=3)
} else {
	p = p + facet_wrap(scenario ~ ., scales="free_x", nrow=3)
}

pdf(paste0("results/", project_name, "_analysis_all.pdf"), width=17, height=5)
print(p)
dev.off()

# average plots
proc = function(x) {
	100*(x[1]-x[3])/x[1]
}

rblc[scenarioClass=="Shift+DROP" & file=="713f5"]
rblc[, c(rep), by=list(scenario, file, metric)]$V1

scores = rblc[,list(value= proc(value), scenarioClass), by=list(scenario, file, metric)]
scoresave = scores[, list(value=mean(value)), by=list(scenarioClass, metric, file)]
scoresub = scoresave[file=="713f5",]

pslope = ggplot(scoresub, aes(y=value, x=scenarioClass, fill=metric)) + geom_bar(stat="identity", position="dodge", size=1.25, alpha=0.75) + 
	facet_grid(. ~ scenarioClass, scales="free_x") + theme_classic() +
	theme(panel.grid.major.y = element_line(colour = "gray"),
	panel.grid.minor.y = element_line(colour = "gray"),
	panel.ontop=FALSE) + ylab ("Similarity change (%)") + theme(strip.background = element_blank())

pdf(paste0("results/", project_name, "_analysis_slope_summary.pdf"), width=7, height=2.5)
print(pslope)
dev.off()



# p = ggplot(scores, aes(y=V1, x=file, fill=metric)) + geom_bar(stat="identity", position="dodge", size=1.25, alpha=0.75) + theme_classic() + 
# 	facet_grid(. ~ set, scales="free_x")

# ggplot(scores, aes(y=V1, x=set, color=metric)) + geom_jitter(stat="identity", size=1.25, alpha=0.75) + theme_classic()

# ggplot(scores, aes(y=V1, x=file, color=metric)) + geom_jitter(stat="identity", size=1.25, alpha=0.75) + theme_classic()

#  + 
# 	facet_grid(. ~ set, scales="free_x")
#  + 
# 	facet_grid(. ~ set + univ + file, scales="free_x") + ylim(0.2, 1) + 
# 	theme(panel.grid.major.y=element_line("#CCCCCC", size=0.2)) +
# 	theme(strip.text.x = element_text(size = 6)) +
# 	scale_x_continuous(breaks=1:3, labels=c("L", "M", "H")) + 
# 	geom_point(alpha=0.75) + theme(strip.background = element_blank()) + xlab(element_blank())
