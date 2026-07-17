### prewd ###

if (!any(ls() %in% "prewd")){prewd = file.path("/data2","xiongl2","Projects")}

### LeftRight function ###

source(file.path(prewd,"CalculationUnit","R","Functions","LeftRight.R"))

### libraries ###

library(foreach)
library(doParallel)
library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)

### mouse.chrs and intron ###
load("/data2/xiongl2/Projects/ProcessedData_Mouse/BasicObjects/mouse.chrs.lengths.RData")
load("/data2/xiongl2/Projects/ProcessedData_Mouse/BasicObjects/mouse.chrs.RData")
load("/data2/xiongl2/Projects/ProcessedData_Mouse/BasicObjects/mouse.gencode.basic.extended.RData")
mouse.gencode.basic.intron.anno.ranges = as(mouse.gencode.basic.extended[which(mouse.gencode.basic.extended[,"type"] == "intron"),],"GRanges")
	
### gather parameters ###

args = commandArgs(trailingOnly=TRUE)
sample = args[1]

### setwd ###

setwd(file.path(prewd,"ProcessedData_Mouse"))

### number of cores ###

mc.cores = 25 #detectCores()

## output ##
dir.create(file.path("BigWigTrack",sample))
  
load(file.path(prewd,"Rawdata",sample,"Bamfiles","size.factors.bam.RData"))
	
bam.files = list.files(pattern = "*.bam$",file.path(prewd,"Rawdata",sample,"Bamfiles"))
names(bam.files) = bam.files

### TT-Seq_Coverage.list ###

for (bam.file in bam.files){
    dir.create(file.path("BigWigTrack",sample,bam.file))	
	registerDoParallel(cores = mc.cores)
		build.Coverage.list = function(which.chr){
			Coverage.list.chr = list()
			param = ScanBamParam(which=GRanges(seqnames = which.chr,ranges = IRanges(0,mouse.chrs.lengths[which.chr])))
			bam = readGAlignmentPairs(file = file.path(prewd,"Rawdata",sample,"Bamfiles",bam.file),param = param,strandMode=2)
			bam = bam[start(left(bam)) <= end(right(bam))]
			bam = bam[(end(right(bam)) - start(left(bam))) <= 500]
			complement.bam = bam[(start(right(bam)) - end(left(bam))) < 2]
			bam = bam[(start(right(bam)) - end(left(bam))) >= 2]
			
			inner.mate.granges = GRanges(seqnames = which.chr,strand = strand(bam),ranges = IRanges(start = end(left(bam)) + 1,end = start(right(bam)) - 1))
            refseq.intron.anno.ranges = mouse.gencode.basic.intron.anno.ranges[seqnames(mouse.gencode.basic.intron.anno.ranges) == which.chr]
            inner.mate.intron.overlaps = findOverlaps(refseq.intron.anno.ranges,inner.mate.granges,maxgap = 0L,minoverlap = 1L,type = "within",select = "all",ignore.strand = FALSE)
			bam = bam[setdiff(1:length(bam),subjectHits(inner.mate.intron.overlaps))]
			
			rle.vec = Rle(0,mouse.chrs.lengths[which.chr])
			coverage.vec = coverage(bam[strand(bam) == "+"])[[which.chr]]
			rle.vec[1:length(coverage.vec)] = coverage.vec
			coverage.vec = coverage(GRanges(seqnames = which.chr,ranges = IRanges(start = end(left(bam[strand(bam) == "+"])) + 1,end = start(right(bam[strand(bam) == "+"])) - 1)))[[which.chr]]
			rle.vec[1:length(coverage.vec)] = rle.vec[1:length(coverage.vec)] + coverage.vec
			
			### complement.bam ###
			if (length(complement.bam[strand(complement.bam) == "+"]) > 0) {
			coverage.vec = coverage(GRanges(seqnames = which.chr,ranges = IRanges(start = start(left(complement.bam[strand(complement.bam) == "+"])),end = end(right(complement.bam[strand(complement.bam) == "+"])))))[[which.chr]]
			rle.vec[1:length(coverage.vec)] = rle.vec[1:length(coverage.vec)] + coverage.vec
			}
			
			Coverage.list.chr[["+"]] = rle.vec
		
			rle.vec = Rle(0,mouse.chrs.lengths[which.chr])
			coverage.vec = coverage(bam[strand(bam) == "-"])[[which.chr]]
			rle.vec[1:length(coverage.vec)] = coverage.vec
			coverage.vec = coverage(GRanges(seqnames = which.chr,ranges = IRanges(start = end(left(bam[strand(bam) == "-"])) + 1,end = start(right(bam[strand(bam) == "-"])) - 1)))[[which.chr]]
			rle.vec[1:length(coverage.vec)] = rle.vec[1:length(coverage.vec)] + coverage.vec
			
			### complement.bam ###
			if (length(complement.bam[strand(complement.bam) == "-"]) > 0) {
			coverage.vec = coverage(GRanges(seqnames = which.chr,ranges = IRanges(start = start(left(complement.bam[strand(complement.bam) == "-"])),end = end(right(complement.bam[strand(complement.bam) == "-"])))))[[which.chr]]
			rle.vec[1:length(coverage.vec)] = rle.vec[1:length(coverage.vec)] + coverage.vec
			}
			
			Coverage.list.chr[["-"]] = rle.vec
			
			save(Coverage.list.chr,file = file.path("BigWigTrack",sample,bam.file,paste0("Coverage.list.",which.chr,".RData")))
			return(Coverage.list.chr)
		}
		
		Coverage.list = foreach(n = mouse.chrs,.noexport = setdiff(ls(),c("mouse.chrs.lengths"))) %dopar% build.Coverage.list(n)
		save(Coverage.list,file = file.path("BigWigTrack",sample,bam.file,paste0("Coverage.list.",bam.file,".RData")))

        names(Coverage.list) = names(mouse.chrs.lengths)
        Coverage.list.plus = lapply(Coverage.list,function(x) x = x[[1]]/size.factors[bam.file])
        coverage.granges.plus = GRanges()
        for (i in mouse.chrs){
            Rle.vec = RleList(Coverage.list.plus[[i]])
            names(Rle.vec) = i
            coverage.granges.plus = c(coverage.granges.plus,as(Rle.vec,"GRanges"))
            }

        Coverage.list.minus = lapply(Coverage.list,function(x) x = x[[2]]*(-1)/size.factors[bam.file])
        coverage.granges.minus = GRanges()
        for (i in mouse.chrs){
             Rle.vec = RleList(Coverage.list.minus[[i]])
             names(Rle.vec) = i
             coverage.granges.minus = c(coverage.granges.minus,as(Rle.vec,"GRanges"))
             }
        
        Sample = strsplit(bam.file,split = "\\.")[[1]][1]
        export.bw(coverage.granges.plus,file.path("BigWigTrack",sample,bam.file,paste0(Sample,".plus.bw")))
        export.bw(coverage.granges.minus,file.path("BigWigTrack",sample,bam.file,paste0(Sample,".minus.bw")))
}