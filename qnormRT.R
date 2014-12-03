qnormRT <-
function( RepliSeqBedGraphs , RepliChipFiles, cores="max" ){
	
	arrayfiles <- RepliChipFiles
	seqfiles <- RepliSeqBedGraphs
	
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	
	
	numseqfiles=length(seqfiles)
	numarrayfiles=length(arrayfiles)
	
	# find the number of samples in each file containing Repli-chip data
	arrayheads<-lapply(1:numarrayfiles,function(x) unlist(read.delim(pipe(paste("head -n 1",arrayfiles[x])), stringsAsFactors=F, header=F) ))
	arraysamplenames<-lapply(arrayheads,"[",-c(1,2))
	numarraycells<-unlist(lapply(arraysamplenames,length))
	
	# create a pool of all scores in all samples in all files
	allscores<-unlist(mclapply(1:(numarrayfiles+numseqfiles), function(x){
		
		if(x<=numseqfiles){
			cat("loading",seqfiles[x],"\n")
			as.numeric(readLines(pipe(paste("cut -f 4",seqfiles[x]))))
		} else{
			a=x-numseqfiles
			cat("loading",arrayfiles[a],"\n")
			d<-unlist(read.delim(pipe(paste("cut -f",paste(3:numarraycells[a],collapse=","),arrayfiles[a]))))
			d<-as.numeric(d[2:length(d)])
			return(d)
		}
		
	}, mc.cores=cores))
	
	outnames1<-paste0("repliseq_",basename(removeext(seqfiles)),"_qnormToPool.bg")
	outnames2<-lapply(1:numarrayfiles,function(x) paste0("replichip_",basename(removeext(arrayfiles[x])),arraysamplenames[[x]],"_qnormToPool.txt"))
	allnames<-c(outnames1,unlist(outnames2))

	# quantile normalize each Repli-seq file and save to a new bedGraph
	cat("normalizing sequencing data\n")
	if(cores>numseqfiles){cores2=numseqfiles} else{cores2=cores}
	mclapply(1:numseqfiles, function(x){
		seqdata<-read.tsv(seqfiles[x])
		curref<-sample(allscores,nrow(seqdata))
		seqdata[,4][order(seqdata[,4])]<-curref[order(curref)]
		write.tsv(seqdata,file=outnames1[x])
	}, mc.cores=cores2)
	
	# quantile normalize each Repli-chip file and save to a new file
	cat("normalizing array data\n")
	if(cores>numarrayfiles){cores2=numarrayfiles} else{cores2=cores}
	mclapply(1:numarrayfiles, function(x){
		for(i in 1:numarraycells[x]){
			arraydata<-read.delim(pipe(paste("cut -f 1,2,",i+2," ",arrayfiles[x],sep="")),header=TRUE)
			curref<-sample(allscores,nrow(arraydata))
			arraydata[,3][order(arraydata[,3])]<-curref[order(curref)]
			write.tsv(arraydata,file=outnames2[[x]][i])
		}
	}, mc.cores=cores2)

	return(allnames)
}
