qnormall<-function( seqfiles , arrayfiles, cores="max" ){
library(parallel)
if(cores=="max"){cores=detectCores()-1}
numseqfiles=length(seqfiles)
numarrayfiles=length(arrayfiles)
arrayheads<-lapply(1:numarrayfiles,function(x) read.delim(pipe(paste("head",arrayfiles[x]))) )
numarraycells<-unlist(lapply(arrayheads,ncol))-2
cat("pooling all data\n")
#seqscores<-unlist(mclapply(1:numseqfiles, function(x) seqdata[[x]][,4] , mc.cores=cores))
#arrayscores<-unlist(mclapply(1:numarrayfiles, function(x) as.vector(arraydata[[x]][,3:ncol(arraydata[[x]])]), mc.cores=cores))
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
if(cores>numseqfiles){cores2=numseqfiles} else{cores2=cores}
cat("normalizing sequencing data\n")
mclapply(1:numseqfiles, function(x){
seqdata<-read.tsv(seqfiles[x])
curref<-sample(allscores,nrow(seqdata))
seqdata[,4][order(seqdata[,4])]<-curref[order(curref)]
outname<-paste(basename(removeext(seqfiles[x])),"_qnormToPool.bg",sep="")
print(outname)
write.tsv(seqdata,file=outname)
}, mc.cores=cores2)
if(cores>numarrayfiles){cores2=numarrayfiles} else{cores2=cores}
cat("normalizing array data\n")
mclapply(1:numarrayfiles, function(x){
for(i in 1:numarraycells[x]){
arraydata<-read.delim(pipe(paste("cut -f 1,2,",i+2," ",arrayfiles[x],sep="")),header=TRUE)
curref<-sample(allscores,nrow(arraydata))
arraydata[,3][order(arraydata[,3])]<-curref[order(curref)]
outname<-paste(basename(removeext(arrayfiles[x])),"_",colnames(arraydata)[3],"_qnormToPool.txt",sep="")
print(outname)
write.tsv(arraydata,file=outname)
}
}, mc.cores=cores2)
}
