## Input set up
inputDir=      # directory of the reads of input samples
sampleDir=     # directory of the reads of ChIP samples
chromsizeDir=  # path to the chromsize folder under the directory where chromHMM is installed
hgcelltable=   # path to the human cell market table
mmcelltable=   # path to the mouse cell market table
rawbinarDir=   # directory of the binarized data
hgnames=       # names of the human samples, should be the same as the names in the cell market table
mmnames=       # names of the mouse samples, should be the same as the names in the cell market table
segDir=        # directory of the chromHMM segmentation outputs
numStates=     # number of states to learn in the model
targetSamples= # the names of the samples that will be used to map states on targets
targets=       # the RT boundaries to map states on, with this header: "chr position strand"
binSize=       # the size of bins in the extension of boundaries
binCount=      # number of bins to extend the boundaries on both sides
########################################

java -jar ChromHMM.jar BinarizeBed -c ${inputDir} ${chromsizeDir}/hg19.txt ${sampleDir} ${hgcelltable} ${binarDir}
java -jar ChromHMM.jar BinarizeBed -c ${inputDir} ${chromsizeDir}/CHROMSIZES/mm9.txt ${sampleDir} ${mmcelltable} ${binarDir}
cd ${rawbinarDir}
binarDir=${rawbinarDir}/modified
mkdir ${binarDir}
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX
do
for ct in ${hgnames}
do
sed "s/chr/hg19chr/g" ${ct}_${chr}_binary.txt >${binarDir}/${ct}_hg19${chr}_binary.txt
done
done

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX
do
for ct in ${mmnames}
do
sed "s/chr/mm9chr/g" ${ct}_${chr}_binary.txt >${binarDir}/${ct}_mm9${chr}_binary.txt
done
done

cd ${chromsizeDir}
sed 's/chr/mm9chr/g' mm9.txt >artificial_mm9.txt
sed 's/chr/hg19chr/g' hg19.txt >artificial_hg19.txt
cat artificial_mm9.txt artificial_hg19.txt >artificial_mm9andhg19.txt
rm artificial_mm9.txt artificial_hg19.txt

mkdir ${segDir}
cd ${segDir}
java -jar ChromHMM.jar LearnModel -l ${chromsizeDir}/artificial_mm9andhg19.txt ${binarDir} ${segDir} ${numStates} hg19

mkdir ./segments
cp *_segments.bed segments/
cd ./segments
for sg in `ls *.bed`
do
sed 's/mm9chr/chr/g' ${sg}|sed 's/hg19chr/chr/g' >returned_${sg}
mv returned_${sg} ${sg}
done

targetbins=${targets}_binned
sed '1d' ${targets} |awk -v sz=${binSize} -v sc=${binCount} '{for (i=-sc;i<=sc;i++) print $1,$2+i*sz,$2+i*sz}' |awk '{OFS="\t";print $0,NR}' >${targetbins}

cd ${segDir}

echo '
args=commandArgs(TRUE)
g=read.table(args[1],header=T)
BinStates=read.table(args[2])
BinStates_arranged_name=paste(as.character(args[2]),"rearranged",sep="_")

arngFun=function(xx)
{
yy=aggregate(xx,by=list(xx[,1]),min)
binNum=nrow(yy)/nrow(g)
zz=yy[sort.list(yy[,1]),]
tt=t(matrix(as.character(zz[,3]),binNum,nrow(g)))
for(i in 1:nrow(g))
{
if(g$strand[i]=="-")
{dd=tt[i,]
for(p in 1:binNum)
{
tt[i,p]=dd[(binNum+1-p)]
}
}
}
tt
}

BinStates_arranged=arngFun(BinStates)
write.table(BinStates_arranged,BinStates_arranged_name,quote=F,sep="\t",row.names=F,col.names=F)
' >rearrange.r

for smp in ${targetSamples}
do
intersectBed -a ${targetbins} -b ${smp}_*_segments.bed -wa -wb >${smp}_hit
intersectBed -v -a ${targetbins} -b ${smp}_*_segments.bed -wa -wb >${smp}_miss
awk <${smp}_miss '{OFS="\t";print $1,$2,$3,$4,"chr0","0","0","E0"}' >${smp}_miss_extended
cat ${smp}_hit ${smp}_miss_extended |cut -f 4,8 |sed "s/E/ /g" >${smp}_states_onbinpoints
Rscript rearrange.r ${targets} ${smp}_states_onbinpoints
rm ${smp}_hit ${smp}_miss ${smp}_miss_extended
done
rm rearrange.r



