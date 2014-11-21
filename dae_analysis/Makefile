

## Mm9mESCTADsClassifiedByRT.csv is a dump of a file sent by B. Pope
## and described as (6/4/2013):
## "... a list of mm9 mESC TADs (columns for Chromosome, Start and End positions
## for each TAD). Even though they're all different sizes (ranging from 160kb to 5Mb)
## most are around 1Mb, so I think it would be best for each TAD to be an anchor instead
## of considering a +/-500kb region about each TAD center."
## then on 6/10/2013
## "... I've attached the list of TADs with another column ("RT") that classifies
## each TAD as either "E" (early), "L" (late) or "TTR" (contains a transition
## between early and late). About half of the TADs fall into the "TTR" category,
## most of which we'd expect to have properties similar to the "L" category,
## but in general these TADs are harder to classify as either "E" or "L". "

## peaks are downloaded from the mouse encode DCC



ANCHOR_RAD = 500000
PEAKS      = $(wildcard peaks/*.bed)
BIN        = 200
ARCHNAME   = W$(ANCHOR_RAD)_B$(BIN)
#DSNAME     = 
#DSPATH     = $(ARCHNAME).h5:$(DSNAME)
TRSIZE     = 200

#all : $(ARCHNAME).tr$(TRSIZE).h5 $(ARCHNAME).h5
all : $(ARCHNAME).h5

.PHONY : train

train : train_dae.cfg $(ARCHNAME).h5 
	dimer_train_dae.py --raw $(ARCHNAME).h5:train$(TRSIZE) $< $(shell ./hiddens.py --nin 7 --nout 5 $(ARCHNAME).h5 4)

$(ARCHNAME).h5 : anchors.tr$(TRSIZE).bed anchors.bed $(PEAKS) 
	@echo "training dataset"
	dimer_dataset_from_peaks.py --fit --bin_size $(BIN) --feature_score anchors.tr$(TRSIZE).bed $(PEAKS) $@:train$(TRSIZE)
	./label_data.py $@:train$(TRSIZE) --labels E TTR L
	@echo "full dataset"
	dimer_dataset_from_peaks.py --fit --bin_size $(BIN) --feature_score anchors.bed $(PEAKS) $@:full
	./label_data.py $@:full --labels E TTR L

#$(ARCHNAME).tr$(TRSIZE).h5 : anchors.tr$(TRSIZE).bed $(PEAKS)
#	dimer_dataset_from_peaks.py --fit --bin_size $(BIN) --feature_score $^ $@:$(DSNAME)
#	./label_data.py $@:$(DSNAME) --labels E TTR L
#
#$(ARCHNAME).h5 : anchors.bed $(PEAKS)
#	dimer_dataset_from_peaks.py --fit --bin_size $(BIN) --feature_score $^ $@:$(DSNAME)
#	./label_data.py $@:$(DSNAME) --labels E TTR L

#chrom start end class_NR class_coded_as_nr strand
anchors.bed : Mm9mESCTADsClassifiedByRT.csv
	cat $< | awk -F "," 'NR-1{mid_pt=($$2+$$3)/2; print $$1,mid_pt-$(ANCHOR_RAD),mid_pt+$(ANCHOR_RAD),$$4}' \
			| awk '{y=0; if($$4=="TTR") y=1; if($$4=="L") y=2; printf("%s\t%d\t%d\t%s\t%d\t+\n", $$1,$$2,$$3,$$4"_"NR,y)}' \
			| tr " " "\t" >$@

anchors.tr$(TRSIZE).bed : Mm9mESCTADsClassifiedByRT.csv
	cat $< | random_lines.py $(TRSIZE) \
			| awk -F "," 'NR-1{mid_pt=($$2+$$3)/2; print $$1,mid_pt-$(ANCHOR_RAD),mid_pt+$(ANCHOR_RAD),$$4}' \
			| awk '{y=0; if($$4=="TTR") y=1; if($$4=="L") y=2; printf("%s\t%d\t%d\t%s\t%d\t+\n", $$1,$$2,$$3,$$4"_"NR,y)}' \
			| tr " " "\t" >$@

#$(ARCHNAME)_$(DSNAME)_X.npy : dimer_dsdump.py $(ARCHNAME).h5
#	python $< $(DSPATH)


.PHONY : clean 

clean : 
	rm anchors.bed $(ARCHNAME).h5 anchors.tr$(TRSIZE).bed 
