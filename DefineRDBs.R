DefineRDBs = function(x, Slope_E = 2.75e-6, Slope_L = 1e-06, Ends = 1e6, Data_Points = 30, Length_Max = 1e6, Length_Min = 2e5, RTU_Min = 0.55, Span = 35, Gap = 8e4, Gap_Dis = 125e3){
	cat("Initializing", "\n")
	require(gtools)
	chrs = levels(x$CHR)		# Define chromosomes names
	
	endPoints = data.frame(CHR=chrs, Start=-1, End=-1); P5 = NULL	 # Create empty variables to collect limts of mapped regions (endPoints) and gaps between data points (P5)
	j = 0
	for (chr in chrs) {
		j = j+1
		RTc = x$POSITION[x$CHR == chr]
		endPoints[j,] = cbind(CHR=chr, Start=(min(RTc)+Ends), End=(max(RTc)-Ends))
		
		L = NULL; L = lapply(2:(length(RTc)-1), function(x) RTc[x]-RTc[x-1]); L <- do.call(rbind, L)
		R = NULL; R = lapply(2:(length(RTc)-1), function(x) RTc[x+1]-RTc[x]); R <- do.call(rbind, R)
		
		P5 = rbind(P5, data.frame(CHR=chr,POSITION=RTc[2:(length(RTc)-1)],DISL=L,DISR=R))	
	}
	
	TTRr = NULL
	for(ct in 3:ncol(x)){
		cat("Mapping TTRs for dataset:", ct-2, "of", ncol(x)-2, "\n")	# Output progress
		highTTR = NULL;	lowTTR = NULL
		for(chr in chrs){
			cat("Current chromosome: ", chr, "\n")
			RTb		   = x[which(x$CHR == chr),c(2,ct)]								# Subset of data table (x) on current chromosome (chr) for current dataset (ct)
			loessProf  = loess(RTb[,2] ~ RTb$POSITION, span=(Span/nrow(RTb)))		# Smooth replication timing data
			loessProf$fitted = loessProf$fitted * 1.59 / IQR(loessProf$fitted)		# Set IQR of smoothed y-values to 1.59 

			loessSlope = NULL
			aboveGap   = (c(0,RTb$POSITION)-c(RTb$POSITION,0)) < -20e3				# Index of positions with > 20kb gaps		
			loessSlope = (diff(loessProf$fitted)/diff(RTb$POSITION))[!aboveGap]		# RT slope for positions with no gaps
	
			df = data.frame(loessSlope, RT=loessProf$fitted[!aboveGap], POS=RTb$POSITION[!aboveGap])
			df = df[1:(nrow(df)-2),]
			
			df$highSlope = df$loessSlope > Slope_E; df$highSlope[1] = FALSE
			df$lowSlope  = df$loessSlope < -Slope_E; df$lowSlope[1] = FALSE
				
			df$highSlopeL = df$loessSlope > Slope_L; df$highSlopeL[1] = FALSE
			df$lowSlopeL = df$loessSlope < -Slope_L; df$lowSlopeL[1] = FALSE
		
			df$highStart = FALSE;	df$highEnd = FALSE
			df$lowStart  = FALSE;	df$lowEnd  = FALSE
				
			# For each class, define start and end positions	
			df$highEnd	<- (c(F, df$highSlope) & !c(df$highSlope, F))[1:nrow(df)]
			df$lowStart  <- (c(df$lowSlope, F) & !c(F, df$lowSlope))[1:nrow(df)]
		
			df$highStart <- (c(df$highSlopeL, F) & !c(F, df$highSlopeL))[1:nrow(df)]
			df$lowEnd	<- (c(F, df$lowSlopeL) & !c(df$lowSlopeL, F))[1:nrow(df)]
		
			HE = NULL; HE = grep(TRUE, df$highEnd); HS = 1e10
			for(i in length(HE):1){
				if(HE[i] < HS[length(HS)]){
					HS = c(HS, max(grep(TRUE, df$highStart)[grep(TRUE, df$highStart) < HE[i]]))
				}else{
					HE[i] = 0
				}
			}
			HE = subset(HE, HE != 0); HS = subset(HS, HS != 1e10); HS = HS[length(HS):1]
		
			LS = NULL; LS = grep(TRUE, df$lowStart); LE = 0
			for(i in 1:length(LS)){
				if(LS[i] > LE[length(LE)]){
					LE = c(LE, min(grep(TRUE, df$lowEnd)[grep(TRUE, df$lowEnd) > LS[i]]))
				}else{
					LS[i] = 0
				}
			}
			LS = subset(LS, LS != 0); LE = subset(LE, LE != 0)
						
			# Find positions of high and low slope boundaries, and the size of the transition regions between them	
			aH = df$POS[HS]
			bH = df$POS[HE]
			aL = df$POS[LS]
			bL = df$POS[LE]	
		
			# Find RT at the positions described above, and deltaRT for boundary start/stops
			aRTh = df$RT[HS]
			bRTh = df$RT[HE]
			aRTl = df$RT[LS]
			bRTl = df$RT[LE]
		
			# Data collection for TTRs
 			z1 = data.frame(CHR=chr, highStart=aH, highEnd=bH, highRTstart=aRTh, highRTend=bRTh)
			z2 = data.frame(CHR=chr, lowStart=aL, lowEnd=bL, lowRTstart=aRTl, lowRTend=bRTl)	
		
			# Filter boundaries too close to chromosome ends
			z1 = z1[z1$highEnd <= as.numeric(endPoints[endPoints$CHR==chr,3])-Ends & z1$highEnd >= as.numeric(endPoints[endPoints$CHR==chr,2])+Ends,]
			z2 = z2[z2$lowStart >= as.numeric(endPoints[endPoints$CHR==chr,2])+Ends & z2$lowStart <= as.numeric(endPoints[endPoints$CHR==chr,3])-Ends,]
		
			Probes = NULL; Probes <- lapply(1:dim(z1)[1], function(i){
				length(RTb$POSITION[RTb$POSITION >= z1$highStart[i] & RTb$POSITION <= z1$highEnd[i]])
			})
			Probes <- do.call(rbind, Probes)
			z1 = cbind(z1,Probes)
		
			Probes = NULL; Probes <- lapply(1:dim(z2)[1], function(i){
				length(RTb$POSITION[RTb$POSITION >= z2$lowStart[i] & RTb$POSITION <= z2$lowEnd[i]])
			})
			Probes <- do.call(rbind, Probes)
			z2 = cbind(z2,Probes)
			
			# Combine TTR data 
			highTTR = rbind(highTTR, z1)
			lowTTR  = rbind(lowTTR, z2)
		}
		
		# Add columns to indicate the source dataset (column number in data table x, SOURCE) and orientation (FLIP = T for "high" or positive slope TTRs) of each TTR
		Rv1 = NULL; Rv1 = data.frame(SOURCE=ct,CHR=highTTR$CHR,POSITION=highTTR$highEnd,FLIP=T,SIZE=(highTTR$highEnd-highTTR$highStart),PROBES=highTTR$Probes,RT_CHANGE=(highTTR$highRTend-highTTR$highRTstart),RT=highTTR$highRTend)
		Fw1 = NULL; Fw1 = data.frame(SOURCE=ct,CHR=lowTTR$CHR,POSITION=lowTTR$lowStart,FLIP=F,SIZE=(lowTTR$lowEnd-lowTTR$lowStart),PROBES=lowTTR$Probes,RT_CHANGE=(lowTTR$lowRTstart-lowTTR$lowRTend),RT=lowTTR$lowRTstart)
		
		# Merge lists of positive and negative slope TTRs and order them by chromosome
		TTRi = NULL; TTRi = rbind(Rv1, Fw1); TTRo = NULL; TTRo = TTRi[order(TTRi$CHR),]
		
		# Generate an equal number and distribution of random positions
		randBounds = NULL
		for(chr in chrs){
			boundSet = subset(TTRo, CHR == chr)
			pointSet = subset(endPoints, CHR == chr)			
			chromSeq = seq.int(from=pointSet$Start, to=pointSet$End, length=((as.numeric(pointSet$End)-as.numeric(pointSet$Start))/50e3))
			randPos	 = sample(chromSeq, nrow(boundSet), replace=F)
			chrBounds = data.frame(RANDOM_POSITION=randPos)
			randBounds = rbind(randBounds, chrBounds)
		}
		
		# Add TTRs for current dataset (ct) to combined list of previous datasets
		TTRr = rbind(TTRr, cbind(TTRo, randBounds))					
	}
	
	cat("Finalizing", "\n")
	# Add column with unique identification number (ID) for each TTR
	TTRu = NULL; TTRu = cbind(1:nrow(TTRr),TTRr); names(TTRu)[1]="ID"
	
	FL = NULL; Pc = NULL; Pc = P5[(P5$DISR > Gap | P5$DISL > Gap),]

	for(chr in chrs){
		PcC = NULL; PcC = droplevels(Pc[Pc$CHR == as.character(chr),])
		mTc = NULL; mTc = TTRu[TTRu$CHR == chr,]
		if(invalid(PcC)){
			FL = rbind(FL, data.frame(ID=mTc$ID, ED=1e7))
		}else{
			ED = NULL; ED = lapply(mTc$POSITION, function(x) min(abs(x-PcC$POSITION))); ED <- do.call(rbind, ED)
			FL = rbind(FL, data.frame(ID=mTc$ID, ED=ED))
		}
	}

	FLo = NULL; FLo = FL[order(FL$ID),]
	mT = NULL; mT = cbind(TTRu,FLo[,2])
	names(mT)[ncol(mT)] = "ED"
	
	# Apply size, rt change, and data point filters to TTR calls
	TTRf = mT[mT$RT_CHANGE>=RTU_Min & mT$SIZE>=Length_Min & mT$SIZE<=Length_Max & mT$ED>Gap_Dis & mT$PROBES>=Data_Points,]

	return(list(TTRu,TTRf)) 			# TTRu is unfiltered and TTRf is filtered
}
