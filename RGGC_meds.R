# load libraries
library(survival)
library(data.table)
library(entropy)
library(pracma)

minCBGperAdmission<-4
ageWindow<-5 # (years)
medianWindow<-1 # (mmol/L)
admissionDurationWindow<-1 # (days)
diabetesDurationWindow<-4 # (years)

 
## functions

drugsPrePost <- function(ID,admissionDate,admissionDuration,timeWindow) {
	# ID<-admissionsWithDrugData[27559,]$ID
	# admissionDate<-admissionsWithDrugData[27559,]$dateplustime1
	# admissionDuration<-admissionsWithDrugData[27559,]$admissionDuration

	timeWindowSeconds<-timeWindow*(60*60*24*(365.25/12))
	dischargeDate<-admissionDate+admissionDuration
	
	drugInfo<-medsDataDT["PatId" == ID]
	drugInfo<-drugInfo[order(drugInfo$dateplustime1),]
	
	preAdmissionSUB<-drugInfo[which(drugInfo$dateplustime1>(admissionDate-timeWindowSeconds) & drugInfo$dateplustime1<admissionDate)]
	if (nrow(preAdmissionSUB)>0) {
		preSU<-ifelse(preAdmissionSUB$BNFCode=="6.1.2.1" | preAdmissionSUB$BNFCode=="06.01.02.01" | substr(preAdmissionSUB$BNFCode,1,7)=="0601021" | preAdmissionSUB$MatchDrugName=="Gliclazide" | preAdmissionSUB$MatchDrugName=="Glibenclamide" | preAdmissionSUB$MatchDrugName=="Glimepiride" | preAdmissionSUB$MatchDrugName=="Glipizide",1,0)
		preIns<-ifelse(preAdmissionSUB$BNFCode=="6.1.1.1" | preAdmissionSUB$BNFCode=="06.01.01.01" | substr(preAdmissionSUB$BNFCode,1,7)=="0601011" | preAdmissionSUB$BNFCode=="6.1.1.2" | preAdmissionSUB$BNFCode=="06.01.01.02" | substr(preAdmissionSUB$BNFCode,1,7)=="0601012" | preAdmissionSUB$MatchDrugName=="Actrapid"  | preAdmissionSUB$MatchDrugName=="Humalog" | preAdmissionSUB$MatchDrugName=="Humulin" | preAdmissionSUB$MatchDrugName=="Hypurin" | preAdmissionSUB$MatchDrugName=="Insulatard" | preAdmissionSUB$MatchDrugName=="Insulin" | preAdmissionSUB$MatchDrugName=="Insuman" | preAdmissionSUB$MatchDrugName=="Lantus",1,0)
		preMF<-ifelse(preAdmissionSUB$BNFCode=="6.1.2.2" | preAdmissionSUB$BNFCode=="06.01.02.02" | substr(preAdmissionSUB$BNFCode,1,7)=="0601022" | preAdmissionSUB$MatchDrugName=="Metformin",1,0)
		preOther<-ifelse(preAdmissionSUB$BNFCode=="6.1.2.3" | preAdmissionSUB$BNFCode=="06.01.02.03" | substr(preAdmissionSUB$BNFCode,1,7)=="0601023",1,0)
		preBB<-ifelse(preAdmissionSUB$BNFCode=="2.4" | preAdmissionSUB$BNFCode=="02.04" | substr(preAdmissionSUB$BNFCode,1,4)=="0204" | preAdmissionSUB$MatchDrugName=="Bisoprolol" | preAdmissionSUB$MatchDrugName=="Propranolol",1,0)
	}
	if (nrow(preAdmissionSUB)==0) {preSU<-0;preIns<-0;preMF<-0;preOther<-0;preBB<-0}
	
	postAdmissionSUB<-drugInfo[which(drugInfo$dateplustime1>=(dischargeDate) & drugInfo$dateplustime1<(dischargeDate+timeWindowSeconds))]
	if (nrow(postAdmissionSUB)>0) {
		postSU<-ifelse(postAdmissionSUB$BNFCode=="6.1.2.1" | postAdmissionSUB$BNFCode=="06.01.02.01" | substr(postAdmissionSUB$BNFCode,1,7)=="0601021" | postAdmissionSUB$MatchDrugName=="Gliclazide" | postAdmissionSUB$MatchDrugName=="Glibenclamide" | postAdmissionSUB$MatchDrugName=="Glimepiride" | postAdmissionSUB$MatchDrugName=="Glipizide",1,0)
		postIns<-ifelse(postAdmissionSUB$BNFCode=="6.1.1.1" | postAdmissionSUB$BNFCode=="06.01.01.01" | substr(postAdmissionSUB$BNFCode,1,7)=="0601011" | postAdmissionSUB$BNFCode=="6.1.1.2" | postAdmissionSUB$BNFCode=="06.01.01.02" | substr(postAdmissionSUB$BNFCode,1,7)=="0601012" | postAdmissionSUB$MatchDrugName=="Actrapid"  | postAdmissionSUB$MatchDrugName=="Humalog" | postAdmissionSUB$MatchDrugName=="Humulin" | postAdmissionSUB$MatchDrugName=="Hypurin" | postAdmissionSUB$MatchDrugName=="Insulatard" | postAdmissionSUB$MatchDrugName=="Insulin" | postAdmissionSUB$MatchDrugName=="Insuman" | postAdmissionSUB$MatchDrugName=="Lantus",1,0)
		postMF<-ifelse(postAdmissionSUB$BNFCode=="6.1.2.2" | postAdmissionSUB$BNFCode=="06.01.02.02" | substr(postAdmissionSUB$BNFCode,1,7)=="0601022" | postAdmissionSUB$MatchDrugName=="Metformin",1,0)
		postOther<-ifelse(postAdmissionSUB$BNFCode=="6.1.2.3" | postAdmissionSUB$BNFCode=="06.01.02.03" | substr(postAdmissionSUB$BNFCode,1,7)=="0601023",1,0)
		postBB<-ifelse(postAdmissionSUB$BNFCode=="2.4" | postAdmissionSUB$BNFCode=="02.04" | substr(postAdmissionSUB$BNFCode,1,4)=="0204" | postAdmissionSUB$MatchDrugName=="Bisoprolol" | postAdmissionSUB$MatchDrugName=="Propranolol",1,0)
	}
	if (nrow(postAdmissionSUB)==0) {postSU<-0;postIns<-0;postMF<-0;postOther<-0;postBB<-0}
	
	preSU<-ifelse(sum(preSU)>0,1,0)
	preIns<-ifelse(sum(preIns)>0,1,0)
	preMF<-ifelse(sum(preMF)>0,1,0)
	preOther<-ifelse(sum(preOther)>0,1,0)
	preBB<-ifelse(sum(preBB)>0,1,0)
	postSU<-ifelse(sum(postSU)>0,1,0)
	postIns<-ifelse(sum(postIns)>0,1,0)
	postMF<-ifelse(sum(postMF)>0,1,0)
	postOther<-ifelse(sum(postOther)>0,1,0)
	postBB<-ifelse(sum(postBB)>0,1,0)
	
	spanSU<-ifelse(preSU==1 & postSU==1,1,0)
	spanIns<-ifelse(preIns==1 & postIns==1,1,0)
	spanMF<-ifelse(preMF==1 & postMF==1,1,0)
	spanOther<-ifelse(preOther==1 & postOther==1,1,0)
	spanBB<-ifelse(preBB==1 & postBB==1,1,0)
	
	outputList<-list(preSU,preIns,preMF,preOther,preBB,spanSU,spanIns,spanMF,spanOther,spanBB)
	return(outputList)
}

# find first admission with more than one CBG
firstAdmissionWithMoreThanNCBG<-function(admissionNumberFlag,CBGinSequencePerAdmission,nCBGperAdmission,nCBGthreshold) {
	# admissionNumberFlag<-DT[which(DT$ID==101325789)]$admissionNumberFlag
	# CBGinSequencePerAdmission<-DT[which(DT$ID==101325789)]$CBGinSequencePerAdmission
	# nCBGperAdmission<-DT[which(DT$ID==101325789)]$nCBGperAdmission
	
	firstMultipleCBGAdmission<-ifelse(CBGinSequencePerAdmission==1 & nCBGperAdmission>nCBGthreshold,1,0)
	cumsulsum<-cumsum(firstMultipleCBGAdmission)
	selectRow<-ifelse(firstMultipleCBGAdmission==1 & cumsulsum==1,1,0)
	
	return(selectRow)
	
}

addSurvivalData<-function(inputFrame,plotTitle,test,threshold) {
	
	SurvivalData<-inputFrame
	
	DaySeconds<-(60*60*24)
	shortCensorPeriodStartDay  <- DaySeconds*90
	shortCensorPeriodEndDay    <- DaySeconds*10000
	
	lastDOD<-max(SurvivalData$deathDateUnix)
	SurvivalData$dateOfDischarge<-SurvivalData$dateplustime1+SurvivalData$admissionDuration
	SurvivalData$timeToDeath<-ifelse(SurvivalData$deathEvent==1,(SurvivalData$deathDateUnix-SurvivalData$dateOfDischarge),0)
#		SurvivalData$timeToDeath<-SurvivalData$timeToDeath/DaySeconds
	SurvivalData$timeToDeathInterval<-ifelse(SurvivalData$deathEvent==0,(lastDOD-SurvivalData$dateOfDischarge),SurvivalData$timeToDeath)
	SurvivalData$timeToDeathInterval[is.na(SurvivalData$timeToDeathInterval)]<-0; SurvivalData<-subset(SurvivalData,timeToDeathInterval>0)
#		SurvivalData$timeToDeathInterval<-SurvivalData$timeToDeathInterval/DaySeconds
		
	SurvivalData$shortDeathEvent <- SurvivalData$deathEvent
	SurvivalData$shortDeathEvent <- ifelse(SurvivalData$deathEvent==1 & SurvivalData$timeToDeath>=(shortCensorPeriodStartDay) & SurvivalData$timeToDeath<(shortCensorPeriodEndDay),1,0)	
	
	SurvivalData$sexDigit<-ifelse(nchar(SurvivalData$charID==9),as.numeric(substr(SurvivalData$charID,8,8)),as.numeric(substr(SurvivalData$charID,9,9)))
	SurvivalData$sexNumber<-ifelse(SurvivalData$sexDigit%%2==0,1,0)
	SurvivalData$sex<-factor(1*(SurvivalData$sexNumber <1),levels=0:1,labels=c("F","M"))
	
#	mfitAge50<-survfit(Surv(timeToDeathInterval, shortDeathEvent) ~ 1, data = SurvivalData)
#	plot(mfitAge50,mark.time=F)
	
	fdata<-SurvivalData
	
	# remove survival data prior to shortcensorperiodstartday to allow a true assessment of post 30-day survival etc
	fdata$timeToDeathInterval<-fdata$timeToDeathInterval-shortCensorPeriodStartDay
	
	quantileIQR<-as.table(quantile(fdata$IQR))
	quantileMedian<-as.table(quantile(fdata$median))
	
	if (test=="minGlu") { fdata$group<-factor(1+ 1*(fdata$minGlu < threshold ),levels=1:2,labels=c(paste("minGlu>=",threshold,sep=""),paste("minGlu<",threshold,sep=""))) }
	if (test=="maxGlu") { fdata$group<-factor(1+ 1*(fdata$maxGlu < threshold ),levels=1:2,labels=c(paste("maxGlu>=",threshold,sep=""),paste("maxGlu<",threshold,sep=""))) }
	if (test=="IQR")    { fdata$group<-factor(1+ 1*(fdata$IQR < quantileIQR[threshold] ),levels=1:2,labels=c(paste("IQR>=",threshold,sep=""),paste("IQR<",threshold,sep=""))) }
	if (test=="median") { fdata$group<-factor(1+ 1*(fdata$IQR < quantileMedian[threshold] ),levels=1:2,labels=c(paste("median>=",threshold,sep=""),paste("median<",threshold,sep=""))) }

	if (test=="BB")     { fdata$group<-factor(1+ 1*(fdata$spanBB == 1 ),levels=1:2,labels=c(paste("spanBB==0",sep=""),paste("spanBB==1",sep=""))) }

	if (test=="Drug")   { fdata$group<-factor(1+ 1*(fdata$caseFlag == 1 ),levels=1:2,labels=c(paste("control",sep=""),paste("case",sep=""))) }
	
		print(fdata[1:5,])
	
	fdata$logAD<-log(fdata$admissionDuration)
	fdata$age2 <- cut(fdata$age, c(10,20, 30,40, 50,60,70, 80, 110), labels = c(paste(c(10,20,30,40,50,60,70), c(19,29,39,49,59,69,79), sep='-'), "80+"))
	fdata$AD2 <- cut(fdata$logAD, c(4,6, 8,10, 12,14,16, 18, 20), labels = c(paste(c(10,20,30,40,50,60,70), c(19,29,39,49,59,69,79), sep='-'), "80+"))
	fdata<-subset(fdata,admissionDuration>0)
	
	sfit1 <- survfit(Surv(timeToDeathInterval, shortDeathEvent) ~ group, fdata)
	
		index <- with(fdata, cbind(as.numeric(age2), as.numeric(sex), as.numeric(group)))
		tab1 <- with(fdata, table(age2, sex))/ nrow(fdata)
		tab2 <- with(fdata, table(age2, sex, group))/nrow(fdata)
		tab3 <- with(fdata, table(group)) / nrow(fdata)
		rwt <- rep(tab1,2)/tab2

		indexAD <- with(fdata, cbind(as.numeric(AD2), as.numeric(sex), as.numeric(group)))
		tab1 <- with(fdata, table(AD2, sex))/ nrow(fdata)
		tab2 <- with(fdata, table(AD2, sex, group))/nrow(fdata)
		tab3 <- with(fdata, table(group)) / nrow(fdata)
		rwtAD <- rep(tab1,2)/tab2
	
		# generate weighting for both admission duration and age
		indexBoth <- with(fdata, cbind(as.numeric(age2), sex, as.numeric(AD2), as.numeric(group)))
		tab1 <- with(fdata, table(age2, sex, AD2))/ nrow(fdata)
		tab2 <- with(fdata, table(age2, sex, AD2, group))/nrow(fdata)
		tab3 <- with(fdata, table(group)) / nrow(fdata)
		rwtBoth <- rep(tab1,2)/tab2
		
	plotName<-paste(plotTitle,"\nn=",nrow(fdata),sep="")
		fdata$rwt <- rwt[index]  # add per subject weights to the data set
		fdata$rwtAD <- rwtAD[indexAD]  # add per subject weights to the data set
		fdata$rwtBoth <- rwtBoth[indexBoth]
	
	sfit3 <- survfit(Surv(timeToDeathInterval, shortDeathEvent) ~ group, data=fdata, weight=rwtBoth)
	plot(sfit1, mark.time=F, col=c(1,2,4), lty=1, lwd=2, xscale=365.25, xlab="time", ylab="Survival",main=plotName)
		
	lines(sfit3,  mark.time=F, col=c(1,2,4), lty=2, lwd=1, xscale=365.25)
	
	sfit4.coxph.old<-coxph(Surv(timeToDeathInterval,  shortDeathEvent) ~ group + cluster(pairNumber), data=fdata)
	
	
		# test for proportions for survival curve weighted for admission duration
		
		chiTest<-prop.test(c(as.numeric(summary(sfit3)$table[, "events"][1]),as.numeric(summary(sfit3)$table[, "events"][2])),c(as.numeric(summary(sfit3)$table[, "n.start"][1]),as.numeric(summary(sfit3)$table[, "n.start"][2])),correct=F)
		# correct chi test - apply corrected survival proportions (from weighted test) to original group sizes
		correctedChiTest<-prop.test(summary(sfit1)$table[,"n.start"] * chiTest$estimate , summary(sfit1)$table[,"n.start"])

		sfit3.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ group, data=fdata, weight=rwtBoth)
			pVal <- summary(sfit3.coxph)$coef[,5]; HR <- round(exp(coef(sfit3.coxph)),2)
			legendText <- paste("weighted surv analysis: p=",pVal," | HR=",HR,"\ncoxph multivariable analysis: p=",summary(sfit4.coxph.old)$coef[,6][1],
			" | HR=",round(exp(coef(sfit4.coxph.old))[1],2),"\nprop.test (unweighted) p val: ",correctedChiTest$p.val,"\nprops: ",as.data.frame(correctedChiTest$estimate)[1,1]," ",as.data.frame(correctedChiTest$estimate)[2,1],sep="")
		summarySurvfit <- summary(sfit3); legendNames <- row.names(summarySurvfit$table)
		legend("bottomleft",c(legendNames),lty=1:6,col=c("black","red","blue","green","orange","purple"),cex=0.8); legend("topright",legendText,cex=0.6,bty='n')
		
	
}

## case control function
controlSelection<-function(cases,controlPool) {
	# 
	cases<-subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & (preSU==0 & preIns==1 & preMF==1))
	# 
	controlPool<-subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & (preSU==1 & preIns==0 & preMF==1))
	
	#cases<-subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & IQR>=4.9)
	#controlPool<-subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & IQR<4.9)
	
	cases<-cases[diagnosisDateUnix!=(-2208988800),];controlPool<-controlPool[diagnosisDateUnix!=(-2208988800),]
	cases$diabetesDuration<-(cases$dateplustime1-cases$diagnosisDateUnix)/(60*60*24*365.25)
	controlPool$diabetesDuration<-(controlPool$dateplustime1-controlPool$diagnosisDateUnix)/(60*60*24*365.25)
	
	# cases<-subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & (spanSU==1 & spanIns==0 & spanMF==1))
	# controlPool<-subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & (spanSU==0 & spanIns==0 & spanMF==1))
	controlPool$used<-0
	
	dim(cases);dim(controlPool)
	
	controlPoolDT<-data.table(controlPool)
	
	# set up flags for survival plotting
	admissionsWithDrugData$caseFlag<-0
	admissionsWithDrugData$controlFlag<-0
	
	testOutput<-as.data.frame(matrix(0,nrow=nrow(cases),ncol=15))
	colnames(testOutput)<-c("caseAge","controlAge","caseMedian","controlMedian","caseIQR","controlIQR","caseMort","controlMort","caseHypoEpisodes","controlHypoEpisodes","caseAdmissionDurationDays","controlAdmissionDurationDays","numberOfMatches","caseDiabetesDuration","controlDiabetesDuration")
	
	#	variable of interest in this loop IQR
	for (i in seq(1,nrow(cases),1)) {
		if (i%%100==0) {print(i)}
		
		availableControlPool<-controlPool[which(controlPool$used==0),]
		testOutput$caseAge[i]<-cases$age[i]
		testOutput$caseMedian[i]<-cases$medianGlu[i]
		testOutput$caseIQR[i]<-cases$IQR[i]
		testOutput$caseMort[i]<-cases$deathEvent[i]
		testOutput$caseHypoEpisodes[i]<-cases$ID_ADMISSIONhypoEpisodes4.60[i]
		testOutput$caseAdmissionDurationDays[i]<-cases$admissionDurationDays[i]
		testOutput$caseDiabetesDuration[i]<-cases$diabetesDuration[i]
		
		# flag as a case in admissionsWithDrugs for Survival Plotting
		# admissionsWithDrugData$caseFlag<-ifelse(cases$ID[i]==admissionsWithDrugData$ID,1,admissionsWithDrugData$caseFlag)
		
		caseAge<-cases$age[i]
		caseMedian<-cases$medianGlu[i]
		caseDuration<-round(cases$admissionDurationDays[i],1)
		case_nCBG<-cases$nCBGperAdmission[i]
		caseDiabetesDuration<-cases$diabetesDuration[i]
		
		# 
		matchingSUB<-availableControlPool[which(availableControlPool$age>(caseAge-(ageWindow/2)) & 
		availableControlPool$age<(caseAge+(ageWindow/2)) & 
		availableControlPool$medianGlu>(caseMedian-(medianWindow/2)) & 
		availableControlPool$medianGlu<(caseMedian+(medianWindow/2)) & 
		round(availableControlPool$admissionDurationDays,1)>(caseDuration-(admissionDurationWindow/2)) & 
		round(availableControlPool$admissionDurationDays,1)<(caseDuration+(admissionDurationWindow/2)) &
		availableControlPool$diabetesDuration>(caseDiabetesDuration-(diabetesDurationWindow/2)) &
		availableControlPool$diabetesDuration<(caseDiabetesDuration+(diabetesDurationWindow/2))),]
		# matchingSUB<-availableControlPool[which(availableControlPool$age>(caseAge-2.5) & availableControlPool$age<(caseAge+2.5) & availableControlPool$medianGlu>(caseMedian-0.1) & availableControlPool$medianGlu<(caseMedian+0.1) & nCBGperAdmission>(case_nCBG-3) & nCBGperAdmission<(case_nCBG+3)),]
		
		if (nrow(matchingSUB)>0) {
			randomRow<-sample(1:nrow(matchingSUB),1)
			testOutput$numberOfMatches[i]<-nrow(matchingSUB)
			testOutput$controlAge[i]<-matchingSUB$age[randomRow]
			testOutput$controlMedian[i]<-matchingSUB$medianGlu[randomRow]
			testOutput$controlIQR[i]<-matchingSUB$IQR[randomRow] # take 1st row ? need to randomise row taken
			testOutput$controlMort[i]<-matchingSUB$deathEvent[randomRow]
			testOutput$controlHypoEpisodes[i]<-matchingSUB$ID_ADMISSIONhypoEpisodes4.60[randomRow]
			testOutput$controlAdmissionDurationDays[i]<-matchingSUB$admissionDurationDays[randomRow]
			testOutput$controlDiabetesDuration[i]<-matchingSUB$diabetesDuration[randomRow]
			
			controlPool$used<-ifelse(controlPool$ID==matchingSUB$ID[randomRow],1,controlPool$used)
			
			admissionsWithDrugData$controlFlag<-ifelse(matchingSUB$ID[randomRow]==admissionsWithDrugData$ID,i,admissionsWithDrugData$controlFlag)
			admissionsWithDrugData$caseFlag<-ifelse(cases$ID[i]==admissionsWithDrugData$ID,i,admissionsWithDrugData$caseFlag)
			
		}
			
	}
	
		# write out data for later plotting/analysis
		fileName<-paste("PRE-includingNameSearch-_SU",cases$preSU[1],"-Ins",cases$preIns[1],"-MF",cases$preMF[1],"_vs_SU",controlPool$preSU[1],"-Ins",controlPool$preIns[1],"-MF",controlPool$preMF[1],"_",minCBGperAdmission+1,"CBGmin","_ageWindow-",ageWindow,"_medianWindow-",medianWindow,"_admissionDurationWindow-",admissionDurationWindow,"_diabetesDurationWindow-",diabetesDurationWindow,sep="")	
		# fileName<-paste("SU",cases$spanSU[1],"-Ins",cases$spanIns[1],"-MF",cases$spanMF[1],"_vs_SU",controlPool$spanSU[1],"-Ins",controlPool$spanIns[1],"-MF",controlPool$spanMF[1],"_",minCBGperAdmission+1,"CBGmin","_ageWindow-",ageWindow,"_medianWindow-",medianWindow,"_admissionDurationWindow-",admissionDurationWindow,sep="")
		testOutputwriteName <- paste("../GlCoSy/GGCdrug/testOutput",fileName,".csv",sep="")
		write.table(testOutput,file=testOutputwriteName,sep=",",append=FALSE,col.names=TRUE)
		admissionsWithDrugDatawriteName <- paste("../GlCoSy/GGCdrug/admissionsWithDrugData",fileName,"-2CBGmin-0.4median.csv",sep="")
		write.table(admissionsWithDrugData,file=admissionsWithDrugDatawriteName,sep=",",append=FALSE,col.names=TRUE)
	
	
	# testOutput<-read.csv(testOutputwriteName, header=TRUE , sep="," , row.names=NULL)
	# admissionsWithDrugData<-read.csv(admissionsWithDrugDatawriteName, header=TRUE , sep="," , row.names=NULL)
	

	testOutputCaseControl<-subset(testOutput,controlIQR>0 & caseIQR>0); dim(testOutputCaseControl)
	ttest<-t.test(log(testOutputCaseControl$caseIQR),log(testOutputCaseControl$controlIQR),paired=T)
		ttest
	quantile(testOutputCaseControl$caseIQR); quantile(testOutputCaseControl$controlIQR)
	mortalityProp<-prop.test(c(sum(testOutputCaseControl$caseMort),sum(testOutputCaseControl$controlMort)),c(nrow(testOutputCaseControl),nrow(testOutputCaseControl)))
		mortalityProp
	hypoProp<-prop.test(c(sum(testOutputCaseControl$caseHypoEpisodes),sum(testOutputCaseControl$controlHypoEpisodes)),c(nrow(testOutputCaseControl),nrow(testOutputCaseControl)))
		hypoProp
	sum(testOutputCaseControl$caseHypoEpisodes)/sum(testOutputCaseControl$caseAdmissionDurationDays)
	sum(testOutputCaseControl$controlHypoEpisodes)/sum(testOutputCaseControl$controlAdmissionDurationDays)
	ttestAdmissionDuration<-t.test(log(testOutputCaseControl$caseAdmissionDurationDays),log(testOutputCaseControl$controlAdmissionDurationDays))
		ttestAdmissionDuration
	quantile(testOutputCaseControl$caseAdmissionDurationDays); quantile(testOutputCaseControl$controlAdmissionDurationDays)
	
	ttestAge<-t.test(testOutputCaseControl$caseAge,testOutputCaseControl$controlAge,paired=T)
		ttestAge
		quantile(testOutputCaseControl$caseAge);quantile(testOutputCaseControl$controlAge)
		
	ttestMedian<-t.test(log(testOutputCaseControl$caseMedian),log(testOutputCaseControl$controlMedian),paired=T)
		ttestMedian
	quantile(testOutputCaseControl$caseMedian);quantile(testOutputCaseControl$controlMedian)
	
	ttestDiabetesDuration<-t.test(log(subset(testOutputCaseControl,caseDiabetesDuration>0 & controlDiabetesDuration>0)$caseDiabetesDuration),log(subset(testOutputCaseControl,caseDiabetesDuration>0 & controlDiabetesDuration>0)$controlDiabetesDuration),paired=T)
		ttestDiabetesDuration
		quantile(testOutputCaseControl$caseDiabetesDuration);quantile(testOutputCaseControl$controlDiabetesDuration)
	
		
	
	# mark pairs matched in order that the cluster argument can be used to perform a matched survival analysis (included in the cox prop hazard p value)
	survivalSUB<-subset(admissionsWithDrugData,caseFlag>0 | controlFlag>0)
	survivalSUB$pairNumber<-ifelse(survivalSUB$caseFlag>0,survivalSUB$caseFlag,ifelse(survivalSUB$controlFlag>0,survivalSUB$controlFlag,0))
	survivalSUB$caseFlag<-ifelse(survivalSUB$caseFlag>0,1,0)
	survivalSUB$controlFlag<-ifelse(survivalSUB$controlFlag>0,1,0)
	
	plotPasteName<-paste(fileName,".pdf",sep="")
	plotfilename <- paste("../GlCoSy/GGCdrug/",plotPasteName,sep="")
	pdf(plotfilename, width=16, height=9)
		plotName<-paste(fileName," | >1 CBG per admission",sep="")
		addSurvivalData(survivalSUB,plotName,"Drug",0)
	dev.off()

}

# attempt to investigate WHEN hypos occur during an admission divided by treatment type - doesnt really work
hypoOccurenceDuringAdmission<-function(testSet,longestAdmission) {
	# needs: admisisonsWithDrugData & DT
	# example testSet:
	# testSet<-subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 1 Diabetes Mellitus" & (preSU==1 | preIns==1 & preMF==1))
	# testSet<-subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 1 Diabetes Mellitus")
	
	summaryOutputName <- paste("../GlCoSy/output/timePointOfHypoCBG.csv",sep="")
	unlink(summaryOutputName)

	whenAreTheHypos<-function(IDtest,admissionNumberFlagtest,hypoThreshold) {
		# IDtest<-testSet$ID[1]; admissionNumberFlagtest<-testSet$admissionNumberFlag[1]
		DTsub<-DT[ID==IDtest & admissionNumberFlag==admissionNumberFlagtest]
		DTsub<-DTsub[order(DTsub$dateplustime1),]
		DTsub$timeToHypoCBG<-ifelse(DTsub$yyyy<hypoThreshold,DTsub$dateplustime1-DTsub$dateplustime1[1],-100)
		reportFrame<-data.frame(DTsub$yyyy,DTsub$timeToHypoCBG,DTsub$admissionDuration,DTsub$ID,DTsub$admissionNumberFlag,DTsub$nCBGperAdmission)
					
			write.table(reportFrame,file=summaryOutputName,sep=",",append=T,col.names=F)
		
		outputList<-list(DTsub$yyyy,DTsub$timeToHypoCBG)
		return(DTsub$ID[1])
		
	}
	# write out list of hypos and time to hypos
	testSet[, c("returnVal") := whenAreTheHypos(ID,admissionNumberFlag,4) , by=.(ID,admissionNumberFlag)]
	
	##### read list back in and plot
	whenHypoData<-read.csv(summaryOutputName, header=T , sep="," , row.names=NULL)
	colnames(whenHypoData)<-c("rowname","yyyy","timeToHypoCBG","admissionDuration","ID","admissionNumberFlag","nCBGperAdmission")
		whenHypoData$admissionDurationDays<-whenHypoData$admissionDuration/(60*60*24); whenHypoData<-whenHypoData[whenHypoData$timeToHypoCBG>=0,]
	
	whenHypoDataDT<-data.table(whenHypoData)
	
		timeToFirstHypoOutputName <- paste("../GlCoSy/output/timeToFirstHypo.csv",sep="")
		unlink(timeToFirstHypoOutputName)
	
			timeToFirstHypoCBG<-function(timeToHypoCBG,admissionDuration) {
				firstHypo<-timeToHypoCBG[1]
				admissionDurationReport<-admissionDuration[1]
			
				outputFrame<-data.frame(firstHypo,admissionDurationReport)
		
				write.table(outputFrame,file=timeToFirstHypoOutputName,sep=",",append=T,col.names=F)
			}
			
		whenHypoDataDT[,(timeToFirstHypoCBG(timeToHypoCBG,admissionDuration)), by=.(ID,admissionNumberFlag)]
		whenFirstHypo<-read.csv(timeToFirstHypoOutputName, header=F , sep="," , row.names=NULL); colnames(whenFirstHypo)<-c("null","firstHypo","admissionDuration")
		plot(whenFirstHypo$admissionDuration,whenFirstHypo$firstHypo,pch=16,cex=0.3)
		plot(whenFirstHypo$admissionDuration/(60*60*24),(whenFirstHypo$firstHypo/(60*60*24))/(whenFirstHypo$admissionDuration/(60*60*24)),pch=16,cex=0.3,xlim=c(0,60),ylim=c(0,1))
		boxplot((whenFirstHypo$firstHypo/(60*60*24))/(whenFirstHypo$admissionDuration/(60*60*24)) ~ cut(whenFirstHypo$admissionDuration/(60*60*24),breaks=seq(0,max(whenFirstHypo$admissionDuration/(60*60*24)),1)),xlim=c(0,30))
	
	# variety of plots ? useful
	for (duration in seq(0,longestAdmission,1)) {
		plotSub<-whenHypoData[whenHypoData$admissionDurationDays>=duration & whenHypoData$admissionDurationDays<(duration+1),]
		plotSub$timeToHypoCBGDays<-plotSub$timeToHypoCBG/(60*60*24)
		histPlotSub<-hist(plotSub$timeToHypoCBGDays,breaks=seq(0,duration+1,0.1),plot=F)
		
		ID_number<-length(unique(plotSub$ID))
		totalEvents<-sum(histPlotSub$counts)
		
		#if (duration==0) { plot(c(1:length(histPlotSub$counts)),(histPlotSub$counts/totalEvents),pch=16,cex=0.2,col=duration+1) }
		#if (duration!=0) { points(c(1:length(histPlotSub$counts)),(histPlotSub$counts/totalEvents),pch=16,cex=0.2,col=duration+1) }
		#lines(c(1:length(histPlotSub$counts)),(histPlotSub$counts/totalEvents),pch=16,cex=0.2,col=duration+1)
	
		#if (duration==0) { plot(c(1:length(histPlotSub$counts)),(histPlotSub$counts/ID_number),pch=16,cex=0.2,col=duration+1) }
		#if (duration!=0) { points(c(1:length(histPlotSub$counts)),(histPlotSub$counts/ID_number),pch=16,cex=0.2,col=duration+1) }
		#lines(c(1:length(histPlotSub$counts)),(histPlotSub$counts/ID_number),pch=16,cex=0.2,col=duration+1)
	
		if (duration==0) { plot(histPlotSub$mids,(cumsum(histPlotSub$counts)/ID_number),pch=16,cex=0.2,col=duration+1,xlim=c(0,longestAdmission),ylim=c(0,10)) }
		if (duration!=0) { points(histPlotSub$mids,(cumsum(histPlotSub$counts)/ID_number),pch=16,cex=0.2,col=duration+1) }
		lines(histPlotSub$mids,(cumsum(histPlotSub$counts)/ID_number),pch=16,cex=0.2,col=duration+1)
	
	
	}
	
	cumPlot<-whenHypoData[order(whenHypoData$timeToHypoCBG),]
	cumPlot$hypoFlag<-1
	cumPlot$cumsum<-cumsum(cumPlot$hypoFlag)
	plot(cumPlot$timeToHypoCBG,cumPlot$cumsum)
	
}

# load pre calculated data.table per admission data
summaryOutputName <- paste("../GlCoSy/output/DTwithPerIDdata.csv",sep=""); DT<-read.csv(summaryOutputName, header=TRUE , sep="," , row.names=NULL)
DT<-as.data.frame(DT)
	# extract unix epoch date for deathDate
	datetime = DT$DeathDate
	dateExtract = substr(datetime,1,10)						# extract time only from date / time information
	dateplustime = strptime(datetime,"%Y-%m-%d") 	# convert date and time to useable format
	dateplustime1 <- as.numeric(dateplustime) # convert date time data to numerical values (absolute value is in seconds)
	DT$deathDateUnix<-dateplustime1; DT$deathDateUnix[is.na(DT$deathDateUnix)]<-0 # ? an hour out - may need adjusted

	# extract unix epoch date of diagnosis
	datetime = DT$DateOfDiagnosisDiabetes_Date
	dateExtract = substr(datetime,1,10)						# extract time only from date / time information
	dateplustime = strptime(datetime,"%Y-%m-%d") 	# convert date and time to useable format
	dateplustime1 <- as.numeric(dateplustime) # convert date time data to numerical values (absolute value is in seconds)
	DT$diagnosisDateUnix<-dateplustime1; DT$diagnosisDateUnix[is.na(DT$diagnosisDateUnix)]<-0 # ? an hour out - may need adjusted
	

DT<-data.table(DT)

inPutName<-paste("../GlCoSy/source/GGC_meds.txt",sep=",")
medsExSUData <- read.csv(inPutName, header=TRUE , sep="," , row.names=NULL)

inPutName<-paste("../GlCoSy/source/GGC_SU.txt",sep=",")
SUData <- read.csv(inPutName, header=TRUE , sep="," , row.names=NULL)

medsData<-rbind(medsExSUData,SUData)

# extract unix epoch date for prescriptions
datetime = medsData$PrescriptionDateTime
dateExtract = substr(datetime,1,10)						# extract time only from date / time information
dateplustime = strptime(datetime,"%Y-%m-%d") 	# convert date and time to useable format
dateplustime1 <- as.numeric(dateplustime) # convert date time data to numerical values (absolute value is in seconds)
medsData$dateplustime1<-dateplustime1 # ? an hour out - may need adjusted

# cut out all prescriptions before 1-1-2009 for now - as linking to admissions
medsData<-subset(medsData,dateplustime1>=1230768000)
# cut out all prescriptions after the last death for now - as linking to admissions
medsData<-subset(medsData,dateplustime1<=1420934400)

## save out / read back medsData
# medsDataName<-paste("../GlCoSy/source/medsData.csv",sep=",")
# write.table(medsData,file=medsDataName,append=F,sep=",",col.names=T)
# medsData <- read.csv(medsDataName, header=TRUE , sep="," , row.names=NULL)

medsDataDT<-data.table(medsData)

## drug analysis

DT[, interestAdmissionFlag := firstAdmissionWithMoreThanNCBG(admissionNumberFlag,CBGinSequencePerAdmission,nCBGperAdmission,minCBGperAdmission) , by=.(ID)]
# select interest admissions as the first admission with more than one CBG
interestAdmissions<-DT[which(DT$interestAdmissionFlag==1)]
# interestAdmissions<-DT[which(DT$admissionNumberFlag==1&DT$CBGinSequence==1&DT$admissionDuration>0)]

match_admissionDrugs <- match(interestAdmissions$ID,medsData$PatId)
match_admissionDrugs[is.na(match_admissionDrugs)]<-0
interestAdmissions$match_admissionDrugs <- ifelse(match_admissionDrugs>1,1,0)

admissionsWithDrugData<-interestAdmissions[which(interestAdmissions$match_admissionDrugs==1)]
admissionsWithDrugData$deathEvent<-ifelse(admissionsWithDrugData$deathDateUnix>0,1,0)

admissionsWithDrugData[, c("preSU","preIns","preMF","preOther","preBB","spanSU","spanIns","spanMF","spanOther","spanBB") := drugsPrePost(ID,dateplustime1,admissionDuration,4) , by=ID]

# write admissionsWithDrugData to file :: 
# OPfilename <- paste("../GlCoSy/output/admissionsWithDrugData.csv",sep=""); write.table(admissionsWithDrugData,file=OPfilename,append=F,sep=",",col.names=T)









### reporting plots/frames
dim(subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & spanSU==1))
dim(subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & (spanSU==1 | spanIns==1))
dim(subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & spanMF==1))
dim(subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & spanIns==1))

dim(subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & spanSU==1 & spanBB==1 & minGlu<4))


x<-boxplot(subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & spanIns==0 & spanSU==1)$minGlu ~ subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & spanIns==0 &spanSU==1)$spanBB)
x<-boxplot(subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & spanIns==0 & spanMF==1)$IQR ~ subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & spanIns==0 & spanMF==1)$spanSU)
x<-boxplot(subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & spanIns==0 & spanMF==1)$timeGateQuartileCoefficientOfDispersion ~ subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & spanIns==0 & spanMF==1)$spanSU)

x<-boxplot(subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & spanIns==0 & spanMF==1)$medianGlu ~ subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & spanIns==0 & spanMF==1)$spanSU)

x<-boxplot(subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & spanIns==0)$maxGlu ~ subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & spanIns==0)$spanSU)

testSet<-subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & spanSU==1)



# admission durations seem similar: quantile(subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & spanIns==0 & spanSU==0)$admissionDuration)/(60*60*24)

## to do
# add date of admission per episode
# pull out admission summary per ID for the admission of interest - ?1st admission or ?1st admission with more than 1 CBG
# (maybe try and pass yyyy for this admission along with dateplustime1/time1 so that characteristics at particular times can be looked at)
# for each admission create columns for each class of agent of interest for pre/post, and then indicate whether on with a 1/0 code. if 1,1 then on pre/post and therefore assume on during admission.




