# this is a subset of the code in RGGC_meds.R to generate mortality plots for specific subgroups with 

# load libraries
library(survival)
library(data.table)
library(entropy)
library(pracma)

minCBGperAdmission<-4
ageWindow<-10 # (years)
medianWindow<-0.5 # (mmol/L)
admissionDurationWindow<-1 # (days)
diabetesDurationWindow<-4 # (years)
IQRwindow<-0.25 # (mmol/L)
eGFRwindow<-20

loadeGFR_minus90_to_plus1=0

controlChartFunction<-function(admissionOfInterest) {
	# match set according to age,median,admissionDuration,diabetesDuration - not IQR intially
	# define 2 and 3 sigma limits for control chart from this match set
	# pull out CBG profile for the admission of interest
	# plot and characterise profile with regard to the 2 and 3 sigma limits
	# associate crossing these limits with outcome, and with probability of readmission/hypoglycaemia and 1y bed use data
	
	# load reference CBGset
	CHIsetCombined <- read.csv("../GlCoSy/source/CHIsetCombined_09-15b.csv",header=TRUE,row.names=NULL)
	CHIsetCombinedDT<-data.table(CHIsetCombined)
	
	# add diabetes duration to the interest set
	admissionsWithDrugData_testSet$diabetesDuration<-(admissionsWithDrugData_testSet$dateplustime1-admissionsWithDrugData_testSet$diagnosisDateUnix)/(60*60*24*365.25)
		
	sigmaLimitsOutput<-as.data.frame(matrix(0,nrow=nrow(admissionsWithDrugData_testSet),ncol=8))
	colnames(sigmaLimitsOutput)<-c("ID","matchSetSize","sd1_bottom","sd1_top","sd2_bottom","sd2_top","sd3_bottom","sd3_top")
	
		for (r in seq(1,nrow(admissionsWithDrugData_testSet),1)) {
		
			if (r%%100==0) {print(r)}
		
				admissionOfInterest<-admissionsWithDrugData_testSet[r,]
			
				matchingSet<-subset(admissionsWithDrugData_testSet,(
					ID!=admissionOfInterest$ID &
					#age>(admissionOfInterest$age)-(ageWindow/2) &
					#age<(admissionOfInterest$age)+(ageWindow/2) &
					#
					medianGlu>(admissionOfInterest$medianGlu)-(medianWindow/2) &
					#
					medianGlu<(admissionOfInterest$medianGlu)+(medianWindow/2) &
					admissionDurationDays>(admissionOfInterest$admissionDurationDays)-(admissionDurationWindow/2) &
					admissionDurationDays<(admissionOfInterest$admissionDurationDays)+(admissionDurationWindow/2)
					#diabetesDuration>(admissionOfInterest$diabetesDuration)-(diabetesDurationWindow/2) &
					#diabetesDuration<(admissionOfInterest$diabetesDuration)+(diabetesDurationWindow/2)
					))
			
			sigmaLimitsOutput$ID[r]<-admissionOfInterest$ID
			sigmaLimitsOutput$matchSetSize[r]<-nrow(matchingSet)
			
			matchingCBGset<-merge(matchingSet,CHIsetCombined,by.x="ID",by.y="ID")
			
			if (nrow(matchingCBGset)>1) {
			
			sigmaLimitsOutput$sd1_bottom[r]<-exp(mean(log(matchingCBGset$yyyy.y)) - (1*(sd(log(matchingCBGset$yyyy.y)))))
			sigmaLimitsOutput$sd1_top[r]<-exp(mean(log(matchingCBGset$yyyy.y)) + (1*(sd(log(matchingCBGset$yyyy.y)))))
			sigmaLimitsOutput$sd2_bottom[r]<-exp(mean(log(matchingCBGset$yyyy.y)) - (2*(sd(log(matchingCBGset$yyyy.y)))))
			sigmaLimitsOutput$sd2_top[r]<-exp(mean(log(matchingCBGset$yyyy.y)) + (2*(sd(log(matchingCBGset$yyyy.y)))))
			sigmaLimitsOutput$sd3_bottom[r]<-exp(mean(log(matchingCBGset$yyyy.y)) - (3*(sd(log(matchingCBGset$yyyy.y)))))
			sigmaLimitsOutput$sd3_top[r]<-exp(mean(log(matchingCBGset$yyyy.y)) + (3*(sd(log(matchingCBGset$yyyy.y)))))
		}
			
		}
		
		# visualise limits
		plot(sigmaLimitsOutput$matchSetSize,sigmaLimitsOutput$sd1_bottom,pch=16,cex=0.3,col="red")
		points(sigmaLimitsOutput$matchSetSize,sigmaLimitsOutput$sd2_bottom,pch=16,cex=0.3,col="green")
		points(sigmaLimitsOutput$matchSetSize,sigmaLimitsOutput$sd3_bottom,pch=16,cex=0.3,col="blue")
		
		# save out file
		summaryOutputName <- paste("../GlCoSy/output/sigmaLimitsOutput_admissionsWithDrugData_testSet_medianAdmissionDuration.csv",sep="")
		write.table(sigmaLimitsOutput,file=summaryOutputName,sep=",",append=FALSE,col.names=TRUE)
		# read back in:
		sigmaLimitsOutput<-read.csv(summaryOutputName, header=TRUE , sep="," , row.names=NULL)
		
		default_sd1_bottom<-exp(mean(log(CHIsetCombined$yyyy)) - (1*(sd(log(CHIsetCombined$yyyy)))))
		default_sd1_top<-exp(mean(log(CHIsetCombined$yyyy)) + (1*(sd(log(CHIsetCombined$yyyy)))))
		default_sd2_bottom<-exp(mean(log(CHIsetCombined$yyyy)) - (2*(sd(log(CHIsetCombined$yyyy)))))
		default_sd2_top<-exp(mean(log(CHIsetCombined$yyyy)) + (2*(sd(log(CHIsetCombined$yyyy)))))
		default_sd3_bottom<-exp(mean(log(CHIsetCombined$yyyy)) - (3*(sd(log(CHIsetCombined$yyyy)))))
		default_sd3_top<-exp(mean(log(CHIsetCombined$yyyy)) + (3*(sd(log(CHIsetCombined$yyyy)))))
		
		ma <- function(x,n=18){filter(x,rep(1/n,n), sides=1)}

		
		for (rr in seq(1,nrow(sigmaLimitsOutput),1)) {
		
			if (rr%%100==0) {print(rr)}
				
				admissionData<-admissionsWithDrugData_testSet[rr,]
				admissionOfInterest<-sigmaLimitsOutput[rr,]
				CBGvalSubset<-CHIsetCombinedDT[(ID==admissionOfInterest$ID) & (dateplustime1>=admissionData$dateplustime1[1]) & dateplustime1<=(admissionData$dateplustime1[1]+admissionData$admissionDuration[1])]
				CBGvalSubset<-CBGvalSubset[order(CBGvalSubset$dateplustime1),]
				
				if (admissionOfInterest$matchSetSize<1) {
					admissionOfInterest$sd1_bottom<-default_sd1_bottom
					admissionOfInterest$sd1_top<-default_sd1_top
					admissionOfInterest$sd2_bottom<-default_sd2_bottom
					admissionOfInterest$sd2_top<-default_sd2_top
					admissionOfInterest$sd3_bottom<-default_sd3_bottom
					admissionOfInterest$sd3_top<-default_sd3_top
				}
				
				scoreSub<-data.frame(CBGvalSubset$yyyy)
				scoreSub$dateplustime1<-CBGvalSubset$dateplustime1
				scoreSub$score<-ifelse(CBGvalSubset$yyyy<admissionOfInterest$sd1_bottom | CBGvalSubset$yyyy>admissionOfInterest$sd1_top,-1,1)
				scoreSub$score<-ifelse(CBGvalSubset$yyyy<admissionOfInterest$sd2_bottom | CBGvalSubset$yyyy>admissionOfInterest$sd2_top,-2,scoreSub$score)

				scoreSub$cumScore<-cumsum(scoreSub$score)
				scoreSub$y2<-0; scoreSub$y2[1:nrow(scoreSub)-1]<-scoreSub$cumScore[2:nrow(scoreSub)]
				scoreSub$x2<-0; scoreSub$x2[1:nrow(scoreSub)-1]<-scoreSub$dateplustime1[2:nrow(scoreSub)]
				
				scoreSub$gradient<-(scoreSub$y2-scoreSub$cumScore)/(scoreSub$x2-scoreSub$dateplustime1)
				scoreSub$rollingGradient<-0; scoreSub$rollingGradient<-ma(scoreSub$gradient)
				

				
				#
				plot(CBGvalSubset$dateplustime1,CBGvalSubset$yyyy,ylim=c(0,28)); lines(CBGvalSubset$dateplustime1,CBGvalSubset$yyyy); abline(admissionOfInterest$sd1_bottom,0);abline(admissionOfInterest$sd1_top,0);abline(admissionOfInterest$sd2_bottom,0,col="blue");abline(admissionOfInterest$sd2_top,0,col="blue");abline(admissionOfInterest$sd3_bottom,0,col="red");abline(admissionOfInterest$sd3_top,0,col="red");abline(4,0,col="red",lty=3)
		 		par(new=TRUE)
		 		plot(CBGvalSubset$dateplustime1,scoreSub$cumScore,cex=0.01,axes=F); lines(CBGvalSubset$dateplustime1,scoreSub$cumScore)
				#axis(side=4)
				par(new=TRUE)
				plot(CBGvalSubset$dateplustime1,scoreSub$gradient,cex=0.2,col="red",axes=F); lines(CBGvalSubset$dateplustime1,scoreSub$gradient,col="red")
				par(new=TRUE)
				plot(CBGvalSubset$dateplustime1,scoreSub$rollingGradient,cex=0.2,col="blue",axes=F); lines(CBGvalSubset$dateplustime1,scoreSub$rollingGradient,col="blue")
				axis(side=4)
				
				# notes - need to manage situations where there are no limits set, as there are no matches. options are to exclude, or set general limits
		
			}
}

reportDrugProportions<-function(cases,controlPool) {
	casesL<-nrow(cases)
	controlPoolL<-nrow(controlPool)
	
	casesSU<-(sum(cases$preSU))/casesL
	casesIns<-(sum(cases$preIns))/casesL
	casesMF<-(sum(cases$preMF))/casesL
	casesOther<-(sum(cases$preOther))/casesL
	
	controlPoolSU<-(sum(controlPool$preSU))/controlPoolL
	controlPoolIns<-(sum(controlPool$preIns))/controlPoolL
	controlPoolMF<-(sum(controlPool$preMF))/controlPoolL
	controlPoolOther<-(sum(controlPool$preOther))/controlPoolL
	
	drugOutput<-as.data.frame(matrix(0,nrow=3,ncol=4))
	colnames(drugOutput)<-c("SU","Ins","MF","Other")
	rownames(drugOutput)<-c("cases","controls","pval")
	
	drugOutput$SU[1]<-casesSU
	drugOutput$SU[2]<-controlPoolSU
	drugOutput$SU[3]<-prop.test(c(sum(cases$preSU),sum(controlPool$preSU)),c(casesL,controlPoolL))$p.value
	
	drugOutput$Ins[1]<-casesIns
	drugOutput$Ins[2]<-controlPoolIns
	drugOutput$Ins[3]<-prop.test(c(sum(cases$preIns),sum(controlPool$preIns)),c(casesL,controlPoolL))$p.value
	
	drugOutput$MF[1]<-casesMF
	drugOutput$MF[2]<-controlPoolMF
	drugOutput$MF[3]<-prop.test(c(sum(cases$preMF),sum(controlPool$preMF)),c(casesL,controlPoolL))$p.value
	
	drugOutput$Other[1]<-casesOther
	drugOutput$Other[2]<-controlPoolOther
	drugOutput$Other[3]<-prop.test(c(sum(cases$preOther),sum(controlPool$preOther)),c(casesL,controlPoolL))$p.value

}


simpleSurvivalPlot_SUcaseVScontrol<-function(inputFrame,ylimMin,postDischargeStartDay) {
  
  SurvivalData<-inputFrame
  
  DaySeconds<-(60*60*24)
  shortCensorPeriodStartDay  <- DaySeconds*postDischargeStartDay
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
  
  SurvivalData$diabetesDurationYears <- (SurvivalData$dateplustime1 - SurvivalData$diagnosisDateUnix) / (DaySeconds * 365.25)
  
  # DaySeconds<-(60*60*24)
  # shortCensorPeriodStartDay  <- DaySeconds*postDischargeStartDay
  # shortCensorPeriodEndDay    <- DaySeconds*10000
  # 
  # lastDOD<-endDateUnix
  # SurvivalData$dateOfDischarge<-sampleDateUnix
  # SurvivalData$timeToDeath<-ifelse(SurvivalData$isDead==1,(SurvivalData$DeathDateUnix-SurvivalData$dateOfDischarge),0)
  # #		SurvivalData$timeToDeath<-SurvivalData$timeToDeath/DaySeconds
  # SurvivalData$timeToDeathInterval<-ifelse(SurvivalData$isDead==0,(lastDOD-SurvivalData$dateOfDischarge),SurvivalData$timeToDeath)
  # SurvivalData$timeToDeathInterval[is.na(SurvivalData$timeToDeathInterval)]<-0; SurvivalData<-subset(SurvivalData,timeToDeathInterval>0)
  # # SurvivalData$timeToDeathInterval<-SurvivalData$timeToDeathInterval/(60*60*24*365.25)
  # 
  # SurvivalData$shortDeathEvent <- SurvivalData$isDead
  # # SurvivalData$shortDeathEvent <- ifelse(SurvivalData$isDead==1 & SurvivalData$timeToDeath>=(shortCensorPeriodStartDay) & SurvivalData$timeToDeath<(shortCensorPeriodEndDay),1,0)	
  # 
  # #  SurvivalData$sexDigit<-ifelse(nchar(SurvivalData$charID==9),as.numeric(substr(SurvivalData$charID,8,8)),as.numeric(substr(SurvivalData$charID,9,9)))
  # # SurvivalData$sexNumber<-ifelse(SurvivalData$sexDigit%%2==0,1,0)
  # #  SurvivalData$sex<-factor(1*(SurvivalData$sexNumber <1),levels=0:1,labels=c("F","M"))
  
  
  # mfitAge50<-survfit(Surv(timeToDeathInterval, shortDeathEvent) ~ (caseFlag == 1), data = SurvivalData)
  mfitAge50<-survfit(Surv(timeToDeathInterval, shortDeathEvent) ~ (case_control == 1), data = SurvivalData)
  shortPlotTitle <- paste("Mortality, time ",round(shortCensorPeriodStartDay)/DaySeconds," to ",round(max(SurvivalData$timeToDeathInterval))/DaySeconds," days\n n= ",nrow(SurvivalData),", threshold: ",quantile(SurvivalData$hba1cIQRinRange)[3],sep="")
  plot(mfitAge50,mark.time=T,lty=1:6,conf.int=F,col=c("black","red","blue","green","orange","purple"),main=shortPlotTitle,xlim=c(shortCensorPeriodStartDay,round(max(SurvivalData$timeToDeathInterval))),lwd=3,ylim=c(ylimMin,1))
 # mfitAge50.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ (caseFlag == 1), data = SurvivalData)
  
  mfitAge50.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ (case_control == 1) + age + medianGlu + admissionDurationDays + diabetesDurationYears + sexNumber, data = SurvivalData)
  
  pVal <- summary(mfitAge50.coxph)$coef[,5]; HR <- round(exp(coef(mfitAge50.coxph)),2)
  legendText <- paste("p = ",pVal," | HR = ",HR,sep="")
  summarySurvfit <- summary(mfitAge50); legendNames <- row.names(summarySurvfit$table)
  legend("bottomleft",c(legendNames),lty=1:6,col=c("black","red","blue","green","orange","purple"),cex=0.8); legend("topright",legendText,cex=0.6)
  
  print(mfitAge50.coxph)
  
}

simpleSurvivalPlot<-function(inputFrame,postDischargeStartDay) {
	
	SurvivalData<-inputFrame
	
	DaySeconds<-(60*60*24)
	shortCensorPeriodStartDay  <- DaySeconds*postDischargeStartDay
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
	
	
		mfitAge50<-survfit(Surv(timeToDeathInterval, shortDeathEvent) ~ (IQR>=quantile(SurvivalData$IQR)[3]), data = SurvivalData)
		shortPlotTitle <- paste("Mortality, time ",round(shortCensorPeriodStartDay)/DaySeconds," to ",round(max(SurvivalData$timeToDeathInterval))/DaySeconds,"days\nTest: ",", threshold: ",quantile(SurvivalData$IQR)[3],sep="")
		plot(mfitAge50,mark.time=F,lty=1:6,conf.int=F,col=c("black","red","blue","green","orange","purple"),main=shortPlotTitle,xlim=c(shortCensorPeriodStartDay,round(max(SurvivalData$timeToDeathInterval))),lwd=3)
		mfitAge50.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ age+admissionDuration+(IQR>=quantile(SurvivalData$IQR)[3]), data = SurvivalData)
			pVal <- summary(mfitAge50.coxph)$coef[,5]; HR <- round(exp(coef(mfitAge50.coxph)),2)
			legendText <- paste("p = ",pVal," | HR = ",HR,sep="")
		summarySurvfit <- summary(mfitAge50); legendNames <- row.names(summarySurvfit$table)
		legend("bottomleft",c(legendNames),lty=1:6,col=c("black","red","blue","green","orange","purple"),cex=0.8); legend("topright",legendText,cex=0.6)

		mfitAge50<-survfit(Surv(timeToDeathInterval, shortDeathEvent) ~ (IQR>=quantile(SurvivalData$IQR)[4]), data = SurvivalData)
		shortPlotTitle <- paste("Mortality, time ",round(shortCensorPeriodStartDay)/DaySeconds," to ",round(max(SurvivalData$timeToDeathInterval))/DaySeconds,"days\nTest: ",", threshold: ",quantile(SurvivalData$IQR)[4],sep="")
		plot(mfitAge50,mark.time=F,lty=1:6,conf.int=F,col=c("black","red","blue","green","orange","purple"),main=shortPlotTitle,xlim=c(shortCensorPeriodStartDay,round(max(SurvivalData$timeToDeathInterval))),lwd=3)
		mfitAge50.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ age+admissionDuration+(IQR>=quantile(SurvivalData$IQR)[4]), data = SurvivalData)
			pVal <- summary(mfitAge50.coxph)$coef[,5]; HR <- round(exp(coef(mfitAge50.coxph)),2)
			legendText <- paste("p = ",pVal," | HR = ",HR,sep="")
		summarySurvfit <- summary(mfitAge50); legendNames <- row.names(summarySurvfit$table)
		legend("bottomleft",c(legendNames),lty=1:6,col=c("black","red","blue","green","orange","purple"),cex=0.8); legend("topright",legendText,cex=0.6)

}

addSurvivalData<-function(inputFrame,plotTitle,test,threshold,postDischargeStartDay) {
	
	SurvivalData<-inputFrame
	
	DaySeconds<-(60*60*24)
	shortCensorPeriodStartDay  <- DaySeconds*postDischargeStartDay
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
# match switch controls matching of median/IQR - 1=median only, 2=both, 3=IQR only
controlSelection<-function(cases,controlPool,plotTitleInsert,matchSwitch,postDischargeStartDay) {
	# cases<-subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & (preSU==0 & preIns==1 & preMF==1))
	# controlPool<-subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" & (preSU==1 & preIns==0 & preMF==1))
	
	# admissionsWithDrugDataT2<-subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus")
	# cases<-subset(admissionsWithDrugDataT2, IQR>=4.5)
	# controlPool<-subset(admissionsWithDrugDataT2, IQR<4.5)
	
	# remove those cases with a diagnosis date of 1/1/1900 - presumably a default within SCI diabetes
	
	# hypoThreshold set for hypo matching. default 4
	hypoThreshold<-4
	
	cases<-cases[cases$diagnosisDateUnix!=(-2208988800),];controlPool<-controlPool[controlPool$diagnosisDateUnix!=(-2208988800),]
	
	# calculate case/control diabetes duration
	cases$diabetesDuration<-(cases$dateplustime1-cases$diagnosisDateUnix)/(60*60*24*365.25)
	controlPool$diabetesDuration<-(controlPool$dateplustime1-controlPool$diagnosisDateUnix)/(60*60*24*365.25)
	controlPool$hypoForMatching<-ifelse(controlPool$minGlu<hypoThreshold,1,0)
	
	controlPool$used<-0
	
	print(dim(cases));print(dim(controlPool))
	
	controlPoolDT<-data.table(controlPool)
	
	# set up flags for survival plotting
	admissionsWithDrugData$caseFlag<-0
	admissionsWithDrugData$controlFlag<-0
	
	testOutput<-as.data.frame(matrix(0,nrow=nrow(cases),ncol=17))
	colnames(testOutput)<-c("caseAge","controlAge","caseMedian","controlMedian","caseIQR","controlIQR","caseMort","controlMort","caseHypoEpisodes","controlHypoEpisodes","caseAdmissionDurationDays","controlAdmissionDurationDays","numberOfMatches","caseDiabetesDuration","controlDiabetesDuration","caseeGFR","controleGFR")
	
	# match median only (original)
	if (matchSwitch==1) {	
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
	}

	# match median _and_ IQR
	if (matchSwitch==2) {	
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
			caseIQR<-cases$IQR[i]
		
			# 
			matchingSUB<-availableControlPool[which(
			availableControlPool$age>(caseAge-(ageWindow/2)) & 
			availableControlPool$age<(caseAge+(ageWindow/2)) & 
			availableControlPool$medianGlu>(caseMedian-(medianWindow/2)) & 
			availableControlPool$medianGlu<(caseMedian+(medianWindow/2)) & 
			round(availableControlPool$admissionDurationDays,1)>(caseDuration-(admissionDurationWindow/2)) & 
			round(availableControlPool$admissionDurationDays,1)<(caseDuration+(admissionDurationWindow/2)) &
			availableControlPool$diabetesDuration>(caseDiabetesDuration-(diabetesDurationWindow/2)) &
			availableControlPool$diabetesDuration<(caseDiabetesDuration+(diabetesDurationWindow/2)) &
			availableControlPool$IQR>(caseIQR-(IQRwindow/2)) & 
			availableControlPool$IQR<(caseIQR+(IQRwindow/2))),]
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
	}

	# match IQR only
	if (matchSwitch==3) {	
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
			caseIQR<-cases$IQR[i]
		
			# 
			matchingSUB<-availableControlPool[which(
			availableControlPool$age>(caseAge-(ageWindow/2)) & 
			availableControlPool$age<(caseAge+(ageWindow/2)) & 
			round(availableControlPool$admissionDurationDays,1)>(caseDuration-(admissionDurationWindow/2)) & 
			round(availableControlPool$admissionDurationDays,1)<(caseDuration+(admissionDurationWindow/2)) &
			availableControlPool$diabetesDuration>(caseDiabetesDuration-(diabetesDurationWindow/2)) &
			availableControlPool$diabetesDuration<(caseDiabetesDuration+(diabetesDurationWindow/2)) &
			availableControlPool$IQR>(caseIQR-(IQRwindow/2)) & 
			availableControlPool$IQR<(caseIQR+(IQRwindow/2))),]
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
	}
	
	# match median and hypoepisodes
	if (matchSwitch==4) {	
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
			caseHypoEpisodes<-cases$ID_ADMISSIONhypoEpisodes4.60[i]
		
			# 
			matchingSUB<-availableControlPool[which(
			availableControlPool$age>(caseAge-(ageWindow/2)) & 
			availableControlPool$age<(caseAge+(ageWindow/2)) & 
			availableControlPool$medianGlu>(caseMedian-(medianWindow/2)) & 
			availableControlPool$medianGlu<(caseMedian+(medianWindow/2)) & 
			round(availableControlPool$admissionDurationDays,1)>(caseDuration-(admissionDurationWindow/2)) & 
			round(availableControlPool$admissionDurationDays,1)<(caseDuration+(admissionDurationWindow/2)) &
			availableControlPool$diabetesDuration>(caseDiabetesDuration-(diabetesDurationWindow/2)) &
			availableControlPool$diabetesDuration<(caseDiabetesDuration+(diabetesDurationWindow/2)) &
			availableControlPool$ID_ADMISSIONhypoEpisodes4.60==caseHypoEpisodes
			),]
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
	}
	
	# match IQR and hypoepisodes
	if (matchSwitch==5) {	
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
			caseIQR<-cases$IQR[i]
			caseHypoEpisodes<-cases$ID_ADMISSIONhypoEpisodes4.60[i]
		
			# 
			matchingSUB<-availableControlPool[which(
			availableControlPool$age>(caseAge-(ageWindow/2)) & 
			availableControlPool$age<(caseAge+(ageWindow/2)) & 
			availableControlPool$IQR>(caseIQR-(IQRwindow/2)) & 
			availableControlPool$IQR<(caseIQR+(IQRwindow/2)) &
			round(availableControlPool$admissionDurationDays,1)>(caseDuration-(admissionDurationWindow/2)) & 
			round(availableControlPool$admissionDurationDays,1)<(caseDuration+(admissionDurationWindow/2)) &
			availableControlPool$diabetesDuration>(caseDiabetesDuration-(diabetesDurationWindow/2)) &
			availableControlPool$diabetesDuration<(caseDiabetesDuration+(diabetesDurationWindow/2)) &
			availableControlPool$ID_ADMISSIONhypoEpisodes4.60==caseHypoEpisodes
			),]
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
	}
	
	# match median and hypoepisodes - as presence or absence of hypoglycemia to avoid over matching for hypos
	if (matchSwitch==6) {	
		#	variable of interest in this loop IQR
				
		for (i in seq(1,nrow(cases),1)) {
			if (i%%100==0) {print(i)}
		
			availableControlPool<-controlPool[which(controlPool$used==0),]
			testOutput$caseAge[i]<-cases$age[i]
			testOutput$caseMedian[i]<-cases$medianGlu[i]
			testOutput$caseIQR[i]<-cases$IQR[i]
			testOutput$caseMort[i]<-cases$deathEvent[i]
			testOutput$caseHypoEpisodes[i]<-ifelse(cases$minGlu[i]<hypoThreshold,1,0)
			testOutput$caseAdmissionDurationDays[i]<-cases$admissionDurationDays[i]
			testOutput$caseDiabetesDuration[i]<-cases$diabetesDuration[i]
		
			# flag as a case in admissionsWithDrugs for Survival Plotting
			# admissionsWithDrugData$caseFlag<-ifelse(cases$ID[i]==admissionsWithDrugData$ID,1,admissionsWithDrugData$caseFlag)
		
			caseAge<-cases$age[i]
			caseMedian<-cases$medianGlu[i]
			caseDuration<-round(cases$admissionDurationDays[i],1)
			case_nCBG<-cases$nCBGperAdmission[i]
			caseDiabetesDuration<-cases$diabetesDuration[i]
			caseHypoEpisodes<-ifelse(cases$minGlu[i]<hypoThreshold,1,0)
		
			# 
			matchingSUB<-availableControlPool[which(
			availableControlPool$age>(caseAge-(ageWindow/2)) & 
			availableControlPool$age<(caseAge+(ageWindow/2)) & 
			availableControlPool$medianGlu>(caseMedian-(medianWindow/2)) & 
			availableControlPool$medianGlu<(caseMedian+(medianWindow/2)) & 
			round(availableControlPool$admissionDurationDays,1)>(caseDuration-(admissionDurationWindow/2)) & 
			round(availableControlPool$admissionDurationDays,1)<(caseDuration+(admissionDurationWindow/2)) &
			availableControlPool$diabetesDuration>(caseDiabetesDuration-(diabetesDurationWindow/2)) &
			availableControlPool$diabetesDuration<(caseDiabetesDuration+(diabetesDurationWindow/2)) &
			availableControlPool$hypoForMatching==caseHypoEpisodes
			),]
			# matchingSUB<-availableControlPool[which(availableControlPool$age>(caseAge-2.5) & availableControlPool$age<(caseAge+2.5) & availableControlPool$medianGlu>(caseMedian-0.1) & availableControlPool$medianGlu<(caseMedian+0.1) & nCBGperAdmission>(case_nCBG-3) & nCBGperAdmission<(case_nCBG+3)),]
		
			if (nrow(matchingSUB)>0) {
				randomRow<-sample(1:nrow(matchingSUB),1)
				testOutput$numberOfMatches[i]<-nrow(matchingSUB)
				testOutput$controlAge[i]<-matchingSUB$age[randomRow]
				testOutput$controlMedian[i]<-matchingSUB$medianGlu[randomRow]
				testOutput$controlIQR[i]<-matchingSUB$IQR[randomRow] # take 1st row ? need to randomise row taken
				testOutput$controlMort[i]<-matchingSUB$deathEvent[randomRow]
				testOutput$controlHypoEpisodes[i]<-matchingSUB$hypoForMatching[randomRow]
				testOutput$controlAdmissionDurationDays[i]<-matchingSUB$admissionDurationDays[randomRow]
				testOutput$controlDiabetesDuration[i]<-matchingSUB$diabetesDuration[randomRow]
			
				controlPool$used<-ifelse(controlPool$ID==matchingSUB$ID[randomRow],1,controlPool$used)
			
				admissionsWithDrugData$controlFlag<-ifelse(matchingSUB$ID[randomRow]==admissionsWithDrugData$ID,i,admissionsWithDrugData$controlFlag)
				admissionsWithDrugData$caseFlag<-ifelse(cases$ID[i]==admissionsWithDrugData$ID,i,admissionsWithDrugData$caseFlag)
			
			}
			
		}
	}
	
	# cases/controls - with or without hypoglycaemia - match other parameters
	if (matchSwitch==7) {	
		#	variable of interest in this loop admissionDurationDays
				
		for (i in seq(1,nrow(cases),1)) {
			if (i%%100==0) {print(i)}
		
			availableControlPool<-controlPool[which(controlPool$used==0),]
			testOutput$caseAge[i]<-cases$age[i]
			testOutput$caseMedian[i]<-cases$medianGlu[i]
			testOutput$caseIQR[i]<-cases$IQR[i]
			testOutput$caseMort[i]<-cases$deathEvent[i]
			testOutput$caseHypoEpisodes[i]<-ifelse(cases$minGlu[i]<hypoThreshold,1,0)
			testOutput$caseAdmissionDurationDays[i]<-cases$admissionDurationDays[i]
			testOutput$caseDiabetesDuration[i]<-cases$diabetesDuration[i]
		
			# flag as a case in admissionsWithDrugs for Survival Plotting
			# admissionsWithDrugData$caseFlag<-ifelse(cases$ID[i]==admissionsWithDrugData$ID,1,admissionsWithDrugData$caseFlag)
		
			caseAge<-cases$age[i]
			caseMedian<-cases$medianGlu[i]
			caseDuration<-round(cases$admissionDurationDays[i],1)
			case_nCBG<-cases$nCBGperAdmission[i]
			caseDiabetesDuration<-cases$diabetesDuration[i]
			caseHypoEpisodes<-ifelse(cases$minGlu[i]<hypoThreshold,1,0)
			caseIQR<-cases$IQR[i]
		
			# 
			matchingSUB<-availableControlPool[which(
				availableControlPool$age>(caseAge-(ageWindow/2)) & 
				availableControlPool$age<(caseAge+(ageWindow/2)) & 
				availableControlPool$medianGlu>(caseMedian-(medianWindow/2)) & 
				availableControlPool$medianGlu<(caseMedian+(medianWindow/2)) & 
				#round(availableControlPool$admissionDurationDays,1)>(caseDuration-(admissionDurationWindow/2)) & 
				#round(availableControlPool$admissionDurationDays,1)<(caseDuration+(admissionDurationWindow/2)) &
				#availableControlPool$diabetesDuration>(caseDiabetesDuration-(diabetesDurationWindow/2)) &
				#availableControlPool$diabetesDuration<(caseDiabetesDuration+(diabetesDurationWindow/2)) &
				#availableControlPool$hypoForMatching==caseHypoEpisodes
				availableControlPool$IQR>(caseIQR-(IQRwindow/2)) & 
				availableControlPool$IQR<(caseIQR+(IQRwindow/2))
			),]
			# matchingSUB<-availableControlPool[which(availableControlPool$age>(caseAge-2.5) & availableControlPool$age<(caseAge+2.5) & availableControlPool$medianGlu>(caseMedian-0.1) & availableControlPool$medianGlu<(caseMedian+0.1) & nCBGperAdmission>(case_nCBG-3) & nCBGperAdmission<(case_nCBG+3)),]
		
			if (nrow(matchingSUB)>0) {
				randomRow<-sample(1:nrow(matchingSUB),1)
				testOutput$numberOfMatches[i]<-nrow(matchingSUB)
				testOutput$controlAge[i]<-matchingSUB$age[randomRow]
				testOutput$controlMedian[i]<-matchingSUB$medianGlu[randomRow]
				testOutput$controlIQR[i]<-matchingSUB$IQR[randomRow] # take 1st row ? need to randomise row taken
				testOutput$controlMort[i]<-matchingSUB$deathEvent[randomRow]
				testOutput$controlHypoEpisodes[i]<-matchingSUB$hypoForMatching[randomRow]
				testOutput$controlAdmissionDurationDays[i]<-matchingSUB$admissionDurationDays[randomRow]
				testOutput$controlDiabetesDuration[i]<-matchingSUB$diabetesDuration[randomRow]
			
				controlPool$used<-ifelse(controlPool$ID==matchingSUB$ID[randomRow],1,controlPool$used)
			
				admissionsWithDrugData$controlFlag<-ifelse(matchingSUB$ID[randomRow]==admissionsWithDrugData$ID,i,admissionsWithDrugData$controlFlag)
				admissionsWithDrugData$caseFlag<-ifelse(cases$ID[i]==admissionsWithDrugData$ID,i,admissionsWithDrugData$caseFlag)
			
			}
			
		}
	}
	
	# match median _and_ IQR _and_ eGFR
	if (matchSwitch==8) {	
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
			testOutput$caseeGFR[i]<-cases$eGFR[i]
		
			# flag as a case in admissionsWithDrugs for Survival Plotting
			# admissionsWithDrugData$caseFlag<-ifelse(cases$ID[i]==admissionsWithDrugData$ID,1,admissionsWithDrugData$caseFlag)
		
			caseAge<-cases$age[i]
			caseMedian<-cases$medianGlu[i]
			caseDuration<-round(cases$admissionDurationDays[i],1)
			case_nCBG<-cases$nCBGperAdmission[i]
			caseDiabetesDuration<-cases$diabetesDuration[i]
			caseIQR<-cases$IQR[i]
			caseeGFR<-cases$eGFR[i]
		
			# 
			matchingSUB<-availableControlPool[which(
			availableControlPool$age>(caseAge-(ageWindow/2)) & 
			availableControlPool$age<(caseAge+(ageWindow/2)) & 
			availableControlPool$medianGlu>(caseMedian-(medianWindow/2)) & 
			availableControlPool$medianGlu<(caseMedian+(medianWindow/2)) & 
			round(availableControlPool$admissionDurationDays,1)>(caseDuration-(admissionDurationWindow/2)) & 
			round(availableControlPool$admissionDurationDays,1)<(caseDuration+(admissionDurationWindow/2)) &
			availableControlPool$diabetesDuration>(caseDiabetesDuration-(diabetesDurationWindow/2)) &
			availableControlPool$diabetesDuration<(caseDiabetesDuration+(diabetesDurationWindow/2)) &
			availableControlPool$eGFR>(caseeGFR-(eGFRwindow/2)) &
			availableControlPool$eGFR<(caseeGFR+(eGFRwindow/2)) &
			availableControlPool$IQR>(caseIQR-(IQRwindow/2)) & 
			availableControlPool$IQR<(caseIQR+(IQRwindow/2))),]
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
				testOutput$controleGFR[i]<-matchingSUB$eGFR[randomRow]
			
				controlPool$used<-ifelse(controlPool$ID==matchingSUB$ID[randomRow],1,controlPool$used)
			
				admissionsWithDrugData$controlFlag<-ifelse(matchingSUB$ID[randomRow]==admissionsWithDrugData$ID,i,admissionsWithDrugData$controlFlag)
				admissionsWithDrugData$caseFlag<-ifelse(cases$ID[i]==admissionsWithDrugData$ID,i,admissionsWithDrugData$caseFlag)
			
			}
			
		}
	}
	
		# write out data for later plotting/analysis
	fileName<-paste(plotTitleInsert,"SU",cases$preSU[1],"-Ins",cases$preIns[1],"-MF",cases$preMF[1],"_vs_SU",controlPool$preSU[1],"-Ins",controlPool$preIns[1],"-MF",controlPool$preMF[1],"_",minCBGperAdmission+1,"CBGmin","_ageWindow-",ageWindow,"_medianWindow-",medianWindow,"_admissionDurationWindow-",admissionDurationWindow,"_diabetesDurationWindow-",diabetesDurationWindow,sep="")	
#		fileName<-paste(plotTitleInsert,"PRE-includingNameSearch-_SU",cases$preSU[1],"-Ins",cases$preIns[1],"-MF",cases$preMF[1],"_vs_SU",controlPool$preSU[1],"-Ins",controlPool$preIns[1],"-MF",controlPool$preMF[1],"_",minCBGperAdmission+1,"CBGmin","_ageWindow-",ageWindow,"_medianWindow-",medianWindow,"_admissionDurationWindow-",admissionDurationWindow,"_diabetesDurationWindow-",diabetesDurationWindow,sep="")	
		# fileName<-paste("SU",cases$spanSU[1],"-Ins",cases$spanIns[1],"-MF",cases$spanMF[1],"_vs_SU",controlPool$spanSU[1],"-Ins",controlPool$spanIns[1],"-MF",controlPool$spanMF[1],"_",minCBGperAdmission+1,"CBGmin","_ageWindow-",ageWindow,"_medianWindow-",medianWindow,"_admissionDurationWindow-",admissionDurationWindow,sep="")
		testOutputwriteName <- paste("~/R/GlCoSy/dataForSubmissions/drugClassPaper/textOutput/",fileName,".csv",sep="")
		write.table(testOutput,file=testOutputwriteName,sep=",",append=FALSE,col.names=TRUE)
		#admissionsWithDrugDatawriteName <- paste("../GlCoSy/GGCdrug/admissionsWithDrugData",fileName,"-2CBGmin-0.4median.csv",sep="")
		#write.table(admissionsWithDrugData,file=admissionsWithDrugDatawriteName,sep=",",append=FALSE,col.names=TRUE)
	
	
	# testOutput<-read.csv(testOutputwriteName, header=TRUE , sep="," , row.names=NULL)
	# admissionsWithDrugData<-read.csv(admissionsWithDrugDatawriteName, header=TRUE , sep="," , row.names=NULL)
	

# write analysis results and quantiles to file
textPasteName<-paste(fileName,".txt",sep="")
textfilename <- paste("~/R/GlCoSy/dataForSubmissions/drugClassPaper/textOutput/",textPasteName,sep="")
sink(textfilename)

	testOutputCaseControl<-subset(testOutput,controlIQR>0 & caseIQR>0)
	print("number of pairs (with more than one CBG): "); print(nrow(testOutputCaseControl))
	
	# IQR test
	IQRttest<-t.test(log(testOutputCaseControl$caseIQR),log(testOutputCaseControl$controlIQR),paired=F)
	print("IQR ttest:");print(IQRttest)
	
	# IQR test
	IQRttest_wilcox<-wilcox.test(testOutputCaseControl$caseIQR,testOutputCaseControl$controlIQR,paired=F)
	print("IQR wilcox:");print(IQRttest_wilcox)
	
	
	print("IQR quantiles: case, control")	
	print(quantile(testOutputCaseControl$caseIQR)); print(quantile(testOutputCaseControl$controlIQR))
	
	# mortality test
	mortalityProp<-prop.test(c(sum(testOutputCaseControl$caseMort),sum(testOutputCaseControl$controlMort)),c(nrow(testOutputCaseControl),nrow(testOutputCaseControl)))
	print("chi-square test overall proportional mortality:"); print(mortalityProp)
	
	# hypo episode event rate
	if (sum(testOutputCaseControl$caseHypoEpisodes) < sum(testOutputCaseControl$controlHypoEpisodes)) { hypoProp<-prop.test(c(sum(testOutputCaseControl$caseHypoEpisodes),sum(testOutputCaseControl$controlHypoEpisodes)),c(nrow(testOutputCaseControl),nrow(testOutputCaseControl)))
	print("chi-square test proportional hypo episodes per admission"); print(hypoProp) }
		
	print("cases: hypo episodes/day")
	print(sum(testOutputCaseControl$caseHypoEpisodes)/sum(testOutputCaseControl$caseAdmissionDurationDays))
	print("controls: hypo episodes/day")
	print(sum(testOutputCaseControl$controlHypoEpisodes)/sum(testOutputCaseControl$controlAdmissionDurationDays))
	
	print('prop test. hypo ep per day')
	print(prop.test(c(sum(testOutputCaseControl$caseHypoEpisodes), sum(testOutputCaseControl$controlHypoEpisodes)), c(sum(testOutputCaseControl$caseAdmissionDurationDays), sum(testOutputCaseControl$controlAdmissionDurationDays))))
	
	ttestAdmissionDuration<-t.test(log(testOutputCaseControl$caseAdmissionDurationDays),log(testOutputCaseControl$controlAdmissionDurationDays))
	print("admissionDuration ttest:");print(ttestAdmissionDuration)
	
	ttestAdmissionDuration_wilcox<-wilcox.test(testOutputCaseControl$caseAdmissionDurationDays,testOutputCaseControl$controlAdmissionDurationDays)
	print("admissionDuration wilcox:");print(ttestAdmissionDuration_wilcox)
	
	print("admissionDuration quantiles: case, control")
	print(quantile(testOutputCaseControl$caseAdmissionDurationDays)); print(quantile(testOutputCaseControl$controlAdmissionDurationDays))
	
	ttestAge<-t.test(testOutputCaseControl$caseAge,testOutputCaseControl$controlAge,paired=F)
	print("age ttest:");print(ttestAge)
	
	ttestAge_wilcox<-wilcox.test(testOutputCaseControl$caseAge,testOutputCaseControl$controlAge,paired=F)
	print("age wilcox:");print(ttestAge_wilcox)
	
	print("age quantiles: case, control")
	print(quantile(testOutputCaseControl$caseAge));print(quantile(testOutputCaseControl$controlAge))
		
	ttestMedian<-t.test(log(testOutputCaseControl$caseMedian),log(testOutputCaseControl$controlMedian),paired=F)
	print("median ttest:");print(ttestMedian)
	print("median quantiles: case, control")
	print(quantile(testOutputCaseControl$caseMedian));print(quantile(testOutputCaseControl$controlMedian))
	
	ttestDiabetesDuration<-t.test(log(subset(testOutputCaseControl,caseDiabetesDuration>0 & controlDiabetesDuration>0)$caseDiabetesDuration),log(subset(testOutputCaseControl,caseDiabetesDuration>0 & controlDiabetesDuration>0)$controlDiabetesDuration),paired=F)
	print("diabetes duration ttest:");print(ttestDiabetesDuration)
	
	ttestDiabetesDuration_wilcox<-wilcox.test(subset(testOutputCaseControl,caseDiabetesDuration>0 & controlDiabetesDuration>0)$caseDiabetesDuration,subset(testOutputCaseControl,caseDiabetesDuration>0 & controlDiabetesDuration>0)$controlDiabetesDuration,paired=F)
	print("diabetes duration wilcox:");print(ttestDiabetesDuration_wilcox)
	
	
	print("diabetes duration quantiles: case, control")
	print(quantile(testOutputCaseControl$caseDiabetesDuration));print(quantile(testOutputCaseControl$controlDiabetesDuration))
	
sink()
		
	
	# mark pairs matched in order that the cluster argument can be used to perform a matched survival analysis (included in the cox prop hazard p value)
	survivalSUB<-subset(admissionsWithDrugData,caseFlag>0 | controlFlag>0)
	survivalSUB$pairNumber<-ifelse(survivalSUB$caseFlag>0,survivalSUB$caseFlag,ifelse(survivalSUB$controlFlag>0,survivalSUB$controlFlag,0))
	survivalSUB$caseFlag<-ifelse(survivalSUB$caseFlag>0,1,0)
	survivalSUB$controlFlag<-ifelse(survivalSUB$controlFlag>0,1,0)
	
	plotPasteName<-paste(fileName,".pdf",sep="")
	plotfilename <- paste("~/R/GlCoSy/dataForSubmissions/drugClassPaper/textOutput/",plotPasteName,sep="")
	modifiedFilename<-paste(plotTitleInsert,"\nPRE-includingNameSearch-_SU",cases$preSU[1],"-Ins",cases$preIns[1],"-MF",cases$preMF[1],"_vs_SU",controlPool$preSU[1],"-Ins",controlPool$preIns[1],"-MF",controlPool$preMF[1],"_",minCBGperAdmission+1,"\nCBGmin","_ageWindow-",ageWindow,"_medianWindow-",medianWindow,"_admissionDurationWindow-",admissionDurationWindow,"_diabetesDurationWindow-",diabetesDurationWindow,sep="")	
		
	pdf(plotfilename, width=16, height=9)
		plotName<-paste(fileName,"",sep="")
		addSurvivalData(survivalSUB,modifiedFilename,"Drug",0,postDischargeStartDay)
	dev.off()

}

regression<-function(inputFrame,deadAtNyears,divisionsOfIQR) {
	#regressionFrame<-admissionsWithDrugData_testSet
	 #
	 regressionFrame<-subset(admissionsWithDrugData_testSet, preIns==1)
	#
	 deadAtNyears<-3
	 #
	  divisionsOfIQR<-5
	 
	 
	 
	lastDOD<-max(regressionFrame$deathDateUnix)
	regressionFrame$dateOfDischarge<-regressionFrame$dateplustime1+regressionFrame$admissionDuration
	regressionFrame$timeToDeath<-ifelse(regressionFrame$deathEvent==1,(regressionFrame$deathDateUnix-regressionFrame$dateOfDischarge),0)
	
	regressionFrame$deadAt1year<-ifelse(regressionFrame$timeToDeath>0 & regressionFrame$timeToDeath<(deadAtNyears*365.25*60*60*24),1,0)
	regressionFrame$deadAt1yearFactor<-as.factor(regressionFrame$deadAt1year)
	regressionFrame$minGlu4<-ifelse(regressionFrame$minGlu<4,1,0)
	regressionFrame$minGlu4Factor<-as.factor(regressionFrame$minGlu4)
		regressionFrame$minGlu3p9<-ifelse(regressionFrame$minGlu<3.9,1,0)
		regressionFrame$minGlu3p9Factor<-as.factor(regressionFrame$minGlu3p9)
			regressionFrame$minGlu3<-ifelse(regressionFrame$minGlu<3,1,0)
			regressionFrame$minGlu3Factor<-as.factor(regressionFrame$minGlu3)
	
	regressionFrame$adequateFollowup<-ifelse((lastDOD-regressionFrame$dateOfDischarge)>=(deadAtNyears*365.25*60*60*24),1,0)
	regressionFrame$diabetesDuration<-(regressionFrame$dateplustime1-regressionFrame$diagnosisDateUnix)/(60*60*24*365.25)
	
	regressionFrame<-subset(regressionFrame,adequateFollowup==1)
	#regressionTestSet<-data.frame(regressionFrame$deadAt1yearFactor,regressionFrame$age,regressionFrame$diabetesDuration,regressionFrame$IQR,regressionFrame$medianGlu,regressionFrame$admissionDurationDays,regressionFrame$minGlu4Factor)
	#colnames(regressionTestSet)<-c("deadAt1yearFactor","age","diabetesDuration","IQR","medianGlu","admissionDurationDays","minGlu4")
	
	regressionTestSet<-data.frame(regressionFrame$deadAt1yearFactor,regressionFrame$age,regressionFrame$diabetesDuration,regressionFrame$IQR,regressionFrame$medianGlu,regressionFrame$admissionDurationDays,regressionFrame$preSU,regressionFrame$preIns,regressionFrame$preMF,regressionFrame$preOther)
	colnames(regressionTestSet)<-c("deadAt1yearFactor","age","diabetesDuration","IQR","medianGlu","admissionDurationDays","preSU","preIns","preMF","preOther")
	
	#regressionTestSet<-data.frame(regressionFrame$deadAt1yearFactor,regressionFrame$age,regressionFrame$IQR,regressionFrame$medianGlu,regressionFrame$admissionDurationDays,regressionFrame$diabetesDuration)
	#colnames(regressionTestSet)<-c("deadAt1yearFactor","age","IQR","medianGlu","admissionDurationDays","diabetesDuration")
	
	#hypoRegressionTestSet<-data.frame(regressionFrame$minGlu4Factor,regressionFrame$age,regressionFrame$diabetesDuration,regressionFrame$IQR,regressionFrame$medianGlu,regressionFrame$admissionDurationDays)
	#colnames(hypoRegressionTestSet)<-c("minGlu4Factor","age","diabetesDuration","IQR","medianGlu","admissionDurationDays")
	
	
	#model <- glm(minGlu4Factor ~.,family=binomial(link='logit'),data=hypoRegressionTestSet)
			
	##### treat IQR as categorical
	divisions<-divisionsOfIQR
	divisionN<-(divisions+1)
	# regressionTestSet$IQRcat1<-cut(regressionTestSet$IQR, seq(0,29,1),labels=c(1:29))
#	regressionTestSet$IQRcat1<-cut(regressionTestSet$IQR, c(quantile(regressionTestSet$IQR)[1],quantile(regressionTestSet$IQR)[2],quantile(regressionTestSet$IQR)[3],quantile(regressionTestSet$IQR)[4],quantile(regressionTestSet$IQR)[5]),labels=c(1:4))
	regressionTestSet$IQRcat1<-cut(regressionTestSet$IQR, quantile(regressionTestSet$IQR,prob = seq(0, 1, length = divisionN), type = 5), labels=c(1:divisions))
	regressionTestSet$IQR<-NULL
	
	regressionTestSetNoDrugData<-data.frame(regressionTestSet$deadAt1yearFactor, regressionTestSet$age, regressionTestSet$diabetesDuration, regressionTestSet$medianGlu, regressionTestSet$admissionDurationDays, regressionTestSet$IQRcat1)
	colnames(regressionTestSetNoDrugData)<-c("deadAt1yearFactor","age","diabetesDuration","medianGlu","admissionDurationDays","IQRcat1")
	
	
	model <- glm(deadAt1yearFactor ~.,family=binomial(link='logit'),data=regressionTestSetNoDrugData)
	exp(cbind(OR = coef(model),confint(model)))
	
	# test the overall effect of IQR by factor
	# library(aod)
	# wald.test(b = coef(model), Sigma = vcov(model), Terms = 6:21)
	
	# calculate the predicted probability of outcome (death at x years) for each mmol/L category of IQR (rankP)
 	newdata1<-with(regressionTestSet, data.frame(age=median(age),medianGlu=median(medianGlu),admissionDurationDays=median(admissionDurationDays),diabetesDuration=median(diabetesDuration),IQRcat1=factor(c(1:divisions))))
 	newdata1$rankP<-predict(model,newdata=newdata1,type="response")
	
	plot(as.numeric(newdata1$IQRcat1),newdata1$rankP,pch=16,cex=1)
	lines(as.numeric(newdata1$IQRcat1),newdata1$rankP,lwd=2)
		# drug prescription
		drugSubPlot<-as.data.frame(matrix(nrow=nlevels(regressionTestSet$IQRcat1),ncol=5))
		colnames(drugSubPlot)<-c("SU","Ins","MF","Other","SUorIns")
		
		for (jj in seq(1,nlevels(regressionTestSet$IQRcat1),1)) {
			drugSub<-subset(regressionTestSet,IQRcat1==jj); 
			nSU<-sum(drugSub$preSU);	propSU<-(nSU/nrow(drugSub)); drugSubPlot$SU[jj]<-propSU
			nIns<-sum(drugSub$preIns);	propIns<-(nIns/nrow(drugSub)); drugSubPlot$Ins[jj]<-propIns
			nMF<-sum(drugSub$preMF);	propMF<-(nMF/nrow(drugSub)); drugSubPlot$MF[jj]<-propMF
			nOther<-sum(drugSub$preOther);	propOther<-(nOther/nrow(drugSub)); drugSubPlot$Other[jj]<-propOther
			
			drugSub$SUorIns<-ifelse((drugSub$preSU | drugSub$preIns),1,0)
			nSUorIns<-sum(drugSub$SUorIns); propSUorIns<-(nSUorIns/nrow(drugSub)); drugSubPlot$SUorIns[jj]<-propSUorIns
			
			}
		
		par(new=TRUE)
		plot(as.numeric(newdata1$IQRcat1),drugSubPlot$SU,col="red",pch=16,cex=0.8,ylim=c(0,1),axes=F);lines(as.numeric(newdata1$IQRcat1),drugSubPlot$SU,col="red")
		axis(side=4)
		points(as.numeric(newdata1$IQRcat1),drugSubPlot$Ins,col="blue",pch=16,cex=0.8);lines(as.numeric(newdata1$IQRcat1),drugSubPlot$Ins,col="blue")
		points(as.numeric(newdata1$IQRcat1),drugSubPlot$MF,col="green",pch=16,cex=0.8);lines(as.numeric(newdata1$IQRcat1),drugSubPlot$MF,col="green")
		points(as.numeric(newdata1$IQRcat1),drugSubPlot$Other,col="gray",pch=16,cex=0.8);lines(as.numeric(newdata1$IQRcat1),drugSubPlot$Other,col="gray")
		points(as.numeric(newdata1$IQRcat1),drugSubPlot$SUorIns,col="orange",pch=16,cex=0.8);lines(as.numeric(newdata1$IQRcat1),drugSubPlot$SUorIns,col="orange",lty=3)
		
		
		
		
	# odds ratio plot with pvals from the logit
	 ODplot<-exp(cbind(OR = coef(model),confint(model)))
	 pvals<-as.data.frame(summary(model)$coefficients[,4])[,1]
	# 
	plot(c(1:length(6:nrow(ODplot))),ODplot[6:nrow(ODplot),1],ylim=c(min(ODplot[6:nrow(ODplot),2]),max(ODplot[6:nrow(ODplot),3])))
	lines(ODplot[6:nrow(ODplot),1])
	#
	 lines(ODplot[6:nrow(ODplot),2],lty=3)
	#
	 lines(ODplot[6:nrow(ODplot),3],lty=3)
	 
	 text(c(1:length(6:nrow(ODplot))),ODplot[6:nrow(ODplot),1],round(pvals[6:nrow(ODplot)],5),pos=3,cex=0.5)
	 
	# for (p in seq(6,6+length(6:nrow(ODplot))-1,1)) {
	#	 print(p)
	#	 # if (pvals[p]<0.05) {text(c(1:length(6:nrow(ODplot))),ODplot[6:nrow(ODplot),1],"*",pos=3,cex=1)}
	# }
	 
	
		plotfilename <- paste("../GlCoSy/GGCdrug/predictedProbability vs medianGlu|age|admissionDuration by IQR quantiles: ",divisions," divisions.pdf",sep="")
		pdf(plotfilename, width=16, height=9)
	
			# plot impact of IQR for range of medianGlu values
			newdata2 <- with(regressionTestSet, data.frame(medianGlu = rep(seq(from = 0, to = 28, length.out = 100),
				divisions), age = median(age), admissionDurationDays=median(admissionDurationDays), diabetesDuration=median(diabetesDuration), IQRcat1 = factor(rep(c(1:divisions), each = 100))))
	
			newdata3 <- cbind(newdata2, predict(model, newdata = newdata2, type = "link",se = TRUE))
	
				newdata3 <- within(newdata3, {
				    PredictedProb <- plogis(fit)
				    LL <- plogis(fit - (1.96 * se.fit))
				    UL <- plogis(fit + (1.96 * se.fit))
				})
		
			ggplot(newdata3, aes(x = medianGlu, y = PredictedProb)) + geom_ribbon(aes(ymin = LL,
			    ymax = UL, fill = IQRcat1), alpha = 0.1) + geom_line(aes(colour = IQRcat1),
			    size = 0.5)
		
		
				# plot impact of IQR for range of age values
				newdata2 <- with(regressionTestSet, data.frame(age = rep(seq(from = 0, to = 100, length.out = 100),
					divisions), medianGlu = median(medianGlu), admissionDurationDays=median(admissionDurationDays), diabetesDuration=median(diabetesDuration), IQRcat1 = factor(rep(c(1:divisions), each = 100))))
	
				newdata3 <- cbind(newdata2, predict(model, newdata = newdata2, type = "link",se = TRUE))
	
					newdata3 <- within(newdata3, {
					    PredictedProb <- plogis(fit)
					    LL <- plogis(fit - (1.96 * se.fit))
					    UL <- plogis(fit + (1.96 * se.fit))
					})
		
				ggplot(newdata3, aes(x = age, y = PredictedProb)) + geom_ribbon(aes(ymin = LL,
				    ymax = UL, fill = IQRcat1), alpha = 0.2) + geom_line(aes(colour = IQRcat1),
				    size = 1)
			
						# plot impact of IQR for range of admissionDuration values
						newdata2 <- with(regressionTestSet, data.frame(admissionDurationDays = rep(seq(from = 0, to = 10, length.out = 100),
							divisions), medianGlu = median(medianGlu), age = median(age), diabetesDuration=median(diabetesDuration), IQRcat1 = factor(rep(c(1:divisions), each = 100))))
	
						newdata3 <- cbind(newdata2, predict(model, newdata = newdata2, type = "link",se = TRUE))
	
							newdata3 <- within(newdata3, {
							    PredictedProb <- plogis(fit)
							    LL <- plogis(fit - (1.96 * se.fit))
							    UL <- plogis(fit + (1.96 * se.fit))
							})
		
						ggplot(newdata3, aes(x = admissionDurationDays, y = PredictedProb)) + geom_ribbon(aes(ymin = LL,
						    ymax = UL, fill = IQRcat1), alpha = 0.2) + geom_line(aes(colour = IQRcat1),
						    size = 1)
							
		dev.off()
	
		
		
		library(ROCR)
		p <- predict(model, newdata=subset(test,select=c(2,3,4,5)), type="response")
		pr <- prediction(p, test$deadAt1yearFactor)
		prf <- performance(pr, measure = "tpr", x.measure = "fpr")
		plot(prf)
		
		auc <- performance(pr, measure = "auc")
		auc <- auc@y.values[[1]]
		auc
	
}


regressionFromAdmission<-function(inputFrame,deadAtNyears) {
  #regressionFrame<-admissionsWithDrugData_testSet
  #regressionFrame<-subset(admissionsWithDrugData_testSet)
  regressionFrame<-inputFrame
  regressionFrame$diabetesDuration<-regressionFrame$dateplustime1-regressionFrame$diagnosisDateUnix
  regressionFrame<-subset(regressionFrame,diagnosisDateUnix>(-2208988800))
  # deadAtNyears<-1
  # divisionsOfIQR<-5
  
  regressionFrame$preSUFactor<-as.factor(regressionFrame$preSU)
  regressionFrame$preInsFactor<-as.factor(regressionFrame$preIns)
  regressionFrame$preMFFactor<-as.factor(regressionFrame$preMF)
  
  
  
  lastDOD<-max(regressionFrame$deathDateUnix)
  
  regressionFrame$timeToDeath<-ifelse(regressionFrame$deathEvent==1,(regressionFrame$deathDateUnix-regressionFrame$dateplustime1),0)
  regressionFrame$deadAtNyears<-ifelse(regressionFrame$timeToDeath>0 & regressionFrame$timeToDeath<(deadAtNyears*365.25*60*60*24),1,0)
  
  regressionFrame$useInRegressionAnalysis<-ifelse(lastDOD-regressionFrame$dateplustime1>=(deadAtNyears*365.25*60*60*24),1,0)
  
  
  logitRegression<-glm(deadAtNyears ~preSUFactor + preInsFactor + preMFFactor + age + diabetesDuration,data=subset(regressionFrame,useInRegressionAnalysis==1))
#  print(sum(regressionFrame$useInRegressionAnalysis))
#  print(summary(logitRegression))
  
  # return estimate and p val - and then plot for differing follow up times.
  outputList<-list(logitRegression$coefficients[["preSUFactor1"]],logitRegression$coefficients[["preInsFactor1"]],logitRegression$coefficients[["preMFFactor1"]],sum(regressionFrame$useInRegressionAnalysis))
  return(outputList)
  
}


regressionFromAdmissionSingleParameter<-function(inputFrame,deadAtNyears,testParameter) {
  #regressionFrame<-admissionsWithDrugData_testSet
  #regressionFrame<-subset(admissionsWithDrugData_testSet)
  regressionFrame<-inputFrame
  regressionFrame$diabetesDuration<-regressionFrame$dateplustime1-regressionFrame$diagnosisDateUnix
  regressionFrame<-subset(regressionFrame,diagnosisDateUnix>(-2208988800))
  # deadAtNyears<-0.2
  # divisionsOfIQR<-5
  
  testFactor<-as.factor(testParameter)
  
#  regressionFrame$preSUFactor<-as.factor(regressionFrame$preSU)
#  regressionFrame$preInsFactor<-as.factor(regressionFrame$preIns)
 # regressionFrame$preMFFactor<-as.factor(regressionFrame$preMF)
  
  
  
  lastDOD<-max(regressionFrame$deathDateUnix)
  
  regressionFrame$timeToDeath<-ifelse(regressionFrame$deathEvent==1,(regressionFrame$deathDateUnix-regressionFrame$dateplustime1),0)
  regressionFrame$deadAtNyears<-ifelse(regressionFrame$timeToDeath>0 & regressionFrame$timeToDeath<(deadAtNyears*365.25*60*60*24),1,0)
  
  regressionFrame$useInRegressionAnalysis<-ifelse(lastDOD-regressionFrame$dateplustime1>=(deadAtNyears*365.25*60*60*24),1,0)
  
  
  logitRegression<-glm(deadAtNyears ~ testFactor + age + diabetesDuration,data=subset(regressionFrame,useInRegressionAnalysis==1))
  #  print(sum(regressionFrame$useInRegressionAnalysis))
  #  print(summary(logitRegression))
  
  # return estimate and p val - and then plot for differing follow up times.
  outputList<-list(logitRegression$coefficients[["testFactor"]],sum(regressionFrame$useInRegressionAnalysis))
  return(outputList)
  
}



process_eGFR<-function(eGFRdata) {
	d<-strptime(eGFRdata$DataItemDate,"%d/%m/%Y")
	dd<-as.numeric(d)
	eGFRdata$dateplustime1<-dd
		
	eGFRdata$DataItemValue[eGFRdata$DataItemValue == '>60'] <- '60'
	eGFRdata$DataItemValue[eGFRdata$DataItemValue == '>60.0'] <- '60'
	eGFRdata$DataItemValue[eGFRdata$DataItemValue == '> 60'] <- '60'
	eGFRdata$DataItemValue[eGFRdata$DataItemValue == '> 60.0'] <- '60'
	eGFRdata$DataItemValue[eGFRdata$DataItemValue == '> 60\n'] <- '60'
	eGFRdata$DataItemValue[eGFRdata$DataItemValue == '>60\n'] <- '60'
	eGFRdata$DataItemValue[eGFRdata$DataItemValue == '>59'] <- '60'
	eGFRdata$DataItemValue[eGFRdata$DataItemValue == '> 59'] <- '60'
	eGFRdata$DataItemValue[eGFRdata$DataItemValue == '>59\n'] <- '60'
	
	eGFRdata$DataItemValue[eGFRdata$DataItemValue == '<60'] <- '59'
	
	
	eGFRdata$eGFR<-as.numeric(levels(eGFRdata$DataItemValue))[eGFRdata$DataItemValue]
	eGFRdata$eGFR<-ifelse(eGFRdata$eGFR>60,60,eGFRdata$eGFR)
	
	eGFRdata$eGFR[is.na(eGFRdata$eGFR)] <- 0
	eGFRdata<-subset(eGFRdata,eGFR>0)

	return(eGFRdata)
}

findRecenteGFR<-function(admissionsWithDrugData_testSet,eGFRdata,preWindowDays,postWindowDays) {
	eGFRdata<-data.table(eGFRdata)
	preWindowSeconds<-preWindowDays*(60*60*24)
	postWindowSeconds<-postWindowDays*(60*60*24)
	output_eGFR<-as.data.frame(matrix(0,ncol=1,nrow=nrow(admissionsWithDrugData_testSet)));colnames(output_eGFR)<-c("eGFR")
	
	for (k in seq(1,nrow(admissionsWithDrugData_testSet),1)) {
		
		if (k%%100==0) {print(k)}
		
		testID<-admissionsWithDrugData_testSet$ID[k]
		testDateOfAdmission<-admissionsWithDrugData_testSet$dateplustime1[k]
		
		eGFRset<-eGFRdata[PatId==testID & (dateplustime1>(testDateOfAdmission-preWindowSeconds)) & (dateplustime1<(testDateOfAdmission+postWindowSeconds))]
		eGFRset<-eGFRset[order(eGFRset$dateplustime1),]   # order by datetime
		
		if (nrow(eGFRset)>0) {output_eGFR$eGFR[k]=tail(eGFRset,1)$eGFR}
	}
	
	return(output_eGFR)
}

# work out male/female	
maleFemaleFunction<-function(IDs) {
	sexDigit<-ifelse(nchar(IDs)==10,substr(IDs,9,9),0)
	sexDigit<-ifelse(nchar(IDs)==9,substr(IDs,8,8),sexDigit)
	sexDigit<-as.numeric(sexDigit)
	sexIdent<-ifelse(sexDigit%%2==0,0,1)
	totalFemale<-length(sexIdent) - sum(sexIdent)
	totalMale<-sum(sexIdent)
	print("total Female:"); print(totalFemale)
	print("total Male"); print(totalMale)
	print("proportion Female"); print(totalFemale / length(sexIdent))
	print("proportion Male"); print(totalMale / length(sexIdent))
}


########

# read data in here
OPfilename <- paste("~/R/GlCoSy/source/admissionsWithDrugData.csv",sep="")
admissionsWithDrugData<-read.csv(OPfilename, header=TRUE , sep="," , row.names=NULL)

# read in eGFR data here
eGFRdataName<-paste("../GlCoSy/source/GGC_eGFR.txt",sep="")
eGFRdata<-read.csv(eGFRdataName, header=TRUE , sep="," , row.names=NULL)
eGFRdata<-process_eGFR(eGFRdata)

# HBdataName<-paste("../GlCoSy/source/GGC_HBA1c.txt",sep="")
# HBdata<-read.csv(HBdataName, header=TRUE , sep="," , row.names=NULL)

output_eGFR<-findRecenteGFR(admissionsWithDrugData,eGFRdata,91,1)
colnames(output_eGFR)<-c("eGFR")
admissionsWithDrugData$eGFR<-output_eGFR$eGFR

## if need to load pre calculated eGFR -90 to +1 day
if (loadeGFR_minus90_to_plus1==1) {
	outputName<-paste("../GlCoSy/source/admissionsWithDrugData_eGFR_91dayRunIn.csv",sep="")
	outputName<-paste("../GlCoSy/source/admissionsWithDrugData_eGFR_367dayRunIn_redo.csv",sep="")
#	write.table(admissionsWithDrugData,file=outputName,sep=",",append=F,col.names=T)
	admissionsWithDrugData<-read.csv(outputName, header=TRUE , sep="," , row.names=NULL)
}

admissionsWithDrugData_testSet<-subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus")
# admissionsWithDrugData_testSet<-subset(admissionsWithDrugData,DiabetesMellitusType_Mapped=="Type 1 Diabetes Mellitus")

admissionsWithDrugData_testSet.full<-admissionsWithDrugData_testSet
admissionsWithDrugData_testSet<-admissionsWithDrugData_testSet[1:1000,]

################################################################################################################################
################################################################################################################################
######## plots for duk 2017 drug abstract

set.seed(42)

## drug test 1
cases<-subset(admissionsWithDrugData_testSet,preSU==1 & preIns==0 & preMF==1)
controlPool<-subset(admissionsWithDrugData_testSet,preSU==0 & preIns==0 & preMF==1)
controlSelection(cases,controlPool,paste("group1",sep=""),1,80)


## drug test 2
cases<-subset(admissionsWithDrugData_testSet,preSU==0 & preIns==1 & preMF==1)
controlPool<-subset(admissionsWithDrugData_testSet,preSU==0 & preIns==0 & preMF==1)
controlSelection(cases,controlPool,paste("group2",sep=""),1,90)

## drug test 3
cases<-subset(admissionsWithDrugData_testSet,preSU==0 & preIns==0 & preMF==1)
controlPool<-subset(admissionsWithDrugData_testSet,preSU==0 & preIns==0 & preMF==0)
controlSelection(cases,controlPool,paste("group3",sep=""),1,90)

## drug test 4
cases<-subset(admissionsWithDrugData_testSet,preSU==0 & preIns==1 & preMF==1)
controlPool<-subset(admissionsWithDrugData_testSet,preSU==1 & preIns==0 & preMF==1)
controlSelection(cases,controlPool,paste("group4",sep=""),1,90)

  # test all cases (su1, ins0, mf1) against controls (su0, ins0, mf1) for survival
  admissionsWithDrugData$su_mf_case <- ifelse(admissionsWithDrugData$preSU == 1 & admissionsWithDrugData$preIns == 0 & admissionsWithDrugData$preMF == 1, 1, 0)
  admissionsWithDrugData$mf_case <- ifelse(admissionsWithDrugData$preSU == 0 & admissionsWithDrugData$preIns == 0 & admissionsWithDrugData$preMF == 1, 1, 0)
  
  survivalTest <- subset(admissionsWithDrugData, su_mf_case == 1 | mf_case == 1)
  
  simpleSurvivalPlot_SUcaseVScontrol(survivalTest,0.6,90)


################################################################################################################################
################################################################################################################################
######## regression analysis

interestFrame<-admissionsWithDrugData_testSet
interestFrame<-subset(admissionsWithDrugData_testSet,age>60)
interestFrame<-subset(admissionsWithDrugData,eGFR>=60 & DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus")



min<-0.002737851 # single day as fraction of a year
max<-3
int<-0.002737851

reportMatrix<-as.data.frame(matrix(0,nrow=length(seq(min,max,int)),ncol=6))
colnames(reportMatrix)<-c("index","interval","SU","Ins","MF","n")

reportMatrix$index<-c(1:length(seq(min,max,int)))
reportMatrix$interval<-seq(min,max,int)

    for (j in seq(1,length(seq(1,length(seq(min,max,int)),1)),1)) {
      testInterval<-subset(reportMatrix,index==j)$interval
      logitEstimates<-regressionFromAdmission(interestFrame,testInterval)
      
      reportMatrix$SU[j]<-logitEstimates[[1]]
      reportMatrix$Ins[j]<-logitEstimates[[2]]
      reportMatrix$MF[j]<-logitEstimates[[3]]
      reportMatrix$n[j]<-logitEstimates[[4]]
    
      
    }

plot(reportMatrix$SU,col="black",pch=16,ylim=c(-0.05,0.04))
points(reportMatrix$Ins,col="red",pch=16)
points(reportMatrix$MF,col="blue",pch=16)
abline(0,0,lwd=2,col="red")

plot(reportMatrix$SU,col="black",pch=16,ylim=c(-0.06,0.06),xlim=c(0,365))


# IQR 50% unmatched (for paper)
ageWindow<- 200 ; medianWindow<-100; admissionDurationWindow<-1000; diabetesDurationWindow<-1000
cases<-subset(admissionsWithDrugData_testSet, IQR>=quantile(admissionsWithDrugData_testSet$IQR)[3])
controlPool<-subset(admissionsWithDrugData_testSet, IQR<quantile(admissionsWithDrugData_testSet$IQR)[3])
controlSelection(cases,controlPool,paste("MortPaper.IQRtop50 vs bottom50. day90postD. unmatched0. cutoff:",quantile(admissionsWithDrugData_testSet$IQR)[3],sep=""),1,90)

# IQR 50% cut
cases<-subset(admissionsWithDrugData_testSet, IQR>=quantile(admissionsWithDrugData_testSet$IQR)[3])
controlPool<-subset(admissionsWithDrugData_testSet, IQR<quantile(admissionsWithDrugData_testSet$IQR)[3])
controlSelection(cases,controlPool,paste("insert",sep=""),1,90)

# IQR 50% cut with no hypos during admission
cases<-subset(admissionsWithDrugData_testSet, IQR>=quantile(admissionsWithDrugData_testSet$IQR)[3] & minGlu>=4)
controlPool<-subset(admissionsWithDrugData_testSet, IQR<quantile(admissionsWithDrugData_testSet$IQR)[3] & minGlu>=4)
controlSelection(cases,controlPool,paste("Mortality Paper Plots: top half IQR vs bottom half IQR. No hypos. cutoff:",quantile(admissionsWithDrugData_testSet$IQR)[3],sep=""),1)

# IQR 50% matched for hypo episodes per admission
cases<-subset(admissionsWithDrugData_testSet, IQR>=quantile(admissionsWithDrugData_testSet$IQR)[3])
controlPool<-subset(admissionsWithDrugData_testSet, IQR<quantile(admissionsWithDrugData_testSet$IQR)[3])
controlSelection(cases,controlPool,paste("Mortality Paper Plots: top IQR vs bottom IQR. matched hypEps. cutoff:",quantile(admissionsWithDrugData_testSet$IQR)[3],sep=""),4,90)

# IQR 50% matched for hypo episodes per admission - presence or absense of hypoglycaemia based on categorical presence/absence of hypoglycaemia
cases<-subset(admissionsWithDrugData_testSet, IQR>=quantile(admissionsWithDrugData_testSet$IQR)[3])
controlPool<-subset(admissionsWithDrugData_testSet, IQR<quantile(admissionsWithDrugData_testSet$IQR)[3])
controlSelection(cases,controlPool,paste("Mortality Paper Plots: top IQR vs bottom IQR. matched hypEps - categoricalHypos. cutoff:",quantile(admissionsWithDrugData_testSet$IQR)[3],sep=""),6,90)


# IQR 75% cut
cases<-subset(admissionsWithDrugData_testSet, IQR>=quantile(admissionsWithDrugData_testSet$IQR)[4])
controlPool<-subset(admissionsWithDrugData_testSet, IQR<quantile(admissionsWithDrugData_testSet$IQR)[4])
controlSelection(cases,controlPool,paste("Mortality Paper Plots: top quartile IQR vs rest IQR. cutoff:",quantile(admissionsWithDrugData_testSet$IQR)[4],sep=""),1,90)

# IQR 75% cut - eGFR>60
cases<-subset(admissionsWithDrugData_testSet, IQR>=quantile(admissionsWithDrugData_testSet$IQR)[4] & (eGFR>=60))
controlPool<-subset(admissionsWithDrugData_testSet, IQR<quantile(admissionsWithDrugData_testSet$IQR)[4] & (eGFR>=60))
controlSelection(cases,controlPool,paste("Mortality Paper Plots: top quartile IQR vs rest IQR. eGFR>60. cutoff:",quantile(admissionsWithDrugData_testSet$IQR)[4],sep=""),1,90)

# IQR 75% cut - eGFR<60
cases<-subset(admissionsWithDrugData_testSet, IQR>=quantile(admissionsWithDrugData_testSet$IQR)[4] & (eGFR<60) & (eGFR>0))
controlPool<-subset(admissionsWithDrugData_testSet, IQR<quantile(admissionsWithDrugData_testSet$IQR)[4] & (eGFR<60)& (eGFR>0))
controlSelection(cases,controlPool,paste("Mortality Paper Plots: top quartile IQR vs rest IQR. eGFR<60 >0. cutoff:",quantile(admissionsWithDrugData_testSet$IQR)[4],sep=""),1,90)

# IQR 75% cut - matching for eGFR using function
cases<-subset(admissionsWithDrugData_testSet, IQR>=quantile(admissionsWithDrugData_testSet$IQR)[4] & eGFR>0)
controlPool<-subset(admissionsWithDrugData_testSet, IQR<quantile(admissionsWithDrugData_testSet$IQR)[4] & eGFR>0)
controlSelection(cases,controlPool,paste("Mortality Paper Plots: top quartile IQR vs rest IQR. eGFR match window20. cutoff:",quantile(admissionsWithDrugData_testSet$IQR)[4],sep=""),8,90)


# IQR top quartile with no hypos during admission
cases<-subset(admissionsWithDrugData_testSet, IQR>=quantile(admissionsWithDrugData_testSet$IQR)[4] & minGlu>=4)
controlPool<-subset(admissionsWithDrugData_testSet, IQR<quantile(admissionsWithDrugData_testSet$IQR)[4])
controlSelection(cases,controlPool,paste("Mortality Paper Plots: top quartile IQR vs rest IQR. No hypos. cutoff:",quantile(admissionsWithDrugData_testSet$IQR)[4],sep=""),1)

# IQR top quartile with no hypos during admission
cases<-subset(admissionsWithDrugData_testSet, IQR>=quantile(admissionsWithDrugData_testSet$IQR)[4] & minGlu>=4)
controlPool<-subset(admissionsWithDrugData_testSet, IQR<quantile(admissionsWithDrugData_testSet$IQR)[4])
controlSelection(cases,controlPool,paste("Mortality Paper Plots: top quartile IQR vs rest IQR. No hypos. cutoff:",quantile(admissionsWithDrugData_testSet$IQR)[4],sep=""),6,90)

# IQR top thirds vs bottom 2 thirds (likely split to demonstrate largest difference)
IQRthirds<-quantile(admissionsWithDrugData_testSet$IQR,prob = seq(0, 1, length = 4), type = 7)
cases<-subset(admissionsWithDrugData_testSet, IQR>=IQRthirds[3])
controlPool<-subset(admissionsWithDrugData_testSet, IQR<IQRthirds[3])
controlSelection(cases,controlPool,paste("Mortality Paper Plots: top thirds IQR vs bottom thirds IQR. cutoff:",IQRthirds[3],sep=""),1,90)

# IQR top thirds vs bottom 2 thirds (likely split to demonstrate largest difference) SU or Ins
IQRthirds<-quantile(admissionsWithDrugData_testSet$IQR,prob = seq(0, 1, length = 4), type = 7)
cases<-subset(admissionsWithDrugData_testSet, (IQR>=IQRthirds[3] & (preSU==1 | preIns==1)))
controlPool<-subset(admissionsWithDrugData_testSet, (IQR<IQRthirds[3] & (preSU==1 | preIns==1)))
controlSelection(cases,controlPool,paste("Ins or SU. Mortality Paper Plots: top thirds IQR vs bottom thirds IQR. cutoff:",IQRthirds[3],sep=""),6,90)

# IQR top thirds vs bottom 2 thirds (likely split to demonstrate largest difference) - no SU or Ins
IQRthirds<-quantile(admissionsWithDrugData_testSet$IQR,prob = seq(0, 1, length = 4), type = 7)
cases<-subset(admissionsWithDrugData_testSet, (IQR>=IQRthirds[3] & (preSU!=1 & preIns!=1)))
controlPool<-subset(admissionsWithDrugData_testSet, (IQR<IQRthirds[3] & (preSU!=1 & preIns!=1)))
controlSelection(cases,controlPool,paste("No Ins or SU. Mortality Paper Plots: top thirds IQR vs bottom thirds IQR. cutoff:",IQRthirds[3],sep=""),6,90)

# IQR top quintile vs bottom quintiles - on SU and MF
IQRquint<-quantile(admissionsWithDrugData_testSet$IQR,prob = seq(0, 1, length = 6), type = 7)
cases<-subset(admissionsWithDrugData_testSet, (IQR>=IQRquint[4] & (preSU==1 & preMF==1)))
controlPool<-subset(admissionsWithDrugData_testSet, (IQR<IQRquint[4] & (preSU==1 & preMF==1)))
controlSelection(cases,controlPool,paste("SU+MF. Mortality Paper Plots: top quintile IQR vs bottom quintiles IQR. cutoff:",IQRquint[4],sep=""),6,90)


# IQR top vs bottom quartile
cases<-subset(admissionsWithDrugData_testSet, IQR>=quantile(admissionsWithDrugData_testSet$IQR)[4])
controlPool<-subset(admissionsWithDrugData_testSet, IQR<quantile(admissionsWithDrugData_testSet$IQR)[2])
controlSelection(cases,controlPool,paste("Mortality Paper Plots: top quart IQR vs bottom quart IQR. cutoff:",quantile(admissionsWithDrugData_testSet$IQR)[3],sep=""),1)

# IQR 2nd vs 3rd quartile cut
cases<-subset(admissionsWithDrugData_testSet, IQR>=quantile(admissionsWithDrugData_testSet$IQR)[3] & IQR<quantile(admissionsWithDrugData_testSet$IQR)[4])
controlPool<-subset(admissionsWithDrugData_testSet, IQR>=quantile(admissionsWithDrugData_testSet$IQR)[2] & IQR<quantile(admissionsWithDrugData_testSet$IQR)[3])
controlSelection(cases,controlPool,paste("Mortality Paper Plots: 3rd vs 2nd IQR quartile. cutoffs:",quantile(admissionsWithDrugData_testSet$IQR)[2]," : ",quantile(admissionsWithDrugData_testSet$IQR)[3]," : ",quantile(admissionsWithDrugData_testSet$IQR)[4],sep=""),1)

# minGlu < 4
cases<-subset(admissionsWithDrugData_testSet, minGlu<4)
controlPool<-subset(admissionsWithDrugData_testSet, minGlu>=4)
controlSelection(cases,controlPool,paste("Mortality Paper Plots: cases: minGlu<4",sep=""),1,90)

# minGlu < 4 matched for median and IQR
cases<-subset(admissionsWithDrugData_testSet, minGlu<4)
controlPool<-subset(admissionsWithDrugData_testSet, minGlu>=4)
controlSelection(cases,controlPool,paste("Mortality Paper Plots: cases: minGlu<4. matched for wIQR-",sep=""),2,90)


# minGlu < 4 matched for median and IQR all eGFR>=60
cases<-subset(admissionsWithDrugData_testSet, minGlu<4 & eGFR>=60)
controlPool<-subset(admissionsWithDrugData_testSet, minGlu>=4 & eGFR>=60)
controlSelection(cases,controlPool,paste("Mortality Paper Plots: cases: minGlu<4. matched for IQR. eGFR>60-",sep=""),2,90)


# minGlu < 3
cases<-subset(admissionsWithDrugData_testSet, minGlu<3)
controlPool<-subset(admissionsWithDrugData_testSet, minGlu>=3)
controlSelection(cases,controlPool,paste("Mortality Paper Plots: cases: minGlu<3",sep=""),1)

# maxGlu > 15
cases<-subset(admissionsWithDrugData_testSet, maxGlu>=15)
controlPool<-subset(admissionsWithDrugData_testSet, maxGlu<15)
controlSelection(cases,controlPool,paste("Mortality Paper Plots: cases: maxGlu>=15",sep=""),1)

# maxGlu > 15 matched for median and IQR
cases<-subset(admissionsWithDrugData_testSet, maxGlu>=15)
controlPool<-subset(admissionsWithDrugData_testSet, maxGlu<15)
controlSelection(cases,controlPool,paste("Mortality Paper Plots: cases: maxGlu>=15. matched for IQR",sep=""),2)

# maxGlu > 20
cases<-subset(admissionsWithDrugData_testSet, maxGlu>=20)
controlPool<-subset(admissionsWithDrugData_testSet, maxGlu<20)
controlSelection(cases,controlPool,paste("Mortality Paper Plots: cases: maxGlu>=20",sep=""),1)

# maxGlu > 25
cases<-subset(admissionsWithDrugData_testSet, maxGlu>=26)
controlPool<-subset(admissionsWithDrugData_testSet, maxGlu<26)
controlSelection(cases,controlPool,paste("Mortality Paper Plots: cases: maxGlu>=26",sep=""),1)

# median 50% cut - matched for IQR
cases<-subset(admissionsWithDrugData_testSet, medianGlu>=quantile(admissionsWithDrugData_testSet$medianGlu)[3])
controlPool<-subset(admissionsWithDrugData_testSet, medianGlu<quantile(admissionsWithDrugData_testSet$medianGlu)[3])
controlSelection(cases,controlPool,paste("Mortality Paper Plots: top half median vs bottom half median (IQR matched).cutoff:",quantile(admissionsWithDrugData_testSet$medianGlu)[3],sep=""),3,90)

# median 50% cut - matched for IQR - no hypos in either group
cases<-subset(admissionsWithDrugData_testSet, medianGlu>=quantile(admissionsWithDrugData_testSet$medianGlu)[3] & minGlu>=4)
controlPool<-subset(admissionsWithDrugData_testSet, medianGlu<quantile(admissionsWithDrugData_testSet$medianGlu)[3] & minGlu>=4)
controlSelection(cases,controlPool,paste("Mortality Paper Plots: top median vs bottom median (IQR matched.no hypos).cutoff:",quantile(admissionsWithDrugData_testSet$medianGlu)[3],sep=""),3)

# median 50% cut - matched for IQR and hypo episodes
cases<-subset(admissionsWithDrugData_testSet, medianGlu>=quantile(admissionsWithDrugData_testSet$medianGlu)[3])
controlPool<-subset(admissionsWithDrugData_testSet, medianGlu<quantile(admissionsWithDrugData_testSet$medianGlu)[3])
controlSelection(cases,controlPool,paste("run2 Mortality Paper Plots: medianHalf IQR hypo matched",quantile(admissionsWithDrugData_testSet$medianGlu)[3],sep=""),5,90)

# median top quartile - matched for IQR
cases<-subset(admissionsWithDrugData_testSet, medianGlu>=quantile(admissionsWithDrugData_testSet$medianGlu)[4])
controlPool<-subset(admissionsWithDrugData_testSet, medianGlu<quantile(admissionsWithDrugData_testSet$medianGlu)[4])
controlSelection(cases,controlPool,paste("Mortality Paper Plots: top quartile medianGlu vs rest (IQR matched). cutoff:",quantile(admissionsWithDrugData_testSet$medianGlu)[4],sep=""),3)

# median bottom quartile - matched for IQR
cases<-subset(admissionsWithDrugData_testSet, medianGlu<quantile(admissionsWithDrugData_testSet$medianGlu)[2])
controlPool<-subset(admissionsWithDrugData_testSet, medianGlu>=quantile(admissionsWithDrugData_testSet$medianGlu)[2])
controlSelection(cases,controlPool,paste("Mortality Paper Plots: bottom quartile medianGlu vs rest (IQR matched). cutoff:",quantile(admissionsWithDrugData_testSet$medianGlu)[2],sep=""),3,90)

# median bottom vs top quartile - matched for IQR
cases<-subset(admissionsWithDrugData_testSet, medianGlu<quantile(admissionsWithDrugData_testSet$medianGlu)[2])
controlPool<-subset(admissionsWithDrugData_testSet, medianGlu>=quantile(admissionsWithDrugData_testSet$medianGlu)[4])
controlSelection(cases,controlPool,paste("Mortality Paper Plots: bottom quart medianGlu vs top (IQR matched). cutoffs:",quantile(admissionsWithDrugData_testSet$medianGlu)[2]," & ",quantile(admissionsWithDrugData_testSet$medianGlu)[4],sep=""),3)

# categorical hypo vs no categorical hypo - matched for all including IQR
cases<-subset(admissionsWithDrugData_testSet, minGlu<4)
controlPool<-subset(admissionsWithDrugData_testSet, minGlu>=4)
controlSelection(cases,controlPool,paste("Mortality Paper Plots: hypo vs no hypo (IQR matched):",sep=""),3,90)


# eGFR >60 vs <60
cases<-subset(admissionsWithDrugData_testSet, eGFR<60)
controlPool<-subset(admissionsWithDrugData_testSet, eGFR>=60)
controlSelection(cases,controlPool,paste("Mortality Paper Plots: eGFF>60 vs eGFR<60. cutoff:",quantile(admissionsWithDrugData_testSet$IQR)[4],sep=""),2,90)



# minGlu < 4 matched for median and IQR - not for diabetes duration. gives a 1.1 day difference for hypo vs no hypo
# easier match for iqr and median with 30 day onward vs 90 day onward 
minCBGperAdmission<-4
ageWindow<-5 # (years)
medianWindow<-0.3 # (mmol/L)
admissionDurationWindow<-1 # (days)
diabetesDurationWindow<-4 # (years)
IQRwindow<-0.1 # (mmol/L)
cases<-subset(admissionsWithDrugData_testSet, minGlu<4 & ID_ADMISSIONhypoEpisodes4.60==1)
controlPool<-subset(admissionsWithDrugData_testSet, minGlu>=4)
controlSelection(cases,controlPool,paste("Mortality Paper Plots: cases: minGlu<4. Investigating LOS",sep=""),7,30)

# testing diabetes care numbers. no drug data
minCBGperAdmission<-4
ageWindow<-5 # (years)
medianWindow<-0.4 # (mmol/L)
admissionDurationWindow<-1 # (days)
diabetesDurationWindow<-4 # (years)
IQRwindow<-0.2 # (mmol/L)
cases<-subset(admissionsWithDrugData_testSet, maxGlu>=10)
	print(nrow(cases))
controlPool<-subset(admissionsWithDrugData_testSet, maxGlu<10)
controlSelection(cases,controlPool,paste("maxglu thresh 10. no drug data. ",sep=""),7,30)


# testing diabetes care numbers. insulin or SU. maxGlu 10-13.3 (180-240) vs maxGlu<10
minCBGperAdmission<-4
ageWindow<-5 # (years)
medianWindow<-1 # (mmol/L)
admissionDurationWindow<-1 # (days)
diabetesDurationWindow<-4 # (years)
IQRwindow<-1 # (mmol/L)
cases<-subset(admissionsWithDrugData_testSet, (maxGlu>=10 & maxGlu<13.3) & (preIns==1 | preSU==1))
	print(nrow(cases))
controlPool<-subset(admissionsWithDrugData_testSet, maxGlu<10 & (preIns==1 & preSU==1))
controlSelection(cases,controlPool,paste("maxglu (>10 & <13.3) vs (<10). preIns==1 OR preSU==1. ",sep=""),7,30)

# diabetes care paper
# 1 regression
testSet<-subset(admissionsWithDrugData_testSet,IQR>0)
plot(log(testSet$maxGlu),log(testSet$admissionDurationDays),pch=16,cex=0.3)
fit<-lm(testSet$admissionDurationDays ~ testSet$maxGlu + testSet$age + testSet$IQR + testSet$medianGlu)
fit<-lm(log(testSet$admissionDurationDays) ~ log(testSet$maxGlu) + testSet$age + log(testSet$IQR) + log(testSet$medianGlu))
abline(fit,col="red")

newdata1<-with(testSet, data.frame(age=median(age),medianGlu=median(medianGlu),admissionDurationDays=median(admissionDurationDays),diabetesDuration=median(diabetesDuration),IQR=median(IQR),maxGlu_cat1=factor(c(1,2,3))))
newdata1$rankP<-predict(model,newdata=newdata1,type="response")


# multiFactor combination plots

	SurvivalData<-admissionsWithDrugData_testSet
	SurvivalData<-plotFrame
	
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
	
		plotfilename <- paste("../GlCoSy/GGCdrug/Mortality Paper Plots: multiFactor Survival Plots.pdf",sep="")
		pdf(plotfilename, width=16, height=9)

			# median/IQR - upper/lower 50%
			sfit4<-survfit(Surv(timeToDeathInterval,  shortDeathEvent) ~ (medianGlu>=quantile(SurvivalData$medianGlu)[3])+(IQR>=quantile(SurvivalData$IQR)[3]), data=SurvivalData)
			plotType<-paste("IQR + medianGlu. upper/lower 50% ",sep="")
			shortPlotTitle <- paste(plotType,"\n, time ",(round(shortCensorPeriodStartDay)/DaySeconds)," to ",(round(max(SurvivalData$timeToDeathInterval))/DaySeconds)," days",sep="")
			plot(sfit4,mark.time=F,lty=1:6,conf.int=F,col=c("black","red","blue","green","orange","purple"),main=shortPlotTitle,xlim=c(shortCensorPeriodStartDay,round(max(SurvivalData$timeToDeathInterval))),lwd=3)
			sfit4.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ age+admissionDuration+(medianGlu>=quantile(admissionsWithDrugData_testSet$medianGlu)[3])+(IQR>=quantile(admissionsWithDrugData_testSet$IQR)[3]), data=SurvivalData)
				pVal <- summary(sfit4.coxph)$coef[,5]; HR <- round(exp(coef(sfit4.coxph)),2)
				legendText <- paste("p = ",pVal," | HR = ",HR,sep="")
			summarySurvfit <- summary(sfit4); legendNames <- row.names(summarySurvfit$table)
			legend("bottomleft",c(legendNames),lty=1:6,col=c("black","red","blue","green","orange","purple"),cex=0.8); legend("topright",legendText,cex=0.6)

			# IQR (upper/lower 50%) with/without hypoglycaemia
			sfit4<-survfit(Surv(timeToDeathInterval,  shortDeathEvent) ~ (minGlu<4)+(IQR>=quantile(admissionsWithDrugData_testSet$IQR)[3]), data=SurvivalData)
			plotType<-paste("IQR (upper/lower 50%) with/without hypoglycaemia ",sep="")
			shortPlotTitle <- paste(plotType,"\n, time ",(round(shortCensorPeriodStartDay)/DaySeconds)," to ",(round(max(SurvivalData$timeToDeathInterval))/DaySeconds)," days",sep="")
			plot(sfit4,mark.time=F,lty=1:6,conf.int=F,col=c("black","red","blue","green","orange","purple"),main=shortPlotTitle,xlim=c(shortCensorPeriodStartDay,round(max(SurvivalData$timeToDeathInterval))),lwd=3)
			sfit4.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ age+admissionDuration+(minGlu<4)+(IQR>=quantile(admissionsWithDrugData_testSet$IQR)[3]), data=SurvivalData)
				pVal <- summary(sfit4.coxph)$coef[,5]; HR <- round(exp(coef(sfit4.coxph)),2)
				legendText <- paste("p = ",pVal," | HR = ",HR,sep="")
			summarySurvfit <- summary(sfit4); legendNames <- row.names(summarySurvfit$table)
			legend("bottomleft",c(legendNames),lty=1:6,col=c("black","red","blue","green","orange","purple"),cex=0.8); legend("topright",legendText,cex=0.6)

			# IQR (upper 25% vs rest) with/without hypoglycaemia
			sfit4<-survfit(Surv(timeToDeathInterval,  shortDeathEvent) ~ (minGlu<4)+(IQR>=quantile(admissionsWithDrugData_testSet$IQR)[4]), data=SurvivalData)
			plotType<-paste("IQR (upper 25% vs rest) with/without hypoglycaemia ",sep="")
			shortPlotTitle <- paste(plotType,"\n, time ",(round(shortCensorPeriodStartDay)/DaySeconds)," to ",(round(max(SurvivalData$timeToDeathInterval))/DaySeconds)," days",sep="")
			plot(sfit4,mark.time=F,lty=1:6,conf.int=F,col=c("black","red","blue","green","orange","purple"),main=shortPlotTitle,xlim=c(shortCensorPeriodStartDay,round(max(SurvivalData$timeToDeathInterval))),lwd=3)
			sfit4.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ age+admissionDuration+(minGlu<4)+(IQR>=quantile(admissionsWithDrugData_testSet$IQR)[4]), data=SurvivalData)
				pVal <- summary(sfit4.coxph)$coef[,5]; HR <- round(exp(coef(sfit4.coxph)),2)
				legendText <- paste("p = ",pVal," | HR = ",HR,sep="")
			summarySurvfit <- summary(sfit4); legendNames <- row.names(summarySurvfit$table)
			legend("bottomleft",c(legendNames),lty=1:6,col=c("black","red","blue","green","orange","purple"),cex=0.8); legend("topright",legendText,cex=0.6)
			
			# overall survival
			sfit4<-survfit(Surv(timeToDeathInterval,  shortDeathEvent) ~ 1, data=SurvivalData)
			plotType<-paste("Overall survival ",sep="")
			shortPlotTitle <- paste(plotType,"\n, time ",(round(shortCensorPeriodStartDay)/DaySeconds)," to ",(round(max(SurvivalData$timeToDeathInterval))/DaySeconds)," days",sep="")
			plot(sfit4,mark.time=F,lty=1:6,conf.int=F,col=c("black","red","blue","green","orange","purple"),main=shortPlotTitle,xlim=c(shortCensorPeriodStartDay,round(max(SurvivalData$timeToDeathInterval))),lwd=3,ylim=c(0.5,1))
			
			

		dev.off()
		
# devloping a death stratification score using age, median, LOS and IQR. scoring weight is given by the odds ratio developed from the training set, and then applied to the test set.

		regressionTestSetTrainTest<-data.frame(regressionTestSet$deadAt1yearFactor, regressionTestSet$age, regressionTestSet$medianGlu, regressionTestSet$admissionDurationDays, regressionTestSet$IQR)
		colnames(regressionTestSetTrainTest)<-c("deadAt1yearFactor","age","medianGlu","admissionDurationDays","IQR")
	
			train <- regressionTestSetTrainTest[1:(nrow(regressionTestSetTrainTest)-1000),]
			test  <- regressionTestSetTrainTest[((nrow(regressionTestSetTrainTest)-1000)+1):nrow(regressionTestSetTrainTest),]
		
			model <- glm(deadAt1yearFactor ~.,family=binomial(link='logit'),data=regressionTestSetTrainTest)
			#
			summary(model)
		
			# to generate odds ratios per unit change in each parameter
			ORs<-exp(cbind(OR = coef(model),confint(model)))
			ageWt<-ORs[2,1]
			medianWt<-ORs[3,1]
			LOSwt<-ORs[4,1]
			IQRwt<-ORs[5,1]
		
			trainScore<-((test$age*ageWt)/(max(test$age)*ageWt))+((test$medianGlu*medianWt)/(max(test$medianGlu)*medianWt))+((test$admissionDurationDays*LOSwt)/(max(test$admissionDurationDays)*LOSwt))+((test$IQR*IQRwt)/(max(test$IQR)*IQRwt))
		
			plotFrame<-regressionFrame[((nrow(regressionFrame)-1000)+1):nrow(regressionFrame),]
			plotFrame$trainScore<-trainScore
			
		SurvivalData<-plotFrame

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


		sfit4<-survfit(Surv(timeToDeathInterval,  shortDeathEvent) ~ (trainScore>median(trainScore)), data=SurvivalData)
		plotType<-paste("IQR + medianGlu. upper/lower 50% ",sep="")
		shortPlotTitle <- paste(plotType,"\n, time ",(round(shortCensorPeriodStartDay)/DaySeconds)," to ",(round(max(SurvivalData$timeToDeathInterval))/DaySeconds)," days",sep="")
		plot(sfit4,mark.time=F,lty=1:6,conf.int=F,col=c("black","red","blue","green","orange","purple"),main=shortPlotTitle,xlim=c(shortCensorPeriodStartDay,round(max(SurvivalData$timeToDeathInterval))),lwd=3)
		sfit4.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ (trainScore>median(trainScore)), data=SurvivalData)
			pVal <- summary(sfit4.coxph)$coef[,5]; HR <- round(exp(coef(sfit4.coxph)),2)
			legendText <- paste("p = ",pVal," | HR = ",HR,sep="")
		summarySurvfit <- summary(sfit4); legendNames <- row.names(summarySurvfit$table)
		legend("bottomleft",c(legendNames),lty=1:6,col=c("black","red","blue","green","orange","purple"),cex=0.8); legend("topright",legendText,cex=0.6)
		
		sfit4<-survfit(Surv(timeToDeathInterval,  shortDeathEvent) ~ (trainScore>quantile(trainScore)[4]), data=SurvivalData)
		plotType<-paste("IQR + medianGlu. upper/lower 50% ",sep="")
		shortPlotTitle <- paste(plotType,"\n, time ",(round(shortCensorPeriodStartDay)/DaySeconds)," to ",(round(max(SurvivalData$timeToDeathInterval))/DaySeconds)," days",sep="")
		plot(sfit4,mark.time=F,lty=1:6,conf.int=F,col=c("black","red","blue","green","orange","purple"),main=shortPlotTitle,xlim=c(shortCensorPeriodStartDay,round(max(SurvivalData$timeToDeathInterval))),lwd=3)
		sfit4.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ (trainScore>quantile(trainScore)[4]), data=SurvivalData)
			pVal <- summary(sfit4.coxph)$coef[,5]; HR <- round(exp(coef(sfit4.coxph)),2)
			legendText <- paste("p = ",pVal," | HR = ",HR,sep="")
		summarySurvfit <- summary(sfit4); legendNames <- row.names(summarySurvfit$table)
		legend("bottomleft",c(legendNames),lty=1:6,col=c("black","red","blue","green","orange","purple"),cex=0.8); legend("topright",legendText,cex=0.6)
		
	
### plot survival by ratio of eAG/hba1c
SurvivalData<-admissionsWithDrugData

	DaySeconds<-(60*60*24)
	shortCensorPeriodStartDay  <- DaySeconds*0
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

SurvivalData_eAG_HbA1c<-merge(reportingDF,SurvivalData,by.x="ID",by.y="ID")


sfit4<-survfit(Surv(timeToDeathInterval,  shortDeathEvent) ~ (SurvivalData_eAG_HbA1c$ratio<quantile(SurvivalData_eAG_HbA1c$ratio)[2]), data=SurvivalData_eAG_HbA1c)
plotType<-paste("IQR + medianGlu. upper/lower 50% ",sep="")
shortPlotTitle <- paste(plotType,"\n, time ",(round(shortCensorPeriodStartDay)/DaySeconds)," to ",(round(max(SurvivalData$timeToDeathInterval))/DaySeconds)," days",sep="")
plot(sfit4,mark.time=F,lty=1:6,conf.int=F,col=c("black","red","blue","green","orange","purple"),main=shortPlotTitle,xlim=c(shortCensorPeriodStartDay,round(max(SurvivalData$timeToDeathInterval))),lwd=3)
sfit4.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ (trainScore>median(trainScore)), data=SurvivalData)
	pVal <- summary(sfit4.coxph)$coef[,5]; HR <- round(exp(coef(sfit4.coxph)),2)
	legendText <- paste("p = ",pVal," | HR = ",HR,sep="")
summarySurvfit <- summary(sfit4); legendNames <- row.names(summarySurvfit$table)
legend("bottomleft",c(legendNames),lty=1:6,col=c("black","red","blue","green","orange","purple"),cex=0.8); legend("topright",legendText,cex=0.6)




		#####
		# IQR 50% cut
		cases<-subset(admissionsWithDrugData_testSet, IQR>=12)
		controlPool<-subset(admissionsWithDrugData_testSet, IQR<12)
		controlSelection(cases,controlPool,paste("TEST",sep=""),1)
		
		# median 50% cut - matched for IQR and hypo episodes
		cases<-subset(admissionsWithDrugData_testSet, medianGlu>=quantile(admissionsWithDrugData_testSet$medianGlu)[4])
		controlPool<-subset(admissionsWithDrugData_testSet, medianGlu>=quantile(admissionsWithDrugData_testSet$medianGlu)[4])
		controlSelection(cases,controlPool,paste("mediabTcutoff:",quantile(admissionsWithDrugData_testSet$medianGlu)[3],sep=""),5)
		
		