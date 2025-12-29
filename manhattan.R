setwd("F:/plot/rice/")


library(data.table)
library(qqman)
#tiff("highresolution.tiff",width=22000, height=12000, units= "px", pointsize =60,res=600)

rm(list=ls())


bqz<-1.3###########-log10p????
bqzq=1.3#########qq??ǩ????
bqzl=0.8###########lod??ǩ????
zv=0.4#########qqtu ????ֵ????

cu=1.0
kel<--0.03
kes<-0.5
biaoq<-0.7

########################
tiff("manhattan7.tiff",width=238, height=102,units = "mm",res=2000)

xia<-3.0
zuo<-2.6
shang<-1.5
you<-1.0

par(mar=c(xia,zuo,shang,you)+0.1)

layout(matrix(c(1,2,3,4,5,6),3,2,byrow = T),widths=c(5.3,1),heights=c(1,1,1))
#####################################
# #
data<-fread("1_intermediate result.csv",header=T)
#data<-data[which(data[,3]=="mrMLM"),]
#data<-data[which(data[,5]!=11),]
data<-as.matrix(data)
data4<-as.matrix(seq(1:1011601))

data8<-as.numeric(data[,8])
pdata8<-10^-data8
locsub<-which(pdata8==0)
pmin<-min(pdata8[pdata8!=0])
subvalue<-10^(1.1*log10(pmin))
pdata8[locsub]<-subvalue
p8change<-as.matrix(pdata8)

data5<-as.matrix(data[,5])
data6<-as.matrix(data[,6])


manresult<-cbind(data5,data6,p8change,data4)


manresult<-apply(manresult,2,as.numeric)

colnames(manresult)<-c("Chromosome","BPnumber","P-value","SNPname")


manresult<-as.data.frame(manresult)

# top_snps = c('rs13895')
# surrounding_snps = list(as.character(manresult$SNPname[1]))

# manhattan(manresult,chr = "Chromosome",bp ="BPnumber",p ="P-value",snp="SNPname",col=c("red","black"),suggestiveline=FALSE,
#           genomewideline = 6.5,highlight = snpOfInterest[,4],ylab=expression('-log'[10]*'(P-value)'),ylim=c(0,20),family="serif")

r1<-as.matrix(fread("all-gene-shu.csv",header=T))

#rr1<-as.matrix(fread("lingwai.csv",header=F))
manresulty2<-manresult

loc1<-numeric()
for(i in 1:nrow(r1)){
  zhi<-as.numeric(r1[i,2])
  loc1[i]<-which(manresulty2[,2]==zhi)
}

snpOfInterest<-manresult[loc1,4]

x<-manresult

chr = "Chromosome"
bp ="BPnumber"
p ="P-value"
snp="SNPname"
col=c("lightgreen","lightskyblue")
#col=c("gray","gray60")

suggestiveline=FALSE
genomewideline = 0

#highlight = snpOfInterest

ylab=expression('-log'[10]*'(P-value)')
#ylim=c(0,20)
#family="serif"

logp=TRUE
annotatePval = NULL
annotateTop = TRUE

highlight<-snpOfInterest


CHR=BP=P=index=NULL

d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]])
if (!is.null(x[[snp]])) d=transform(d, SNP=x[[snp]])
d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
d <- d[order(d$CHR, d$BP), ]
if (logp) {
  d$logp <- -log10(d$P)
} else {
  d$logp <- d$P
}
d$pos=NA
d$index=NA
ind = 0
for (i in unique(d$CHR)){
  ind = ind + 1
  d[d$CHR==i,]$index = ind
}


nchr = length(unique(d$CHR))
if (nchr==1) { ## For a single chromosome
  ## Uncomment the next two linex to plot single chr results in Mb
  #options(scipen=999)
  #d$pos=d$BP/1e6
  d$pos=d$BP
  ticks=floor(length(d$pos))/2+1
  xlabel = paste('Chromosome',unique(d$CHR),'position')
  labs = ticks
} else { ## For multiple chromosomes
  lastbase=0
  ticks=NULL
  ticks1=NULL
  for (i in unique(d$index)) {
    if (i==1) {
      d[d$index==i, ]$pos=d[d$index==i, ]$BP
    } else {
      lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
      d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
    }
    # Old way: assumes SNPs evenly distributed
    # ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
    # New way: doesn't make that assumption
    ticks1 = c(ticks1, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)
    ticks = c(ticks,c(min(d[d$index == i,]$pos),max(d[d$index == i,]$pos)))
  }
  ticks=ticks-800000
  ticks[length(ticks)]=ticks[length(ticks)]+800000*2
  xlabel = ''
  #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
  labs <- unique(d$CHR)
}

xmax = ceiling(max(d$pos) * 1.03)
xmin = floor(max(d$pos) * -0.03)


def_args <- list(xaxt='n',yaxt="n",bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                 xlim=c(xmin,xmax), ylim=c(0,50),
                 xlab=xlabel,ylab="",mgp=c(1.3,0.3,0),cex.lab=biaoq)

dotargs <- list(NULL)
do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))



if (!is.null(chrlabs)) {
  if (is.character(chrlabs)) {
    if (length(chrlabs)==length(labs)) {
      labs <- chrlabs
    } else {
      warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
    }
  } else {
    warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
  }
}

# Add an axis.
if (nchr==1) { #If single chromosome, ticks and labels automatic.
  axis(1, ...)
} else { # if multiple chrs, use the ticks and labels you created above.
  axis(1, at=ticks, labels=rep("",24),lwd=cu,tck=kel,mgp=c(3,0.1,0.5),cex.axis=kes+0.2)
}

mtext("Chromosomes in rice",side=1,line=bqzq+0.6,cex=biaoq-0.1,font=2)


axis(2, at=seq(0,50,10),lwd=cu,tck=kel,mgp=c(bqz,0.5,0),cex.axis=biaoq-0.1,las=1)
#mtext(expression(-log[10](P)),side=2,line=bqz,cex=biaoq-0.2,font=1)

mtext(expression(paste(-Log[10]," ",italic("P")," ","value",sep="")),side=2,line=bqzq,cex=biaoq-0.2,font=1)


# Create a vector of alternatiting colors
col=rep(col, max(d$CHR))
#############################
#############################
#Add points to the plot
if (nchr==1) {
  with(d, points(pos, logp, pch=20, col=col[1], ...))
} else {
  # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
  icol=1
  for (i in unique(d$index)) {
    with(d[d$index==unique(d$index)[i], ], points(pos, logp, col=col[icol], pch=20,cex=0.7))
    icol=icol+1
  }
}
#############################
#############################


# Add suggestive and genomewide lines
# if (suggestiveline) abline(h=suggestiveline, col="black")
# if (genomewideline) abline(h=genomewideline, col="red")

# Highlight snps from a character vector

d.highlight=d[which(d$SNP %in% highlight), ]
#d.highlight<-manresult[loc1,]

LOD<-as.numeric(r1[,3])

# pp<-numeric(27)
# pp[9]<--log10(4.6e-07)
# pp[10]<--log10(1.1e-06)
# pp[11]<--log10(8.5e-07)
# pp[13]<--log10(9.4e-07)



d.highlight<-as.data.frame(cbind(d.highlight,LOD))


#
# with(d.highlight, points(pos[c(9,10,11,13)], pp[c(9,10,11,13)], col="gray50", pch=20,type="h"))
# with(d.highlight, points(pos[c(9,10,11,13)], pp[c(9,10,11,13)], col="gray50", pch=20))
#


#with(d.highlight, points(pos, LOD, col="red", pch=20))

par(new=T)

def_args <- list(xaxt='n', yaxt='n',bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                 xlim=c(xmin,xmax), ylim=c(0,16),xlab="",ylab="")

dotargs <- list(NULL)
do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))

with(d.highlight, points(pos[c(1,2+1,6+1,13+1,14+1,15+1,17+2,19+2,20+2,22+2,23+2)], LOD[c(1,2+1,6+1,13+1,14+1,15+1,17+2,19+2,20+2,22+2,23+2)], col="red", pch=20,type="h",lwd=0.25))

with(d.highlight, points(pos[c(16+2,24+2)], LOD[c(16+2,24+2)], col="black", pch=20,type="h",lwd=0.25))

with(d.highlight, points(pos[c(2,17)], LOD[c(2,17)], col="blue", pch=20,type="h",lwd=0.25))

axis(4,mgp=c(bqz,zv,0),at=seq(0,16,4),col="black",col.ticks="black",col.axis="black",lwd=cu,tck=kel,cex.axis=kes,las=1)


mtext("LOD score",side=4,line=bqzl,cex=biaoq-0.2,font=1,col="black")

with(d.highlight, points(pos[c(1,2+1,6+1,13+1,14+1,15+1,17+2,19+2,20+2,22+2,23+2)], LOD[c(1,2+1,6+1,13+1,14+1,15+1,17+2,19+2,20+2,22+2,23+2)], col="red", pch=20))

with(d.highlight, points(pos[c(16+2,24+2)], LOD[c(16+2,24+2)], col="black", pch=20))

with(d.highlight, points(pos[c(2,17)], LOD[c(2,17)], col="blue", pch=20))

abline(h=3,col="red",lty=5,lwd=0.15)

zis<-0.4

text(-18e6,18.5,"A  Rice",col="black",cex=biaoq+0.2,xpd=T,font=2)
text(ticks1,-3.5,c("1","2","3","4","5","6","7","8","9","10","11","12"),col="black",cex=biaoq-0.1,xpd=T)


#######################################
######################################

text(25e6,6,"OsOFP2",col="red",cex=zis,font=3)
text(40e6,4.8,"Os05g0182500",col="blue",cex=zis)
text(50e6,6.2,"FUWA",col="red",cex=zis,font=3)
# text(63e6,5.5,"SDG725",col="red",cex=zis,font=3)
# text(70e6,4.6,"OsNST1",col="red",cex=zis,font=3)
# text(75e6,6.0,"OsMKK4",col="red",cex=zis,font=3)
text(83e6,7.7,"BG1",col="red",cex=zis,font=3)
# text(86e6,8.0,"OsMADS1",col="red",cex=zis,font=3)
# text(92e6,6.0,"OsAPC6",col="red",cex=zis,font=3)
# text(94e6,4.8,expression(paste(italic("OspPLAIII"),italic(alpha))),col="red",cex=zis)
# text(96e6,13,"GS3",col="red",cex=zis,font=3)
# text(116e6,4.4,"OsNaPRT1",col="red",cex=zis,font=3)
# text(125e6,5.4,"qGwt4a",col="red",cex=zis,font=3)
text(147e6,4.3,"flo2",col="red",cex=zis,font=3)
text(150e6,7.5,"OsMKP1",col="red",cex=zis,font=3)
text(154e6,5.0,"GS5",col="red",cex=zis,font=3)
text(148e6,9,"Os05g0182500",col="blue",cex=zis)
text(157e6,13.9,"GW5",col="black",cex=zis,font=3)
text(180e6,7,"OsACS6",col="red",cex=zis,font=3)
#text(188e6,5.8,"qGw-6",col="red",cex=zis,font=3)
text(228e6,3.8,"OsBZR1",col="red",cex=zis,font=3)
text(236e6,5.0,"qGL7",col="red",cex=zis,font=3)
#text(244e6,4.6,"OsFIE2",col="red",cex=zis,font=3)
text(263e6,4.3,"WTG1",col="red",cex=zis,font=3)
text(270e6,5.0,"SLG",col="red",cex=zis,font=3)
text(292e6,5.2,"OsFD1",col="black",cex=zis,font=3)
#text(336e6,8.6,"IKU2",col="red",cex=zis,font=3)

legend(235e6,19.0,c("","","",""),
       col=c("red","black","blue","gray50"),pch=20,
       bty="n",y.intersp=0.8,x.intersp=0.2,pt.cex=0.8,cex=0.7,horiz=T,xpd=T)

text(269e6,17.2,"LOD score",col="black",cex=zis+0.5,font=1,xpd=T)

legend(294e6,19.0,c("",""),
       col=c("lightgreen","lightskyblue"),pch=20,
       bty="n",y.intersp=0.8,x.intersp=0.2,pt.cex=0.8,cex=0.7,horiz=T,xpd=T)

text(324.5e6,17.2,expression(paste(-Log[10]," ",italic("P")," ","value",sep="")),col="black",cex=zis+0.5,font=2,xpd=T)

dev.off()