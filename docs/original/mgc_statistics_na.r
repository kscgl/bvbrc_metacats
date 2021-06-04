temp1=c(as.matrix(msaTable))
temp2=replace(temp1, which(temp1=="TRUE"),"T")
temp3=matrix(temp2,ncol=ncol(msaTable))
msaTableNew=as.data.frame(temp3)
n=ncol(msaTableNew)-2
names(msaTableNew)=c("mgcId","group",paste("Site",1:n, sep=""))
for(pos in 1:n) {
    contable = table(msaTableNew[[paste("Site",pos,sep="")]],msaTableNew[["group"]])
    if(any(rownames(contable)=='#')) {contable=contable[-(which(rownames(contable)=='#')),,drop=FALSE]}
    if(nrow(contable)==0) {next}
    if(nrow(contable)==1) {
        resultCST=paste(resultCST,paste("Site", pos, sep=""), "\tNA\t1\tNA\tN\tNA\n")
    } else {
        if(any(contable[contable[,]>0] < 5)) {flagsparse = "Y"} else {flagsparse = "N"}
        residueDiv = residueDiversity(contable)
        rname=rownames(contable)
        ACGT=c("A","C","G","T")
        matresult=match(ACGT, rname)
        A=rep(0,ncol(contable))
        C=rep(0,ncol(contable))
        G=rep(0,ncol(contable))
        T=rep(0,ncol(contable))
        if(is.na(matresult[1])) {contable=rbind(contable,A)}
        if(is.na(matresult[2])) {contable=rbind(contable,C)}
        if(is.na(matresult[3])) {contable=rbind(contable,G)}
        if(is.na(matresult[4])) {contable=rbind(contable,T)}
        fit=chisq.test(contable + 0.001)
        if(fit$p.value < pvalcutoff) {
            sigpvals <- c(sigpvals, fit$p.value)
            positions <- c(positions, pos)
        }
        resultCST=paste(resultCST, paste("Site", pos, sep=""), "\t", fit$statistic, "\t", fit$p.value, "\t", fit$parameter, "\t", flagsparse, "\t", residueDiv, "\n")
    }
    mcMat = mc(contable + 0.001)
    contRowCount = nrow(mcMat)
    for(j in 1:(contRowCount - 1)) {
        for(k in (j + 1):contRowCount) {
            if(j < k) {
                resultMC=paste(resultMC, paste("Site", pos, sep=""), "\t", mcMat[j,k], "\t", j, "\t", k, "\n")
            }
        }
    }
}