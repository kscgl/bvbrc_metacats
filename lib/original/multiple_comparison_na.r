library(gregmisc)
mc=function(mat) {
    nrow=nrow(mat)
    ncol=ncol(mat)
    result.pvalue=matrix(rep(NA,ncol*ncol),ncol=ncol)
    for(ii in 1 : (ncol - 1)) {
        for(jj in (ii+1):ncol) {
            obs=rbind(mat[,ii],mat[,jj])
            temp1=rep(NA,nrow)
            temp2=rep(NA,nrow)
            for(ss in 1:nrow) {
                temp1[ss]=sum(mat[,ii])*sum(mat[ss,ii]+mat[ss,jj])/(sum(mat[,ii])+sum(mat[,jj]))
                temp2[ss]=sum(mat[,jj])*sum(mat[ss,ii]+mat[ss,jj])/(sum(mat[,ii])+sum(mat[,jj]))
            }
            exp=rbind(temp1,temp2)
            stat=sum((obs-exp)^2/exp)
            result.pvalue[jj,ii]=pchisq(stat, df=(nrow-1), lower.tail = FALSE)
            result.pvalue[ii,jj]=result.pvalue[jj,ii]
        }
    }
    diag(result.pvalue) = 1
    return(result.pvalue)
}