residueDiversity = function(contable) {
    retval <- c("")
    nresidue = nrow(contable)
    ngroup = ncol(contable)
    for(i in 1 : ngroup) {
        retval = paste(retval, "group", i, "(", sep="")
        for(j in 1 : nresidue) {
            residue = rownames(contable)[j]
            count = contable[rownames(contable)==residue,i]
            if(count > 0) {
                retval = paste(retval, count, "_", residue, ",_", sep="")
            }
        }
        retval = paste(retval, ")", sep="")
        if(i < ngroup) {
            retval = paste(retval, "|", sep="")
        }
    }
    retval = gsub(",_)", ")", retval, fixed=TRUE)
    return (retval)
}