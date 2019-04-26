##' @title Blast search against kegg database
##' @name seq2fun
##' @description seq2fun provide function to annotate sequences
##' @importFrom magrittr %>%
##' @importFrom dplyr select_
##' @importFrom dplyr filter_
##' @importFrom dplyr left_join
##' @importFrom dplyr inner_join
##' @importFrom dplyr mutate
##' @importFrom utils read.table
##' @importFrom dplyr distinct_
##' @param query query sequences
##' @param dbfile database
##' @param type blastx, blastn or tblast
##' @param evalue e-value
##' @param num_threads num_threads
##' @param anno kegg annotation file
##' @param path working directory
##' @param bidirectional use bidirectional blast
##' @author Kai Guo
##' @export
seq2fun <- function(query, db, type = 'blastp', evalue = 1e-10,
                    num_threads = 1, path=NULL, bidirectional = FALSE){
    if(is.null(path)){
        path=getwd()
    }
    seqcode <- db$species
    filepath = paste(path,"/",seqcode,".fasta",sep="")
    outfmt <- 10
    tmpdir <- tempdir()
    tmp <- basename(tempfile(tmpdir = tmpdir))
    seq_in <- paste(tmpdir, "/", tmp, ".fasta", sep="")
    blastr <- paste(tmpdir, "/", tmp, ".txt", sep="")
    writeXStringSet(query, seq_in, format = "fasta", append = F)
    tool <- .find.tools(tool = type)
    params <- paste("-query", seq_in, "-db", filepath, "-evalue", evalue,
            "-outfmt", outfmt, "-num_threads", num_threads,
            "-out", blastr, sep=" ")
    system(paste(tool, params, sep=" "))
    res <- read.table(blastr, sep=",", quote="", header=F)
    colnames(res)[1:2] <- c("Query",  "Subject")
    res <- res%>%distinct_(~Query, ~Subject)
    if(isTRUE(bidirectional)){
        .makeblastdb(filepath = seq_in, dbtype = "prot")
        blastrq <- paste(tmpdir, "/", tmp, "q.txt", sep="")
        paramsq <- paste("-query", filepath, "-db", seq_in, "-evalue", evalue,
                         "-outfmt", outfmt, "-num_threads", num_threads,
                         "-out", blastrq, sep=" ")
        system(paste(tool, paramsq, sep=" "))
        resq <- read.table(blastrq, sep=",", quote="", header=F)
        colnames(resq)[1:2] <- c("Queryq", "Subjectq")
        resq <- resq%>%distinct_(~Queryq, ~Subjectq)%>%
            mutate(seq = paste(Subjectq, Queryq, sep="_"))
        res <- res%>%mutate(seq = paste(Query, Subject, sep="_"))%>%
            inner_join(resq, by = c("seq" = "seq"))%>%select_(~Query, ~Subject)
    }
    res <- left_join(res, db$anno, by = c("Subject" = "keggid"))
    res <- res%>%select_(~Query, ~pathway, ~annotation)%>%distinct_(~Query,
            ~pathway, ~annotation)
    return(res)
}
##' make blast database
##' @importFrom Biostrings writeXStringSet
##' @param db KEGG db
##' @param dbtype character "nucl" or "prot"
##' @param args other arguments for makeblastdb
##' @author Kai Guo
##' @export
makeblastdb <- function(db, dbtype = "prot", args="", path=NULL) {
    if(is.null(path)){
        path=getwd()
    }
    seqcode <- db$species
    filepath = paste(path,"/",seqcode,".fasta",sep="")
    writeXStringSet(db$db, filepath, format="fasta",append=F)
    .makeblastdb(filepath, dbtype = "prot")
}

##' makeblastdb command
##' @param filepath path of fasta file
##' @param dbtype character "nucl" or "prot"
##' @param args other arguments for makeblastdb
##' @export
##' @author Kai Guo
.makeblastdb <- function(filepath, dbtype = "prot", args = ""){
    system(paste(.find.tools("makeblastdb"), "-in", filepath,
                 "-dbtype", dbtype, args))
}
