##' @title Blast search against kegg database
##' @name seq2fun
##' @description seq2fun provide function to annotate sequences
##' @importFrom magrittr %>%
##' @importFrom dplyr select_
##' @importFrom dplyr filter_
##' @importFrom dplyr left_join
##' @importFrom utils read.table
##' @importFrom dplyr distinct_
##' @param query query sequences
##' @param dbfile database
##' @param type blastx, blastn or tblast
##' @param evalue e-value
##' @param num_threads num_threads
##' @param anno kegg annotation file
##' @author Kai Guo
##' @export
seq2fun <- function(query, dbfile, anno, type = 'blastp', evalue = 1e-10, num_threads = 1){
    outfmt <- 10
    tmpdir <- tempdir()
    tmp <- basename(tempfile(tmpdir = tmpdir))
    seq_in <- paste(tmpdir, "/", tmp, ".fasta", sep="")
    blastr <- paste(tmpdir, "/", tmp, ".txt", sep="")
    writeXStringSet(query, seq_in, format = "fasta", append = F)
    tool <- .find.tools(tool = type)
    params <- paste("-query", seq_in, "-db", dbfile, "-evalue", evalue,
            "-outfmt", outfmt, "-num_threads", num_threads,
            "-out", blastr, sep=" ")
    system(paste(tool, params, sep=" "))
    res <- read.table(blastr, sep=",", quote="", header=F)
    colnames(res)[1:2] <- c("Query",  "Subject")
    res <- res%>%distinct_(~Query,~Subject)
    res <- left_join(res, anno, by = c("Subject" = "keggid"))
    res <- res%>%select_(~Query, ~pathway, ~annotation)%>%distinct_(~Query,
            ~pathway, ~annotation)
    return(res)
}
##' make blast database
##' @param file sequence file
##' @param dbtype character "nucl" or "prot"
##' @param args other arguments for makeblastdb
##' @author Kai Guo
##' @export
makeblastdb <- function(file, dbtype = "prot", args="") {
    system(paste(.find.tools("makeblastdb"), "-in", file,
                 "-dbtype", dbtype, args))
}

