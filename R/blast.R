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
##' @param db database
##' @param type blastx, blastn or tblast
##' @param evalue e-value
##' @param num_threads num_threads
##' @param path working directory
##' @param bidirectional use bidirectional blast
##' @param runblast use blast+ or diamond
##' @author Kai Guo
##' @export
seq2fun <- function(query, db, type = 'blastp', evalue = 1e-10,
                    num_threads = 1, path=NULL, bidirectional = FALSE,
                    runblast = TRUE){
    if(is.null(path)){
        path=getwd()
    }
    seqcode <- db$species
    filepath = paste(path,"/",seqcode,".fasta",sep="")
    outfmt <- 6
    tmpdir <- tempdir()
    tmp <- basename(tempfile(tmpdir = tmpdir))
    seq_in <- paste(tmpdir, "/", tmp, ".fasta", sep="")
    blastr <- paste(tmpdir, "/", tmp, ".txt", sep="")
    writeXStringSet(query, seq_in, format = "fasta", append = F)
    if(isTRUE(runblast)){
    runblast(seq_in, dbfile = filepath,tool = type, evalue = evalue,
            outfmt = outfmt,num_threads = num_threads, blastr = blastr)
    }else{
    rundiamond(seq_in, dbfile = filepath,tool = type, evalue = evalue,
            outfmt = outfmt,num_threads = num_threads, blastr = blastr)
    }
    res <- read.table(blastr, sep="\t", quote="", header=F)
    colnames(res)[1:2] <- c("Query",  "Subject")
    res <- res%>%distinct_(~Query, ~Subject)
    if(isTRUE(bidirectional)){
        .makeblastdb(filepath = seq_in, dbtype = "prot")
        blastrq <- paste(tmpdir, "/", tmp, "q.txt", sep="")
#       paramsq <- paste("-query", filepath, "-db", seq_in, "-evalue", evalue,
#                         "-outfmt", outfmt, "-num_threads", num_threads,
#                         "-out", blastrq, sep=" ")
#        system(paste(tool, paramsq, sep=" "))
        if(isTRUE(runblast)){
        runblast(filepath, dbfile = seq_in, tool = type, evalue = evalue,
                outfmt = outfmt,
                num_threads = num_threads, blastrq = blastrq)
        }else{
        rundiamond(filepath, dbfile = seq_in, tool = type, evalue = evalue,
                outfmt = outfmt,
                num_threads = num_threads, blastr = blastrq)
        }
        resq <- read.table(blastrq, sep="\t", quote="", header=F)
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
##' @param db KEGG db or sequences
##' @param seqcode speceis name or KEGG code (without spaces)
##' @param dbtype character "nucl" or "prot"
##' @param args other arguments for makeblastdb
##' @param path working directory
##' @param runblast use blast+ or diamond
##' @author Kai Guo
##' @export
makeblastdb <- function(db, seqcode = NULL, dbtype = "prot", args="", path=NULL,
        runblast = TRUE) {
    if(is.null(path)){
        path=getwd()
    }
    if(is.null(seqcode)){
        seqcode=db$species
    }
    if(is.list(db)){
        seqs = db$db
    }else{
        seqs = db
    }
    filepath = paste(path,"/",seqcode,".fasta",sep="")
    writeXStringSet(seqs, filepath, format="fasta",append=F)
    if(isTRUE(runblast)){
    .makeblastdb(filepath, dbtype = dbtype)
    }else{
    .makediamonddb(filepath = filepath)
    }
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

##' diamond makedb command
##' @param filepath path of fasta file
##' @param args other arguments for diamond makedb
##' @export
##' @author Kai Guo
.makediamonddb <- function(filepath, args = ""){
    system(paste(.find.tools("diamond"), "makedb --in",
        filepath, "-d", filepath))
}

##' run blast analysis
##' @param seq_in query sequences
##' @param dbfile database
##' @param tool blastx, blastn or blastp
##' @param evalue e-value
##' @param num_threads num_threads
##' @param outfmt blast out format
##' @param blastr blast output result file
##' @export
##' @author Kai Guo
runblast <- function(seq_in, dbfile, tool = 'blastp', evalue = 1e-10,
                    num_threads = 1, outfmt = 10, blastr){
    params <- paste("-query", seq_in, "-db", dbfile, "-evalue", evalue,
                    "-outfmt", outfmt, "-num_threads", num_threads,
                    "-out", blastr, sep=" ")
    tool <- .find.tools(tool = tool)
    system(paste(tool, params, "> /dev/null", sep=" "))
    cat("Finish blast ",tool," ......\n")
}

##' run dimond blast analysis
##' @param seq_in query sequences
##' @param dbfile database
##' @param tool blastx, blastn or blastp
##' @param evalue e-value
##' @param num_threads num_threads
##' @param outfmt blast out format
##' @param blastr blast output result file
##' @export
##' @author Kai Guo
rundiamond <- function(seq_in, dbfile, tool = 'blastp', evalue = 1e-10,
                        num_threads = 1, outfmt = 6, blastr){
    params <- paste("-q", seq_in, "-d", dbfile, "-e", evalue,
                    "-f", outfmt, "-p", num_threads,
                    "-o", blastr, "--sensitive",sep=" ")
    loc <- .find.tools(tool = "diamond")
    system(paste(loc, tool, params, "> /dev/null", sep=" "))
    cat("Finish dimaond ",tool," ......\n")
}
