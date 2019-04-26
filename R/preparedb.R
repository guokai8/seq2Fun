##' @title Prepare pathway sequences database from KEGG
##' @importFrom KEGGREST keggGet
##' @importFrom KEGGREST keggList
##' @importFrom Biostrings writeXStringSet
##' @param species species name for KEGG parse
##' @param seqtype sequence type: AA or DNA (default: AA)
##' @param savedb save the kegg database for future using.
##' @param path path for saving kegg database (savedb = TRUE,default="./")
##' @return list with sequences and annotation for pathway
##' @author Kai Guo
##' @export
preparedb <- function(species = "human", seqtype = "AA", savedb = TRUE,
                    path = NULL){
    if(is.null(path)){
        path=getwd()
    }
    if(seqtype == "AA"){
        seqquery = "aaseq"
    }else{
        seqquery = "ntseq"
    }
    specode <- .get.spe.code(species = species)
    cat("Find ",species, "with kegg code: ",specode,"\n")
    keggentry <- get.kegg.entry(specode = specode)
    if(file.exists(paste(path,"/",specode, ".fasta", sep=""))){
        cat("Find ",paste(path,"/",specode, ".fasta", sep="")," in ",path,"\n")
        if(seqtype == "AA"){
            query=readAAStringSet(paste(path, "/", specode, ".fasta", sep=""))
        }else{
            query=readDNAStringSet(paste(path, "/", specode, ".fasta", sep=""))
        }
    }else{
    nr <- nrow(keggentry)
    nn <- ceiling(nr/10 + 1)
    lhs<- .vsplit(keggentry$keggid,nn)
    cat("Downloading sequences from KEGG database ......\n")
    query <- lapply(lhs, function(x) keggGetm(x, seqquery))
    query <- do.call(c,query)
    cat("All sequences were downloaded\n")
    if(isTRUE(savedb)){
        writeXStringSet(query,paste(path,"/",specode,".fasta",sep=""),format="fasta",append = F)
    }
    }
    out<-list(db = query, anno = keggentry, species = specode, type = seqquery)
    return(out)
}
##' find the species code for kegg
##' @param species character, either the kegg code, scientific name or the
##' common name of the target species. Default species="hsa", it is
##' equivalent to use either "Homo sapiens" (scientific name) or
##' "human" (common name).
.get.spe.code <- function(species="hsa"){
    species.data = get.kegg.code(species)
    species = species.data[1]
    return(species)
}
##' collect all the pathway name of species from KEGG
##' @importFrom KEGGREST keggLink
##' @importFrom KEGGREST keggList
##' @importFrom dplyr left_join
##' @param specode KEGG species code
##' @return a data frame with pathways information
##' @export
##' @author Kai Guo
get.kegg.entry <- function(specode){
    query <- keggLink(specode, "pathway")
    query <- data.frame(keggid = as.vector(query),pathway = sub('path:','',
            names(query)))
    pathways <- keggList('pathway', specode)
    pathway <- sub('path:', '', names(pathways))
    annot <- data.frame(cbind(pathways))
    annot$pathway <- pathway
    colnames(annot)[1] <- 'annotation'
    res <- left_join(query,annot,by=c('pathway'='pathway'))
    return(res)
}

##' split vector into chunk with provide length
##' @param v vector
##' @param n number of chunks
.vsplit <- function(v, n) {
    l = length(v)
    r = l/n
    return(lapply(1:n, function(i) {
        s = max(1, round(r*(i-1))+1)
        e = min(l, round(r*i))
        return(v[s:e])
    }))
}
##' Function modified from KEGGREST
##' @importFrom Biostrings readAAStringSet
##' @importFrom Biostrings readDNAStringSet
##' @param dbentries One or more (up to a maximum of 10) KEGG identifiers.
##' @param option Optional. Option governing the format of the output.
##' aaseq is an amino acid sequence, ntseq is a nucleotide sequence. image
##' returns an object which can be written to a PNG file, kgml returns a
##' KGML document.
keggGetm <- function(dbentries,
                    option=c("aaseq", "ntseq", "mol", "kcf", "image", "kgml"))
{
    if (length(dbentries) > 10)
        warning(paste("More than 10 inputs supplied, only the first",
                      "10 results will be returned."))
    dbentries <- paste(dbentries, collapse="+")
    url <- sprintf("%s/get/%s", "http://rest.kegg.jp", dbentries)
    if (!missing(option))
    {
        url <- sprintf("%s/%s", url, option)
        if (option %in% c("aaseq", "ntseq"))
        {
            t <- tempfile()
            tryCatch(cat(.getUrl(url, .textParser), file=t),error = function(cond){return("")})
            if (option == "aaseq")
                return(readAAStringSet(t))
            else if (option == "ntseq")
                return(readDNAStringSet(t))
        }
    }
}
##' extract kegg code
##' modified from pathview kegg.species.code
##' @param species
get.kegg.code <- function (species = "hsa")
{
    nspec = length(species)
    if (!exists("korg"))
        data(korg)
    ridx = match(species, korg[, 1:5])%%nrow(korg)
    nai = is.na(ridx)
    if (sum(nai) > 0) {
        si = try(load(url("https://pathview.uncc.edu/data/korg.1.rda")))
        if (class(si) != "try-error") {
            ridx.1 = match(species, korg.1[, 1:5])%%nrow(korg.1)
            nai.1 = is.na(ridx.1)
            if (sum(nai.1) < sum(nai)) {
                korg = korg.1
                ridx = ridx.1
                nai = nai.1
            }
        }
        if (sum(nai) > 0) {
            na.msg = sprintf("Unknown species '%s'!", paste(species[nai],
                                                            sep = "", collapse = "', '"))
            message("Note: ", na.msg)
        }
        if (sum(nai) == nspec) {
            stop.msg = "All species are invalid!"
            stop(stop.msg)
        }
    }
    if (any(ridx[!nai] == 0))
        ridx[!nai & ridx == 0] = nrow(korg)
    species.info = korg[ridx, 3]
    return(species.info)
}

