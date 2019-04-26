##' find excute blast or diamond
##' @param tool tools: blastp, blastx, blastn
##' @export
##' @author Kai Guo
.find.tools<-function(tool="blastp"){
    path <- Sys.which(tool)
    if(all(path=="")) {
        stop("No ",tool," not found! Please install ",tool,
        " manually or use get.tools to install ")
        return(character(0))
    }
    path[which(path!="")[1]]
}

##' @title Print the help function
##' @name blast_help
##' @param tool blast tools: blastn, blastp, blastx
##' @examples
##' blast_help()
##' @author Kai Guo
##' @export
blast_help <- function(tool = "blastx") {
    system(paste(.find.tools(c(tool)), "-help"))
}
##' @importFrom KEGGREST keggList
##' @title list organism names
##' @export
##' @author Kai Guo
listspecies <- function(){
    org <- keggList("organism")
    return(org)
}


