## collected from pathview package
.printf <- function(...) message(noquote(sprintf(...)))

.cleanUrl <- function(url)
{
    url <- gsub(" ", "%20", url, fixed=TRUE)
    url <- gsub("#", "%23", url, fixed=TRUE)
    url <- gsub(":", "%3a", url, fixed=TRUE)
    sub("http(s)*%3a//", "http\\1://", url)
}
##' @importFrom httr GET
##' @importFrom httr http_status
##' @importFrom httr content
##' @importFrom httr stop_for_status
.getUrl <- function(url, parser, ...)
{
    url <- .cleanUrl(url)
    debug <- getOption("KEGGREST_DEBUG", FALSE)
    if (debug)
        .printf("url == %s", url)
    response <- GET(url)
    stop_for_status(response)
    content <- .strip(content(response, "text"))
    if (nchar(content) == 0)
        return(character(0))
    do.call(parser, list(content, ...))
}

.strip <- function(str)
{
    gsub("^\\s+|\\s+$", "", str)
}


.textParser <- function(txt)
{
    txt
}

