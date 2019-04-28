.onLoad <- function(libname, pkgname) {
  pnames=rownames(installed.packages())
  if("seq2Fun" %in% pnames){
    data(korg, package ="pathview")
  }
disclaimer="##############################################################################\n
seq2Fun is an open source software package. The seq2Fun downloads and uses KEGG data.
Non-academic uses may require a KEGG license agreement (details at http://www.kegg.jp/kegg/legal.html).\n
##############################################################################\n\n"
packageStartupMessage(wordwrap(disclaimer, 80))
}
