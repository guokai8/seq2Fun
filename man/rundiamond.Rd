% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blast.R
\name{rundiamond}
\alias{rundiamond}
\title{run dimond blast analysis}
\usage{
rundiamond(seq_in, dbfile, tool = "blastp", evalue = 1e-10,
  num_threads = 1, outfmt = 6, blastr)
}
\arguments{
\item{seq_in}{query sequences}

\item{dbfile}{database}

\item{tool}{blastx, blastn or blastp}

\item{evalue}{e-value}

\item{num_threads}{num_threads}

\item{outfmt}{blast out format}

\item{blastr}{blast output result file}
}
\description{
run dimond blast analysis
}
\author{
Kai Guo
}
