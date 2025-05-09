% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/receptor.R
\name{prepare_receptor}
\alias{prepare_receptor}
\title{Prepare receptor PDB files for molecular docking}
\usage{
prepare_receptor(
  pdb_files,
  output_dir = "./",
  python_path = NULL,
  prepare_script = NULL,
  add_opts = "-A hydrogens"
)
}
\arguments{
\item{pdb_files}{Character vector of paths to PDB files}

\item{output_dir}{Character, output directory for prepared files, default "./"}

\item{python_path}{Character, path to Python executable, default "python"}

\item{prepare_script}{Character, path to prepare_receptor4.py script}

\item{add_opts}{Character, additional options to pass to prepare_receptor4.py, default "-A hydrogens -e True"}
}
\value{
Character vector of paths to successfully prepared PDBQT files
}
\description{
This function processes PDB files to prepare them for molecular docking by
converting them to PDBQT format using AutoDock Tools scripts.
}
\examples{
\dontrun{
pdbqt_files <- prepare_receptor(pdb_files,
                                output_dir = "prepared",
                                prepare_script = "path/to/prepare_receptor4.py")
}

}
