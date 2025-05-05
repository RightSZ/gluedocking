split_alt <- function(inputs,
                      out_dir = "split_alt",
                      python_path = NULL,
                      script_path = NULL,
                      keep_label = "A") {

  if (is.null(python_path)) {
    python_path <- Sys.getenv("GLUEDOCK_PYTHON_PATH", unset = NA)
    if (is.na(python_path)) {
      stop("Python path not provided and GLUEDOCK_PYTHON_PATH environment variable not set")
    }
  }

  if (is.null(script_path)) {
    script_path <- Sys.getenv("GLUEDOCK_PREPARE_SPLIT_ALT", unset = NA)
    if (is.na(script_path)) {
      stop("Script path not provided and GLUEDOCK_PREPARE_SPLIT_ALT environment variable not set")
    }
  }

  if (!file.exists(script_path)) stop("Split script not found: ", script_path)
  if (!nzchar(Sys.which(python_path)) && !file.exists(python_path)) {
    stop("Python executable not found: ", python_path)
  }
  python_path <- normalizePath(python_path, winslash = "/", mustWork = TRUE)
  script_path <- normalizePath(script_path, winslash = "/", mustWork = TRUE)

  pdbs <- collect_pdb_files(inputs)
  if (!length(pdbs)) stop("No PDB files found in inputs.")

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
    message("Created output directory: ", out_dir)
  }

  kept <- process_pdb_files(pdbs, out_dir, python_path, script_path, keep_label)

  invisible(kept)
}

#' Collect PDB files from input paths
#'
#' @param inputs Character vector of file paths or directories
#' @return Character vector of PDB file paths
collect_pdb_files <- function(inputs) {
  unlist(lapply(inputs, function(x) {
    if (dir.exists(x)) {
      list.files(x, "\\.pdb$", full.names = TRUE)
    } else if (file.exists(x)) {
      if (grepl("\\.pdb$", x, ignore.case = TRUE)) {
        x
      } else {
        warning("Skipping non-PDB file: ", x)
        NULL
      }
    } else {
      warning("Skipping invalid path: ", x)
      NULL
    }
  }), use.names = FALSE)
}

#' Check if PDB file contains alternate locations
#'
#' @param pdb_file Path to PDB file
#' @return Logical indicating if alternate locations exist
has_alternate_locations <- function(pdb_file) {
  lines <- readLines(pdb_file)
  atom_lines <- lines[grepl("^(ATOM|HETATM)", lines)]
  altlocs <- substring(atom_lines, 17, 17)

  tags <- setdiff(unique(altlocs), " ")
  length(tags) > 0
}

#' Process PDB files for alternate location splitting
#'
#' @param pdbs Character vector of PDB file paths
#' @param out_dir Output directory
#' @param python_path Path to Python executable
#' @param script_path Path to splitting script
#' @param keep_label Label to keep (default "A")
#' @return Character vector of processed file paths
process_pdb_files <- function(pdbs, out_dir, python_path, script_path, keep_label) {
  kept <- character(0)

  for (pdb in pdbs) {
    stem <- tools::file_path_sans_ext(basename(pdb))
    message("Processing: ", basename(pdb))

    out_file <- file.path(out_dir, paste0(stem, "_", keep_label, ".pdb"))

    if (!has_alternate_locations(pdb)) {
      file.copy(pdb, out_file, overwrite = TRUE)
      message("  No alternate locations found; copied as: ", basename(out_file))
      kept <- c(kept, out_file)
      next
    }

    kept <- c(kept, split_pdb_with_script(pdb, out_dir, stem, python_path, script_path, out_file, keep_label))
  }

  return(kept)
}

#' Split PDB file with external Python script
#'
#' @param pdb_file Path to PDB file
#' @param out_dir Output directory
#' @param stem File stem name
#' @param python_path Path to Python executable
#' @param script_path Path to splitting script
#' @param out_file Expected output file path
#' @param keep_label Label to keep
#' @return Character vector with path to processed file if successful
split_pdb_with_script <- function(pdb_file, out_dir, stem, python_path, script_path, out_file, keep_label) {
  prefix <- file.path(out_dir, stem)
  prefix_norm <- normalizePath(prefix, winslash = "/", mustWork = FALSE)
  args <- c("-r", normalizePath(pdb_file, winslash = "/", mustWork = TRUE),
            "-o", prefix_norm)

  message("  Splitting alternate locations via script...")
  tryCatch({
    res <- system2(python_path, args = c(script_path, args),
                  stdout = TRUE, stderr = TRUE)
    cat(res, sep = "\n")

    if (file.exists(out_file)) {
      message("  Found split file: ", basename(out_file))
      return(out_file)
    } else {
      warning("  Expected file not found: ", basename(out_file))
      return(character(0))
    }
  }, error = function(e) {
    warning("  Error running script: ", e$message)
    return(character(0))
  })
}
