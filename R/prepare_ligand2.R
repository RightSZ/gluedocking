
#' Prepare ligand files for molecular docking
#'
#' This function processes molecule files to prepare them for molecular docking by
#' converting them to PDBQT format using AutoDock Tools scripts.
#'
#' @param mol_files Character vector of paths to molecule files (SDF, MOL, etc.)
#' @param out_dir Character, output directory for prepared files, default "./"
#' @param python_path Character, path to Python executable, default "python"
#' @param prepare_script Character, path to prepare_ligand4.py script
#' @param add_opts Character, additional options to pass to prepare_ligand4.py, default "-A hydrogens"
#'
#' @return Character vector of paths to successfully prepared PDBQT files
#' @importFrom tools file_path_sans_ext
#' @examples
#' \dontrun{
#' pdbqt_files <- prepare_ligand(sdf_files,
#'                               out_dir = "prepared",
#'                               prepare_script = "path/to/prepare_ligand4.py")
#' }
#'
#' @export
prepare_ligand2 <- function(mol_files, out_dir = "./",
                           python_path = NULL,
                           prepare_script = NULL,
                           add_opts = "-A hydrogens") {

  if (length(mol_files) == 0) {
    stop("No molecule files provided")
  }

  mol_files <- as.character(mol_files)

  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
  }

  if (is.null(prepare_script)) {
    prepare_script <- Sys.getenv("GLUEDOCK_PREPARE_LIGAND", unset = NA)
    if (is.na(prepare_script)) {
      stop("prepare_script parameter is required and GLUEDOCK_PREPARE_LIGAND environment variable not set. Please provide prepare_script parameter or set GLUEDOCK_PREPARE_LIGAND environment variable.")
    }
  }

  if (!file.exists(prepare_script)) {
    stop("Specified prepare_script does not exist: ", prepare_script)
  }

  original_wd <- getwd()
  on.exit(setwd(original_wd), add = TRUE)
  is_absolute_path <- function(path) {
    # On Unix/Mac, absolute paths start with "/"
    # On Windows, they start with a drive letter followed by ":" or with "\\"
    grepl("^(/|[A-Za-z]:|\\\\\\\\)", path)
  }

  process_one <- function(mol_file) {
    if (!file.exists(mol_file)) {
      warning(sprintf("Input molecule file does not exist: %s", mol_file))
      return(NULL)
    }

    mol_dir <- dirname(mol_file)
    mol_basename <- basename(mol_file)
    base_name <- tools::file_path_sans_ext(mol_basename)


    message(sprintf("Changing to directory: %s", mol_dir))
    setwd(mol_dir)


    if (is.null(out_dir)) {
      output_pdbqt <- paste0(base_name, ".pdbqt")
      abs_output_pdbqt <- file.path(mol_dir, output_pdbqt)
    } else {

      if (!is_absolute_path(out_dir)) {
        abs_out_dir <- file.path(original_wd, out_dir)
      } else {
        abs_out_dir <- out_dir
      }

      if (!dir.exists(abs_out_dir)) {
        dir.create(abs_out_dir, recursive = TRUE)
      }

      output_pdbqt <- paste0(base_name, ".pdbqt")
      abs_output_pdbqt <- file.path(abs_out_dir, output_pdbqt)
    }

    tool_path <- normalizePath(prepare_script, winslash = "/", mustWork = FALSE)

    if (is.null(python_path)) {
      python_path <- Sys.getenv("GLUEDOCK_PYTHON_PATH", unset = NA)
      if (is.na(python_path)) {
        stop("Python path not provided and GLUEDOCK_PYTHON_PATH environment variable not set. Please provide python_path parameter or set GLUEDOCK_PYTHON_PATH environment variable.")
      }
    }

    py_path <- normalizePath(python_path, winslash = "/", mustWork = FALSE)

    if (!file.exists(py_path)) {
      warning(sprintf("Python executable not found: %s", python_path))
      return(NULL)
    }

    # Use the local filename since we're in the file's directory
    args <- c(
      shQuote(tool_path),
      "-l", shQuote(mol_basename),
      "-o", shQuote(abs_output_pdbqt),
      add_opts
    )

    message("Running:\n  ", py_path, " ", paste(args, collapse = " "))

    res <- suppressWarnings(
      system2(
        command = py_path,
        args = args,
        stdout = TRUE,
        stderr = TRUE
      )
    )

    res_clean <- res[!grepl("swig/python detected a memory leak|SWIG detected a memory leak|memory leak of type 'BHtree'", res, ignore.case = TRUE)]

    if (length(res_clean)) message(paste(res_clean, collapse = "\n"))

    setwd(original_wd)

    if (!file.exists(abs_output_pdbqt)) {
      warning(sprintf("PDBQT not generated, please check output messages for file: %s", mol_file))
      return(NULL)
    }

    message(sprintf("Successfully prepared ligand: %s", abs_output_pdbqt))
    return(abs_output_pdbqt)
  }

  out_files <- lapply(mol_files, process_one)
  valid_files <- out_files[!sapply(out_files, is.null)]

  total <- length(mol_files)
  success <- length(valid_files)
  failed <- total - success

  if (failed > 0) {
    message(sprintf("Total: %d files, Success: %d, Failed: %d", total, success, failed))

    failed_files <- mol_files[sapply(out_files, is.null)]
    message("Failed files: ", paste(basename(failed_files), collapse = ", "))
  }

  return(unlist(valid_files))
}
