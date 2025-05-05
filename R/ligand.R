#' Convert molecule file formats using OpenBabel
#'
#' This function converts molecular structure files between different formats
#' using the OpenBabel chemical toolbox. It supports converting a single file,
#' multiple files, or all files in a directory.
#'
#' @param input_file Character, path to the input molecule file, or a character vector of file paths,
#'                   or a directory containing molecule files to convert
#' @param input_format Character, format of the input file, default NULL (auto-detect from file extension)
#' @param output_file Character, path to the output file, default NULL (auto-generate based on input file name)
#' @param output_format Character, format for the output file, default NULL (auto-detect from output file extension)
#' @param output_dir Character, directory to save output files, default converted
#' @param pattern Character, file pattern to match when input_file is a directory, default NULL (all files)
#' @param recursive Logical, whether to search recursively in directories, default FALSE
#' @param obabel_path Character, path to OpenBabel executable, default NULL (uses GLUEDOCK_OBABEL_PATH environment variable)
#'
#' @return Character vector of paths to successfully converted files, or NULL if conversion failed
#' @importFrom tools file_path_sans_ext file_ext
#' @examples
#' \dontrun{
#' # Convert a single file
#' convert_molecule("ligand.sdf", output_format = "mol2")
#' }
#' @export
convert_molecule <- function(input_file, input_format = NULL,
                             output_file = NULL, output_format = NULL,
                             output_dir = "converted", pattern = NULL, recursive = FALSE,
                             obabel_path = NULL) {

  if (is.null(obabel_path)) {
    obabel_path <- Sys.getenv("GLUEDOCK_OBABEL_PATH", unset = NA)
    if (is.na(obabel_path)) {
      stop("obabel_path parameter is required and GLUEDOCK_OBABEL_PATH environment variable not set. Please provide obabel_path parameter or set GLUEDOCK_OBABEL_PATH environment variable.")
    }
  }

  if (!file.exists(obabel_path)) {
    stop("OpenBabel executable not found: ", obabel_path)
  }

  if (length(input_file) == 1 && dir.exists(input_file)) {
    files_to_convert <- list.files(
      path = input_file,
      pattern = pattern,
      full.names = TRUE,
      recursive = recursive
    )

    if (length(files_to_convert) == 0) {
      warning("No files found in directory: ", input_file)
      return(NULL)
    }

    message(sprintf("Found %d files to convert in directory: %s", length(files_to_convert), input_file))
  } else {
    files_to_convert <- input_file

    missing_files <- files_to_convert[!file.exists(files_to_convert)]
    if (length(missing_files) > 0) {
      warning("The following input files do not exist: ", paste(missing_files, collapse = ", "))
      files_to_convert <- files_to_convert[file.exists(files_to_convert)]

      if (length(files_to_convert) == 0) {
        stop("No valid input files to convert")
      }
    }
  }

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      message(sprintf("Created output directory: %s", output_dir))
    }
  }

  converted_files <- character(0)

  for (file_path in files_to_convert) {
    file_name <- tools::file_path_sans_ext(basename(file_path))
    file_ext <- tools::file_ext(file_path)

    current_input_format <- input_format
    if (is.null(current_input_format)) {
      current_input_format <- file_ext
      message(sprintf("Input format not specified for %s, using file extension: %s", basename(file_path), current_input_format))
    }

    current_output_format <- output_format
    current_output_file <- output_file

    if (is.null(current_output_format)) {
      if (is.null(current_output_file)) {
        stop("Either output_file or output_format must be specified")
      }
      current_output_format <- tools::file_ext(current_output_file)
      message(sprintf("Output format not specified for %s, using file extension: %s", basename(file_path), current_output_format))
    }

    if (is.null(current_output_file)) {
      if (is.null(output_dir)) {
        current_output_file <- paste0(file_name, ".", current_output_format)
      } else {
        current_output_file <- file.path(output_dir, paste0(file_name, ".", current_output_format))
      }
      message(sprintf("Output file not specified for %s, using: %s", basename(file_path), current_output_file))
    }

    output_file_dir <- dirname(current_output_file)
    if (!dir.exists(output_file_dir) && output_file_dir != ".") {
      dir.create(output_file_dir, recursive = TRUE)
    }

    args <- c(
      shQuote(file_path),
      "-i", current_input_format,
      "-o", current_output_format,
      "-O", shQuote(current_output_file)
    )

    message(sprintf("Converting %s to %s", basename(file_path), basename(current_output_file)))
    message("Running:\n  ", obabel_path, " ", paste(args, collapse = " "))

    res <- suppressWarnings(
      system2(
        command = obabel_path,
        args = args,
        stdout = TRUE,
        stderr = TRUE
      )
    )

    if (length(res) > 0) {
      message(paste(res, collapse = "\n"))
    }

    if (!file.exists(current_output_file)) {
      warning(sprintf("Conversion failed for %s: output file not created: %s", basename(file_path), current_output_file))
    } else {
      message(sprintf("Successfully converted %s to %s", basename(file_path), current_output_file))
      converted_files <- c(converted_files, current_output_file)
    }
  }

  total <- length(files_to_convert)
  success <- length(converted_files)
  failed <- total - success

  if (failed > 0) {
    message(sprintf("Total: %d files, Success: %d, Failed: %d", total, success, failed))
  } else {
    message(sprintf("Successfully converted all %d files", success))
  }

  if (success > 0) {
    return(converted_files)
  } else {
    return(NULL)
  }
}


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
prepare_ligand <- function(mol_files, out_dir = "./",
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

  return(invisible(unlist(valid_files)))
}
