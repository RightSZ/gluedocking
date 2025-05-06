
#' Prepare receptor PDB files for molecular docking
#'
#' This function processes PDB files to prepare them for molecular docking by
#' converting them to PDBQT format using AutoDock Tools scripts.
#'
#' @param pdb_files Character vector of paths to PDB files
#' @param output_dir Character, output directory for prepared files, default "./"
#' @param python_path Character, path to Python executable, default "python"
#' @param prepare_script Character, path to prepare_receptor4.py script
#' @param add_opts Character, additional options to pass to prepare_receptor4.py, default "-A hydrogens -e True"
#'
#' @return Character vector of paths to successfully prepared PDBQT files
#' @importFrom tools file_path_sans_ext
#' @examples
#' \dontrun{
#' pdbqt_files <- prepare_receptor(pdb_files,
#'                                 output_dir = "prepared",
#'                                 prepare_script = "path/to/prepare_receptor4.py")
#' }
#'
#' @export
prepare_receptor <- function(pdb_files, output_dir = "./",
                             python_path = NULL,
                             prepare_script = NULL,
                             add_opts = "-A hydrogens") {

  if (length(pdb_files) == 0) {
    stop("No PDB files provided")
  }

  pdb_files <- as.character(pdb_files)

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }

  if (is.null(prepare_script)) {
    prepare_script <- Sys.getenv("GLUEDOCK_PREPARE_RECEPTOR", unset = NA)
    if (is.na(prepare_script)) {
      stop("prepare_script parameter is required and GLUEDOCK_PREPARE_RECEPTOR environment variable not set. Please provide prepare_script parameter or set GLUEDOCK_PREPARE_RECEPTOR environment variable.")
    }
  }

  if (!file.exists(prepare_script)) {
    stop("Specified prepare_script does not exist: ", prepare_script)
  }

  process_one <- function(pdb_file) {

    if (!file.exists(pdb_file)) {
      warning(sprintf("Input PDB file does not exist: %s", pdb_file))
      return(NULL)
    }

    base_name <- tools::file_path_sans_ext(basename(pdb_file))

    if (is.null(output_dir)) {
      output_pdbqt <- paste0(base_name, ".pdbqt")
    } else {
      output_pdbqt <- file.path(output_dir, paste0(base_name, ".pdbqt"))
    }

    tool_path <- normalizePath(prepare_script, winslash = "\\", mustWork = FALSE)

    if (is.null(python_path)) {
      python_path <- Sys.getenv("GLUEDOCK_PYTHON_PATH", unset = NA)
      if (is.na(python_path)) {
        stop("Python path not provided and GLUEDOCK_PYTHON_PATH environment variable not set. Please provide python_path parameter or set GLUEDOCK_PYTHON_PATH environment variable.")
      }
    }

    py_path <- normalizePath(python_path, winslash = "\\", mustWork = FALSE)

    if (!file.exists(py_path)) {
      warning(sprintf("Python executable not found: %s", python_path))
      return(NULL)
    }

    args <- c(
      shQuote(tool_path),
      "-r", shQuote(pdb_file),
      "-o", shQuote(output_pdbqt),
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
    if (!file.exists(output_pdbqt)) {
      warning(sprintf("PDBQT not generated, please check output messages for file: %s", pdb_file))
      return(NULL)
    }

    message(sprintf("Successfully prepared receptor: %s", output_pdbqt))
    return(output_pdbqt)
  }

  out_files <- lapply(pdb_files, process_one)
  valid_files <- out_files[!sapply(out_files, is.null)]

  total <- length(pdb_files)
  success <- length(valid_files)
  failed <- total - success

  if (failed > 0) {
    message(sprintf("Total: %d files, Success: %d, Failed: %d", total, success, failed))

    failed_files <- pdb_files[sapply(out_files, is.null)]
    message("Failed files: ", paste(basename(failed_files), collapse = ", "))
  }

  return(invisible(unlist(valid_files)))
}


#' Calculate docking box parameters from PDB structures
#'
#' This function analyzes PDB files to determine appropriate docking box
#' dimensions and center coordinates based on the protein structure.
#'
#' @param pdb_files Character vector of paths to PDB files
#' @param padding Numeric, additional padding to add around the protein in Angstroms, default 5
#'
#' @return Data frame containing file names, center coordinates (x,y,z) and box dimensions (x,y,z)
#' @import bio3d
#' @examples
#' \dontrun{
#' box_params <- calculate_box(pdb_files, padding = 5)
#' }
#'
#' @export
calculate_box <- function(pdb_files, padding = 5) {

  pdb_files <- as.character(pdb_files)

  if (length(pdb_files) == 0) {
    stop("No PDB files provided")
  }

  process_one <- function(pdb_file) {

    if (!file.exists(pdb_file)) {
      warning(sprintf("PDB file does not exist: %s", pdb_file))
      return(NULL)
    }

    tryCatch({
      pdb <- bio3d::read.pdb(pdb_file)

      at <- pdb$atom[pdb$atom$type == "ATOM", c("x","y","z")]

      center <- colMeans(at)

      ranges <- apply(at, 2, range)

      size <- (ranges[2, ] - ranges[1, ]) + 2 * padding

      result <- data.frame(
        file = basename(pdb_file),
        filepath = normalizePath(pdb_file, winslash = "/", mustWork = FALSE),
        center_x = center[1],
        center_y = center[2],
        center_z = center[3],
        size_x = size[1],
        size_y = size[2],
        size_z = size[3]
      )

      return(result)
    }, error = function(e) {
      warning(sprintf("Error processing file %s: %s", pdb_file, e$message))
      return(NULL)
    })
  }

  results <- lapply(pdb_files, process_one)

  valid_results <- results[!sapply(results, is.null)]

  if (length(valid_results) == 0) {
    warning("No valid PDB files processed")
    return(data.frame(
      file = character(),
      filepath = character(),
      center_x = numeric(),
      center_y = numeric(),
      center_z = numeric(),
      size_x = numeric(),
      size_y = numeric(),
      size_z = numeric()
    ))
  }

  result_df <- do.call(rbind, valid_results)

  total <- length(pdb_files)
  success <- length(valid_results)
  failed <- total - success

  if (failed > 0) {
    message(sprintf("Total: %d files, Success: %d, Failed: %d", total, success, failed))

    failed_files <- pdb_files[sapply(results, is.null)]
    message("Failed files: ", paste(basename(failed_files), collapse=", "))
  }
  return(result_df)
}

#' Generate qvina or qvina-w config files for docking
#'
#' @param receptor_paths Character vector, full paths to receptor PDBQT files
#' @param ligand_paths   Character vector, full paths to ligand PDBQT files
#' @param box_df         Data.frame generated by calculate_box(), must contain columns:
#'                       file (basename without path), center_x/center_y/center_z,
#'                       size_x/size_y/size_z
#' @param output_dir        Output directory for configuration files
#' @param exhaustiveness Integer, search intensity, default 8
#' @return               All written configuration file paths (invisible return)
#' @importFrom tools file_path_sans_ext
#' @examples
#' \dontrun{
#' # Create configuration files for docking
#' config_files <- write_configs(
#'   receptor_paths = receptor_pdbqt,
#'   ligand_paths = ligand_pdbqt,
#'   box_df = box_params,
#'   output_dir = "configs"
#' )
#' }
#' @export
write_configs <- function(receptor_paths,
                          ligand_paths,
                          box_df,
                          output_dir = "configs",
                          exhaustiveness = 8) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  sans_ext <- function(x) {
    tolower(tools::file_path_sans_ext(basename(x)))
  }
  required_cols <- c("file","center_x", "center_y", "center_z", "size_x", "size_y", "size_z")
  missing_cols <- setdiff(required_cols, colnames(box_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in box_df: ", paste(missing_cols, collapse = ", "))
  }

  pairs <- expand.grid(receptor = receptor_paths, ligand = ligand_paths, stringsAsFactors = FALSE)

  conf_paths <- vector("character", nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    rec <- pairs$receptor[i]
    lig <- pairs$ligand[i]

    if (!file.exists(rec)) {
      warning(sprintf("Receptor file does not exist: %s", rec))
      next
    }
    if (!file.exists(lig)) {
      warning(sprintf("Ligand file does not exist: %s", lig))
      next
    }

    rec <- normalizePath(rec, winslash = "/", mustWork = FALSE)
    lig <- normalizePath(lig, winslash = "/", mustWork = FALSE)

    rec_id <- sans_ext(rec)
    lig_id <- sans_ext(lig)

    row_idx <- which(tolower(box_df$file) == rec_id)
    if (length(row_idx) == 0) {
      if ("filepath" %in% colnames(box_df)) {
        row_idx <- which(tolower(sans_ext(box_df$filepath)) == rec_id)
      }

      if (length(row_idx) == 0) {
        warning(sprintf("Cannot find center/size information for receptor '%s' in box_df, check the 'file' column.", rec_id))
        next
      }
    }

    row <- box_df[row_idx[1], , drop = FALSE]

    cfg <- c(
      sprintf("receptor = %s", rec),
      sprintf("ligand = %s", lig),
      sprintf("center_x = %f", row$center_x),
      sprintf("center_y = %f", row$center_y),
      sprintf("center_z = %f", row$center_z),
      sprintf("size_x = %f", row$size_x),
      sprintf("size_y = %f", row$size_y),
      sprintf("size_z = %f", row$size_z),
      sprintf("exhaustiveness = %d", exhaustiveness)
    )

    conf_file <- file.path(
      output_dir,
      sprintf("%s_%s_config.txt", rec_id, lig_id)
    )
    writeLines(cfg, conf_file)
    conf_paths[i] <- normalizePath(conf_file, mustWork = FALSE)
    message(sprintf("[#%d] Config written: %s", i, conf_paths[i]))
  }

  valid_paths <- conf_paths[nzchar(conf_paths)]

  total <- nrow(pairs)
  success <- length(valid_paths)
  failed <- total - success

  if (failed > 0) {
    message(sprintf("Total: %d pairs, Success: %d, Failed: %d", total, success, failed))
  } else {
    message(sprintf("Successfully generated %d configuration files", success))
  }

  invisible(valid_paths)
}


