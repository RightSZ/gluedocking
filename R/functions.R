
#' Download PDB files from RCSB PDB database
#'
#' This function downloads PDB files from the RCSB Protein Data Bank using their PDB IDs.
#'
#' @param pdb_ids Character vector of PDB IDs to download
#' @param out_dir Character, output directory for downloaded files, default "./"
#' @param verify_ssl Logical, whether to verify SSL certificates, default TRUE
#'
#' @return List of successfully downloaded file paths
#' @import httr
#' @examples
#' \dontrun{
#' files <- download_receptor(c("1iep", "4hg7"), out_dir = "receptors")
#' }
#'
#' @export
download_receptor <- function(pdb_ids, out_dir = "./", verify_ssl = TRUE) {
  pdb_ids <- as.character(pdb_ids)

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  download_one <- function(pdb_id) {
    pdb_id <- tolower(pdb_id)
    url <- sprintf("https://files.rcsb.org/download/%s.pdb", pdb_id)
    dest <- file.path(out_dir, paste0(pdb_id, ".pdb"))

    config <- httr::config(ssl_verifypeer = verify_ssl)
    resp <- tryCatch(
      httr::GET(url, config),
      error = function(e) {
        if (verify_ssl) {
          message(sprintf("SSL error occurred, retrying without SSL verification...", pdb_id))
          return(httr::GET(url, httr::config(ssl_verifypeer = FALSE)))
        } else {
          stop(e)
        }
      }
    )
    if (httr::status_code(resp) != 200 && verify_ssl) {
      message(sprintf("Download failed with SSL verification, retrying without verification...",
                      pdb_id, httr::status_code(resp)))
      resp <- httr::GET(url, httr::config(ssl_verifypeer = FALSE))
    }

    status_code <- httr::status_code(resp)
    if (status_code != 200) {
      return(NULL)
    } else {
      writeBin(httr::content(resp, "raw"), dest)
      message(sprintf("[%s] Download: %s", pdb_id, dest))
      return(dest)
    }
  }

  out_files <- lapply(pdb_ids, download_one)
  valid_files <- out_files[!sapply(out_files, is.null)]

  total <- length(pdb_ids)
  success <- length(valid_files)
  failed <- total - success

  if (failed > 0) {
    message(sprintf("Total: %d IDs, Success: %d, Failed: %d", total, success, failed))

    failed_ids <- pdb_ids[sapply(out_files, is.null)]
    message("Failed IDs: ", paste(failed_ids, collapse=", "))
  }

  return(valid_files)
}


#' Prepare receptor PDB files for molecular docking
#'
#' This function processes PDB files to prepare them for molecular docking by
#' converting them to PDBQT format using AutoDock Tools scripts.
#'
#' @param pdb_files Character vector of paths to PDB files
#' @param out_dir Character, output directory for prepared files, default "./"
#' @param python_path Character, path to Python executable, default "python"
#' @param prepare_script Character, path to prepare_receptor4.py script
#' @param add_opts Character, additional options to pass to prepare_receptor4.py, default "-A hydrogens -e True"
#'
#' @return Character vector of paths to successfully prepared PDBQT files
#' @importFrom tools file_path_sans_ext
#' @examples
#' \dontrun{
#' pdbqt_files <- prepare_receptor(pdb_files,
#'                                 out_dir = "prepared",
#'                                 prepare_script = "path/to/prepare_receptor4.py")
#' }
#'
#' @export
prepare_receptor <- function(pdb_files, out_dir = "./",
                             python_path = NULL,
                             prepare_script = NULL,
                             add_opts = "-A hydrogens") {

  if (length(pdb_files) == 0) {
    stop("No PDB files provided")
  }

  pdb_files <- as.character(pdb_files)

  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
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

    if (is.null(out_dir)) {
      output_pdbqt <- paste0(base_name, ".pdbqt")
    } else {
      output_pdbqt <- file.path(out_dir, paste0(base_name, ".pdbqt"))
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

  return(unlist(valid_files))
}

#' Download ligand structures from PubChem
#'
#' This function downloads chemical compound structures from PubChem database using their CIDs.
#'
#' @param cids Character vector of PubChem Compound IDs (CIDs) to download
#' @param out_dir Character, output directory for downloaded files, default "ligands"
#' @param verify_ssl Logical, whether to verify SSL certificates, default TRUE
#'
#' @return List of successfully downloaded file paths
#' @import httr
#' @examples
#' \dontrun{
#' files <- download_ligand(c("2244", "5090"), out_dir = "compounds")
#' }
#'
#' @export
download_ligand <- function(cids, out_dir = "ligands", verify_ssl = TRUE) {
  cids <- as.character(cids)

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  download_one <- function(cid) {
    url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/%s/SDF", cid)
    dest <- file.path(out_dir, paste0(cid, ".sdf"))

    config <- httr::config(ssl_verifypeer = verify_ssl)
    resp <- tryCatch(
      httr::GET(url, config),
      error = function(e) {
        if (verify_ssl) {
          message(sprintf("SSL error occurred, retrying without SSL verification...", cid))
          return(httr::GET(url, httr::config(ssl_verifypeer = FALSE)))
        } else {
          stop(e)
        }
      }
    )

    if (httr::status_code(resp) != 200 && verify_ssl) {
      message(sprintf("Download failed with SSL verification, retrying without verification...",
                      cid, httr::status_code(resp)))
      resp <- httr::GET(url, httr::config(ssl_verifypeer = FALSE))
    }

    status_code <- httr::status_code(resp)
    if (status_code != 200) {
      if (status_code == 404) {
        warning(sprintf("CID %s does not exist or is invalid (status: 404)", cid))
      } else {
        warning(sprintf("Failed to download ligand CID %s (status %s)",
                        cid, status_code))
      }
      return(NULL)
    } else {
      writeBin(httr::content(resp, "raw"), dest)
      message(sprintf("[%s] Download: %s", cid, dest))
      return(dest)
    }
  }

  out_files <- lapply(cids, download_one)
  valid_files <- out_files[!sapply(out_files, is.null)]

  total <- length(cids)
  success <- length(valid_files)
  failed <- total - success

  if (failed > 0) {
    message(sprintf("Total: %d CIDs, Success: %d, Failed: %d", total, success, failed))

    failed_cids <- cids[sapply(out_files, is.null)]
    message("Failed CIDs: ", paste(failed_cids, collapse=", "))
  }
  return(valid_files)
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

  return(unlist(valid_files))
}


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
#'
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
#' @param out_dir        Output directory for configuration files
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
#'   out_dir = "configs"
#' )
#' }
#' @export
write_configs <- function(receptor_paths,
                          ligand_paths,
                          box_df,
                          out_dir = "configs",
                          exhaustiveness = 8) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
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
      out_dir,
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


#' Run molecular docking with AutoDock Vina
#'
#' This function runs molecular docking using AutoDock Vina.
#'
#' @param receptor Path to receptor PDBQT file (only for direct mode)
#' @param ligand Path to ligand PDBQT file (only for direct mode)
#' @param center Numeric vector of length 3 for box center (x,y,z) (only for direct mode)
#' @param size Numeric vector of length 3 for box size (x,y,z) (only for direct mode)
#' @param config_paths Character vector of paths to config files or directories containing config files
#' @param exhaustiveness Integer, search intensity, default 8 (only for direct mode)
#' @param out Character, output directory for docked files, default "docked"
#' @param logs Character, output directory for log files, default "logs"
#' @param vina_path Character, path to vina executable, default NULL
#' @param seed Integer, random seed for reproducibility, default 12345
#' @param cpu Integer, number of CPU cores to use for parallel computation, default NULL (auto-detect)
#' @importFrom tools file_path_sans_ext
#' @examples
#' \dontrun{
#' # Run molecular docking using AutoDock Vina
#' run_vina(
#'   config_paths = "configs",
#'   out = "docked",
#'   logs = "logs"
#' )
#' }
#' @export
run_vina <- function(receptor = NULL, ligand = NULL, center = NULL, size = NULL,
                     config_paths = NULL, exhaustiveness = 8, out = "docked",
                     logs = "logs", vina_path = NULL, seed = 12345, cpu = NULL) {

  if (!is.null(cpu)) {
    if (!is.numeric(cpu) || cpu != as.integer(cpu) || cpu <= 0) {
      warning(sprintf("Provided CPU value '%s' is not a valid positive integer, using auto-detection", cpu))
      cpu <- NULL
    }
  }

  if (is.null(vina_path)) {
    vina_path <- Sys.getenv("GLUEDOCK_VINA_PATH", unset = NA)
    if (is.na(vina_path)) {
      stop("vina_path parameter is required and GLUEDOCK_VINA_PATH environment variable not set. Please provide vina_path parameter or set GLUEDOCK_VINA_PATH environment variable.")
    }
  }

  if (!file.exists(vina_path)) {
    stop("AutoDock Vina executable does not exist: ", vina_path)
  }

  if (is.null(config_paths)) {
    if (is.null(receptor) || is.null(ligand) || is.null(center) || is.null(size)) {
      stop("In direct mode, receptor, ligand, center, and size must be provided")
    }

    if (!file.exists(receptor)) {
      stop("Receptor file does not exist: ", receptor)
    }

    if (!file.exists(ligand)) {
      stop("Ligand file does not exist: ", ligand)
    }

    if (length(center) != 3) {
      stop("Center must be a numeric vector of length 3")
    }

    if (length(size) != 3) {
      stop("Size must be a numeric vector of length 3")
    }

    receptor <- normalizePath(receptor, winslash = "/", mustWork = FALSE)
    ligand <- normalizePath(ligand, winslash = "/", mustWork = FALSE)

    if (!is.null(seed)) {
      if (!is.numeric(seed) || seed != as.integer(seed) || seed <= 0) {
        warning(sprintf("Provided seed value '%s' is not a valid positive integer, using default value %d", seed, 12345))
        seed <- 12345
      }
    }

    if (!dir.exists(out)) {
      dir.create(out, recursive = TRUE)
    }

    if (!dir.exists(logs)) {
      dir.create(logs, recursive = TRUE)
    }

    rec_id <- tolower(tools::file_path_sans_ext(basename(receptor)))
    lig_id <- tolower(tools::file_path_sans_ext(basename(ligand)))

    out_file <- file.path(out, paste0(rec_id, "_", lig_id, ".pdbqt"))
    log_file <- file.path(logs, paste0(rec_id, "_", lig_id, "_log.txt"))

    cmd <- sprintf(
      "%s --receptor %s --ligand %s --center_x %f --center_y %f --center_z %f --size_x %f --size_y %f --size_z %f --exhaustiveness %d --out %s --seed %d --log %s",
      shQuote(vina_path),
      shQuote(receptor),
      shQuote(ligand),
      center[1], center[2], center[3],
      size[1], size[2], size[3],
      exhaustiveness,
      shQuote(out_file),
      seed,
      shQuote(log_file)
    )

    if (!is.null(cpu)) {
      cmd <- paste0(cmd, sprintf(" --cpu %d", as.integer(cpu)))
    }

    message("Running direct docking command:\n  ", cmd)

    tryCatch(
      system(cmd, intern = TRUE),
      error = function(e) {
        warning("Error executing vina: ", e$message)
      }
    )

  } else {

    config_files <- character(0)

    for (path in config_paths) {
      if (dir.exists(path)) {

        files <- list.files(
          path = path,
          pattern = ".*_config\\.txt$",
          full.names = TRUE
        )
        config_files <- c(config_files, files)
      } else if (file.exists(path)) {

        config_files <- c(config_files, path)
      } else {
        warning("Config path does not exist: ", path)
      }
    }

    if (length(config_files) == 0) {
      stop("No valid config files found")
    }

    if (!is.null(seed)) {
      if (!is.numeric(seed) || seed != as.integer(seed) || seed <= 0) {
        warning(sprintf("Provided seed value '%s' is not a valid positive integer, using default value %d", seed, 12345))
        seed <- 12345
      }
    }

    if (!dir.exists(out)) {
      dir.create(out, recursive = TRUE)
    }

    if (!dir.exists(logs)) {
      dir.create(logs, recursive = TRUE)
    }

    results <- list()

    for (i in seq_along(config_files)) {
      config_file <- config_files[i]

      config_base <- tools::file_path_sans_ext(basename(config_file))
      config_base <- sub("_config$", "", config_base)  # Remove _config suffix if present

      out_file <- file.path(out, paste0(config_base, ".pdbqt"))
      log_file <- file.path(logs, paste0(config_base, "_log.txt"))

      cmd <- sprintf(
        "%s --config %s --out %s --seed %d --log %s",
        shQuote(vina_path),
        shQuote(config_file),
        shQuote(out_file),
        seed,
        shQuote(log_file)
      )

      if (!is.null(cpu)) {
        cmd <- paste0(cmd, sprintf(" --cpu %d", as.integer(cpu)))
      }

      message(sprintf("[#%d/%d] Running config docking command:\n  %s",
                      i, length(config_files), cmd))

      tryCatch(
        system(cmd, intern = TRUE),
        error = function(e) {
          warning("Error executing vina with config ", config_file, ": ", e$message)
        }
      )
    }

    message(sprintf("Processed %d config files", length(config_files)))
  }
}



#' Check docking log files and rerun failed docking jobs
#'
#' @param logs Character, directory containing log files, default "logs"
#' @param out Character, output directory for docked files, default "docked"
#' @param config_paths Character vector of paths to config files or directories containing config files
#' @param vina_path Character, path to vina executable, default "vina"
#' @param seed Integer, random seed for reproducibility, default 12345
#' @param cpu Integer, number of CPU cores to use for parallel computation, default NULL (auto-detect)
#' @param force Logical, whether to force rerun all jobs regardless of log status, default FALSE
#' @return Character vector of rerun log files
#' @importFrom tools file_path_sans_ext
#' @examples
#' \dontrun{
#' # Check log files and rerun any failed docking jobs
#' check_logs(
#'   logs = "logs",
#'   out = "docked",
#'   config_paths = "configs"
#' )
#' }
#' @export
check_logs <- function(logs = "logs", out = "docked", config_paths = NULL,
                       vina_path = NULL, seed = 12345, cpu = NULL, force = FALSE) {
  if (!dir.exists(logs)) {
    warning(sprintf("Log directory does not exist: %s", logs))
    return(invisible(NULL))
  }

  log_files <- list.files(path = logs, pattern = "\\.txt$", full.names = TRUE)

  if (length(log_files) == 0) {
    message("No log files found in directory: ", logs)
    return(invisible(NULL))
  }

  incomplete_logs <- c()

  for (file in log_files) {
    if (force) {
      incomplete_logs <- c(incomplete_logs, file)
      next
    }

    lines <- tryCatch(
      readLines(file, warn = FALSE),
      error = function(e) {
        warning(sprintf("Error reading log file %s: %s", file, e$message))
        return(character(0))
      }
    )

    if (length(lines) == 0) {
      incomplete_logs <- c(incomplete_logs, file)
      next
    }

    header_found <- any(grepl("mode \\|.*affinity", lines))
    separator_found <- any(grepl("-----\\+------------\\+----------\\+----------", lines))

    if (!(header_found && separator_found)) {
      incomplete_logs <- c(incomplete_logs, file)
    }
  }

  total <- length(log_files)
  incomplete <- length(incomplete_logs)

  if (incomplete == 0) {
    message("All ", total, " log files contain docking results.")
    return(invisible(character(0)))
  }

  message(sprintf("Found %d/%d incomplete log files that need to be rerun.", incomplete, total))

  rerun_configs <- c()

  for (log_file in incomplete_logs) {
    base_name <- tools::file_path_sans_ext(basename(log_file))
    base_name <- sub("_log$", "", base_name)

    if (!is.null(config_paths)) {
      for (config_path in config_paths) {
        if (dir.exists(config_path)) {
          config_file <- file.path(config_path, paste0(base_name, "_config.txt"))
          if (file.exists(config_file)) {
            rerun_configs <- c(rerun_configs, config_file)
            break
          }
        } else if (file.exists(config_path) && grepl(base_name, config_path)) {
          rerun_configs <- c(rerun_configs, config_path)
          break
        }
      }
    }
  }

  if (length(rerun_configs) == 0 && incomplete > 0) {
    warning("Could not find config files for incomplete logs. Please provide valid config_paths.")
    return(invisible(incomplete_logs))
  }

  if (length(rerun_configs) > 0) {
    message(sprintf("Rerunning %d docking jobs...", length(rerun_configs)))
    run_vina(
      config_paths = rerun_configs,
      out = out,
      logs = logs,
      vina_path = vina_path,
      seed = seed,
      cpu = cpu
    )

    message("Completed rerunning docking jobs.")
  }

  return(invisible(incomplete_logs))
}

#' Parse QVina docking log files
#'
#' This function extracts docking results from QVina log files, including
#' binding affinities and RMSD values for each docking pose.
#'
#' @param file_path Character, path to the QVina log file
#'
#' @return Data frame containing docking results with columns for mode, affinity (kcal/mol),
#'         and RMSD values, or NULL if parsing fails
#' @importFrom utils read.table
#' @examples
#' \dontrun{
#' results <- parse_qvina_log("docking_results.txt")
#' }
#'
#' @export
parse_qvina_log <- function(file_path) {
  lines <- readLines(file_path, warn = FALSE)

  header_idx <- grep("mode \\|.*affinity", lines)
  separator_idx <- grep("-----\\+------------\\+----------\\+----------", lines)

  if(length(header_idx) == 0 || length(separator_idx) == 0) {
    return(NULL)
  }

  data_start <- separator_idx[1] + 1
  data_lines <- lines[data_start:length(lines)]

  result_lines <- data_lines[grepl("^\\s*\\d+", data_lines)]

  if(length(result_lines) == 0) {
    return(NULL)
  }

  df <- read.table(text = result_lines, header = FALSE, fill = TRUE)

  names(df)[1:4] <- c("mode", "affinity_kcalmol", "rmsd_lb", "rmsd_ub")

  return(df)
}

#' Parse all log files in a directory and combine results
#'
#' @param log_folder Character, path to the directory containing log files, default "logs"
#' @param pattern Character, pattern to match log files, default "\\.txt$"
#' @param verbose Logical, whether to print detailed messages, default TRUE
#' @return A data frame containing combined docking results from all log files
#' @examples
#' \dontrun{
#' # Parse all log files and combine results
#' results <- parse_logs(log_folder = "logs")
#' }
#' @export
parse_logs <- function(log_folder = "logs", pattern = "\\.txt$", verbose = TRUE) {
  if (!dir.exists(log_folder)) {
    stop(sprintf("Log directory does not exist: %s", log_folder))
  }

  log_files <- list.files(path = log_folder, pattern = pattern, full.names = TRUE)

  if (length(log_files) == 0) {
    warning(sprintf("No log files found in directory: %s", log_folder))
    return(NULL)
  }

  if (verbose) {
    message(sprintf("Found %d log files to process", length(log_files)))
  }

  all_results <- list()
  failed_files <- character(0)

  for (f in log_files) {
    df <- parse_qvina_log(f)
    if (!is.null(df)) {
      fname <- basename(f)
      name_without_suffix <- sub("_log\\.txt$", "", fname)
      parts <- strsplit(name_without_suffix, "_")[[1]]
      if (length(parts) >= 2) {
        protein <- parts[1]
        ligand <- paste(parts[-1], collapse = "_")
      } else {
        protein <- name_without_suffix
        ligand <- NA
      }
      df$protein <- protein
      df$ligand <- ligand
      df$file <- fname
      df$path <- f

      all_results[[f]] <- df
    } else {
      failed_files <- c(failed_files, f)
      if (verbose) {
        message(sprintf("No docking results found in file: %s", basename(f)))
      }
    }
  }

  total <- length(log_files)
  success <- length(all_results)
  failed <- length(failed_files)

  if (verbose) {
    message(sprintf("Total: %d files, Success: %d, Failed: %d", total, success, failed))

    if (failed > 0) {
      message("Failed files: ", paste(basename(failed_files), collapse = ", "))
    }
  }

  if (length(all_results) == 0) {
    warning("No valid docking results found in any log file")
    return(NULL)
  }
  combined_df <- do.call(rbind, all_results)
  rownames(combined_df) <- NULL
  return(combined_df)
}
