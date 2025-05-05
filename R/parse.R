
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
#' @keywords internal
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
