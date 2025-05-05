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
