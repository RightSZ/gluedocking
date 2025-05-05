#' Configure environment for GlueDocking
#'
#' This function sets up the necessary environment variables for GlueDocking by checking
#' and storing paths to required external tools in the R environment file.
#'
#' @param python_path Character, path to Python executable
#' @param prepare_receptor_script Character, path to prepare_receptor4.py script
#' @param prepare_ligand_script Character, path to prepare_ligand4.py script
#' @param prepare_split_alt_script Character, path to prepare_pdb_split_alt_confs.py script
#' @param obabel_path Character, path to OpenBabel executable
#' @param vina_path Character, path to AutoDock Vina executable
#' @param force Logical, whether to overwrite existing environment variables, default FALSE
#'
#' @return Invisibly returns TRUE if setup was successful
#' @importFrom utils file.edit
#' @importFrom tools file_path_sans_ext file_ext
#' @examples
#' \dontrun{
#' prepare_for_gluedocking(
#'   python_path = "C:/Python39/python.exe",
#'   prepare_receptor_script = "C:/MGLTools/scripts/prepare_receptor4.py",
#'   prepare_ligand_script = "C:/MGLTools/scripts/prepare_ligand4.py",
#'   prepare_split_alt_script = "C:/MGLTools/scripts/prepare_pdb_split_alt_confs.py",
#'   obabel_path = "C:/OpenBabel/bin/obabel.exe",
#'   vina_path = "C:/AutoDock_Vina/vina.exe"
#' )
#' }
#'
#' @export
prepare_for_gluedocking <- function(python_path = NULL,
                                    prepare_receptor_script = NULL,
                                    prepare_ligand_script = NULL,
                                    prepare_split_alt_script = NULL,
                                    obabel_path = NULL,
                                    vina_path = NULL,
                                    force = FALSE) {

  check_file_exists <- function(path, name) {
    if (!is.null(path)) {
      if (!file.exists(path)) {
        warning(sprintf("The specified %s file does not exist: %s", name, path))
        return(FALSE)
      }
      return(TRUE)
    }
    return(NA)
  }

  python_exists <- check_file_exists(python_path, "Python executable")
  receptor_script_exists <- check_file_exists(prepare_receptor_script, "prepare_receptor4.py script")
  ligand_script_exists <- check_file_exists(prepare_ligand_script, "prepare_ligand4.py script")
  split_alt_script_exists <- check_file_exists(prepare_split_alt_script, "prepare_pdb_split_alt_confs.py script")
  obabel_exists <- check_file_exists(obabel_path, "OpenBabel executable")
  vina_exists <- check_file_exists(vina_path, "AutoDock Vina executable")

  env_vars <- Sys.getenv()

  env_python <- Sys.getenv("GLUEDOCK_PYTHON_PATH", unset = NA)
  env_receptor <- Sys.getenv("GLUEDOCK_PREPARE_RECEPTOR", unset = NA)
  env_ligand <- Sys.getenv("GLUEDOCK_PREPARE_LIGAND", unset = NA)
  env_split_alt <- Sys.getenv("GLUEDOCK_PREPARE_SPLIT_ALT", unset = NA)
  env_obabel <- Sys.getenv("GLUEDOCK_OBABEL_PATH", unset = NA)
  env_vina <- Sys.getenv("GLUEDOCK_VINA_PATH", unset = NA)

  vars_to_write <- character()

  if (!is.null(python_path)) {
    if (!is.na(env_python)) {
      message(sprintf("GLUEDOCK_PYTHON_PATH already exists in R environment: %s", env_python))
      if (force && python_exists) {
        vars_to_write <- c(vars_to_write, sprintf('GLUEDOCK_PYTHON_PATH="%s"', python_path))
        message("Will overwrite existing Python path")
      }
    } else if (python_exists) {
      vars_to_write <- c(vars_to_write, sprintf('GLUEDOCK_PYTHON_PATH="%s"', python_path))
    }
  }

  if (!is.null(prepare_receptor_script)) {
    if (!is.na(env_receptor)) {
      message(sprintf("GLUEDOCK_PREPARE_RECEPTOR already exists in R environment: %s", env_receptor))
      if (force && receptor_script_exists) {
        vars_to_write <- c(vars_to_write, sprintf('GLUEDOCK_PREPARE_RECEPTOR="%s"', prepare_receptor_script))
        message("Will overwrite existing receptor script path")
      }
    } else if (receptor_script_exists) {
      vars_to_write <- c(vars_to_write, sprintf('GLUEDOCK_PREPARE_RECEPTOR="%s"', prepare_receptor_script))
    }
  }

  if (!is.null(prepare_ligand_script)) {
    if (!is.na(env_ligand)) {
      message(sprintf("GLUEDOCK_PREPARE_LIGAND already exists in R environment: %s", env_ligand))
      if (force && ligand_script_exists) {
        vars_to_write <- c(vars_to_write, sprintf('GLUEDOCK_PREPARE_LIGAND="%s"', prepare_ligand_script))
        message("Will overwrite existing ligand script path")
      }
    } else if (ligand_script_exists) {
      vars_to_write <- c(vars_to_write, sprintf('GLUEDOCK_PREPARE_LIGAND="%s"', prepare_ligand_script))
    }
  }
  
  if (!is.null(prepare_split_alt_script)) {
    if (!is.na(env_split_alt)) {
      message(sprintf("GLUEDOCK_PREPARE_SPLIT_ALT already exists in R environment: %s", env_split_alt))
      if (force && split_alt_script_exists) {
        vars_to_write <- c(vars_to_write, sprintf('GLUEDOCK_PREPARE_SPLIT_ALT="%s"', prepare_split_alt_script))
        message("Will overwrite existing split alt script path")
      }
    } else if (split_alt_script_exists) {
      vars_to_write <- c(vars_to_write, sprintf('GLUEDOCK_PREPARE_SPLIT_ALT="%s"', prepare_split_alt_script))
    }
  }

  if (!is.null(obabel_path)) {
    if (!is.na(env_obabel)) {
      message(sprintf("GLUEDOCK_OBABEL_PATH already exists in R environment: %s", env_obabel))
      if (force && obabel_exists) {
        vars_to_write <- c(vars_to_write, sprintf('GLUEDOCK_OBABEL_PATH="%s"', obabel_path))
        message("Will overwrite existing OpenBabel path")
      }
    } else if (obabel_exists) {
      vars_to_write <- c(vars_to_write, sprintf('GLUEDOCK_OBABEL_PATH="%s"', obabel_path))
    }
  }

  if (!is.null(vina_path)) {
    if (!is.na(env_vina)) {
      message(sprintf("GLUEDOCK_VINA_PATH already exists in R environment: %s", env_vina))
      if (force && vina_exists) {
        vars_to_write <- c(vars_to_write, sprintf('GLUEDOCK_VINA_PATH="%s"', vina_path))
        message("Will overwrite existing Vina path")
      }
    } else if (vina_exists) {
      vars_to_write <- c(vars_to_write, sprintf('GLUEDOCK_VINA_PATH="%s"', vina_path))
    }
  }

  if (length(vars_to_write) > 0) {
    message("The following environment variables will be written to R environment file:")
    message(paste(vars_to_write, collapse = "\n"))

    renviron_path <- file.path(Sys.getenv("HOME"), ".Renviron")

    if (file.exists(renviron_path)) {
      existing_content <- readLines(renviron_path, warn = FALSE)
    } else {
      existing_content <- character()
    }

    for (var_entry in vars_to_write) {
      var_name <- sub("=.*", "", var_entry)

      existing_idx <- grep(paste0("^", var_name, "="), existing_content)

      if (length(existing_idx) > 0) {

        existing_content[existing_idx[1]] <- var_entry
        message(sprintf("Updated existing entry for %s", var_name))
      } else {

        existing_content <- c(existing_content, var_entry)
        message(sprintf("Added new entry for %s", var_name))
      }
    }

    writeLines(existing_content, renviron_path)
    message(sprintf("Environment variables written to %s", renviron_path))
    message("Please restart R session to apply the environment variables")
  } else {
    message("No new environment variables to add")
  }
  invisible(TRUE)
}