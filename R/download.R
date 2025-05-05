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

  return(invisible(valid_files))
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
  return(invisible(valid_files))
}
