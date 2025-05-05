#' Trim receptor PDB files to keep only protein atoms
#'
#' @param input_paths Character vector of file paths or directories containing PDB files
#' @param output_dir Character string specifying the output directory, default is "trimmed"
#' @return Character vector of paths to the trimmed PDB files
#' @export
#' @examples
#' \dontrun{
#' trim_receptor("./test/1a7c.pdb")
#' trim_receptor(c("./test/1a7c.pdb", "./test/1abc.pdb"), output_dir = "my_output")
#' trim_receptor("./test", output_dir = "processed")
#' }
trim_receptor <- function(input_paths, output_dir = "trimmed") {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }


  process_pdb <- function(pdb_file) {
    tryCatch({
      # Read PDB file
      message("Processing: ", pdb_file)
      pdb <- bio3d::read.pdb(pdb_file)

      # Select only protein atoms
      sel <- bio3d::atom.select(pdb, "protein")

      # Trim PDB to keep only protein atoms
      clean <- bio3d::trim.pdb(pdb, sel)

      # Generate output filename
      file_basename <- basename(pdb_file)
      output_filename <- file.path(output_dir, paste0(tools::file_path_sans_ext(file_basename), ".trimmed.pdb"))

      # Write trimmed PDB
      bio3d::write.pdb(clean, file = output_filename)

      message("Saved trimmed PDB to: ", output_filename)
      return(output_filename)
    }, error = function(e) {
      warning("Error processing file ", pdb_file, ": ", e$message)
      return(NULL)
    })
  }

  # Process all input paths
  output_files <- c()

  for (path in input_paths) {
    if (dir.exists(path)) {
      # If path is a directory, process all PDB files in it
      pdb_files <- list.files(path, pattern = "\\.pdb$", full.names = TRUE)
      if (length(pdb_files) == 0) {
        warning("No PDB files found in directory: ", path)
        next
      }

      for (pdb_file in pdb_files) {
        result <- process_pdb(pdb_file)
        if (!is.null(result)) {
          output_files <- c(output_files, result)
        }
      }
    } else if (file.exists(path) && grepl("\\.pdb$", path)) {
      # If path is a PDB file, process it
      result <- process_pdb(path)
      if (!is.null(result)) {
        output_files <- c(output_files, result)
      }
    } else {
      warning("Path does not exist or is not a PDB file: ", path)
    }
  }

  if (length(output_files) == 0) {
    warning("No PDB files were processed successfully")
  } else {
    message("Successfully processed ", length(output_files), " PDB files")
  }

  return(invisible(output_files))
}
