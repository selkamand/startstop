
#' Find start and stop codons
#'
#' Searches a fasta sequence for all start / stop codons
#'
#'
#' @param fasta_path path to a fasta file (DNA) (string)
#' @param start  a vector of 'start' codons to search for (character)
#' @param stop a vector of 'stop' codons to search for (character)
#' @param write_bedfiles write bed files describing positions of start/stop codons in each frame (boolean)
#'
#' @return a dataframe describing the positions of start/stop-codons found. Will also write bed files describing start/stop positions if \code{write_bedfiles == true}.
#' @export
#'
find_start_and_stop_codons <- function(fasta_path, start = c('TTG', 'CTG', 'ATG'), stop=c('TAA', 'TAG', 'TGA'), write_bedfiles=TRUE){
  fasta_path <- fasta_path[1]
  assertthat::assert_that(file.exists(fasta_path))
  sequence_to_search <- seqinr::read.fasta(file = fasta_path, seqtype = "DNA")

  message("Searching for Start Codons: [",  paste0(start, collapse = ","), "]", " and Stop Codons: [", paste0(stop, collapse = ",") ,"]")

  if(length(sequence_to_search) > 1){
    message("[Warning] Only identifying start and stop codons in the first FASTA entry")
  }

  seqname=names(sequence_to_search)[1]
  seqlength = length(sequence_to_search[[1]])

  #Start Codons
  start_codons_df <- purrr::map_dfr(c(-3, -2, -1, 1, 2,3), .f = function(frame){
    assertthat::assert_that(frame %in% c(-3, -2, -1, 1, 2, 3))
    if (frame < 0) {
      seq_to_search <- rev(seqinr::comp(sequence_to_search[[1]]))
    }
    else
    {
      seq_to_search <- sequence_to_search[[1]]
    }

    seqinr_frame = abs(frame)-1

    seq = seqinr::splitseq(seq_to_search, frame = seqinr_frame)
    dplyr::tibble(
      codon_number=seq %>%
        grep(paste0(start, collapse="|"), x = ., ignore.case = TRUE),
      frame=frame,
      codon=seq[codon_number],
      start=codon_number*3-3+seqinr_frame,
      class="start"
    )
  })

  stop_codons_df <- purrr::map_dfr(c(-3, -2, -1, 1, 2,3), .f = function(frame){
    assertthat::assert_that(frame %in% c(-3, -2, -1, 1, 2, 3))
    if (frame < 0) {
      seq_to_search <- rev(seqinr::comp(sequence_to_search[[1]]))
    }
    else
    {
      seq_to_search <- sequence_to_search[[1]]
    }

    seqinr_frame = abs(frame)-1

    seq = seqinr::splitseq(seq_to_search, frame = seqinr_frame)
    dplyr::tibble(
      codon_number=seq %>%
        grep(paste0(stop, collapse="|"), x = ., ignore.case = TRUE),
      frame=frame,
      codon=seq[codon_number],
      start=codon_number*3-3+seqinr_frame,
      class="stop"
    )
  })

  start_and_stop_codons_df <- rbind(start_codons_df, stop_codons_df)

  #put coords in terms of positive strand, even for -ve
  message("seqlength: ", seqlength)
  start_and_stop_codons_df <- start_and_stop_codons_df %>%
    dplyr::mutate(
      start = ifelse(frame > 0 , yes = start, no=seqlength-start-3),
      end = start+3
      )

  start_and_stop_codons_df <- start_and_stop_codons_df %>%
    dplyr::mutate(
      seqname = seqname,
      name = paste0(class, "_", codon),
      itemRGB = ifelse(class=="start", yes = "0,230,0", no = "255,0,0"),
      thick_start = start,
      thick_end = end,
      score = ".",
      strand = ifelse(frame > 0, yes = "+", no = "-")
      )


  if(write_bedfiles){
    message("Writing each frames start/stop codons to bedfiles. Try visualising in IGV")
    start_and_stop_codons_df %>%
      dplyr::group_by(frame) %>%
      dplyr::group_walk( ~ .x %>%
          dplyr::select(seqname, start, end, name, score, strand, thick_start, thick_end, itemRGB) %>%
          write.table(file = paste0(seqname, "_start_stop_", .y$frame, ".bed"), col.names = FALSE, row.names = FALSE, sep = "\t", quote=FALSE)
      )
  }

  start_and_stop_codons_df %>%
    dplyr::select(-thick_start, -thick_end, -score) %>%
    return()
}
