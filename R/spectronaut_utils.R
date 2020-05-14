# MaxQuant data import utils
#
# Author: astukalov
###############################################################################

Spectronaut_NAs <- c("", "Filtered")


#' @export
read.Spectronaut.ProteinsReport <- function(file, nrows = Inf, import_data = c(), guess_max = min(10000L, nrows), delim=delim)
{
    proteinGroups.df <- readr::read_delim(file, delim = delim, n_max = nrows,
                                          #col_types = readr::cols(`Fasta headers` = "c", `id` = "i"),
                                          na = Spectronaut_NAs, guess_max = guess_max)
    col_renames <- c("protgroup_sn_id" = "PG.ProteinGroups",
                     #"protein_acs" = "PG.UniProtIds",
                     "majority_protein_acs" = "PG.ProteinAccessions",
                     "gene_names" = "PG.Genes", "protein_names" = "PG.ProteinNames",
                     "protein_descriptions" = "PG.ProteinDescriptions",
                     "q_value" = "PG.Qvalue")
    res.df <- dplyr::select(proteinGroups.df, !!col_renames[col_renames %in% colnames(proteinGroups.df)]) %>%
        dplyr::mutate(protgroup_id = row_number() - 1L,
                      protein_acs = majority_protein_acs)
    col_info <- list(protgroup = colnames(res.df))
    if ('quantity' %in% import_data) {
        quantities.df <- proteinGroups.df %>% dplyr::select(tidyselect::matches("\\.PG\\.(Quantity|RunEvidenceCount|NrOfStrippedSequencesMeasured|NrOfStippedSequencesIdentified|NrOfPrecursorsIdentified)"))
        res.df <- dplyr::bind_cols(res.df, quantities.df)
        col_info$quantity <- colnames(quantities.df)
    }
    attr(res.df, "column_groups") <- col_info
    return (res.df)
}

#' @export
read.Spectronaut.PepmodstatesReport <- function(file, nrows = Inf, import_data = c(), guess_max = min(10000L, nrows), delim=',')
{
    proteinGroups.df <- readr::read_delim(file, delim = delim, n_max = nrows,
                                          #col_types = readr::cols(`Fasta headers` = "c", `id` = "i"),
                                          na = Spectronaut_NAs, guess_max = guess_max)
    col_renames <- c("protgroup_sn_id" = "PG.ProteinGroups",
                     #"protein_acs" = "PG.UniProtIds",
                     "majority_protein_acs" = "PG.ProteinAccessions",
                     "gene_names" = "PG.Genes", "protein_names" = "PG.ProteinNames",
                     "protein_descriptions" = "PG.ProteinDescriptions",
                     "q_value" = "PG.Qvalue",
                     "pepmodstate_seq" = "EG.PrecursorId")
    res.df <- dplyr::select(proteinGroups.df, !!!col_renames[col_renames %in% colnames(proteinGroups.df)]) %>%
        dplyr::mutate(pepmodstate_id = row_number() - 1L,
                      protein_acs = majority_protein_acs) %>%
        tidyr::extract(pepmodstate_seq, c("pepmod_seq", "charge"), "_([^.]+)_\\.(\\d+)", remove=FALSE) %>%
        dplyr::mutate(charge = parse_integer(charge),
                      peptide_seq = str_remove_all(pepmod_seq, "\\[[^]]+\\]"))
    col_info <- list(pepmodstate = colnames(res.df))
    if ('pg_stats' %in% import_data) {
        pg_stats.df <- proteinGroups.df %>% dplyr::select(matches("\\.PG\\.(RunEvidenceCount|NrOfModifiedSequencesMeasured|NrOfModifiedSequencesIdentified|NrOfStrippedSequencesIdentified)"))
        res.df <- dplyr::bind_cols(res.df, pg_stats.df)
        col_info$quantity <- colnames(pg_stats.df)
    }
    if ('quantity' %in% import_data) {
        intensities.df <- proteinGroups.df %>% dplyr::select(matches("\\.PEP\\.(RunEvidenceCount|Quantity)"))
        res.df <- dplyr::bind_cols(res.df, intensities.df)
        col_info$quantity <- colnames(intensities.df)
    }
    if ('eg_quantity' %in% import_data) { # FIXME what's that?
        intensities.df <- proteinGroups.df %>% dplyr::select(matches("\\EG\\.(TotalQuantity \\(Settings\\))"))
        res.df <- dplyr::bind_cols(res.df, intensities.df)
        col_info$eg_quantity <- colnames(intensities.df)
    }
    attr(res.df, "column_groups") <- col_info
    return (res.df)
}

SpectronautMetrics <- c("RunEvidenceCount" = "nevidences",
                        "IsSingleHit" = "is_singlehit",
                        "Quantity" = "intensity",
                        "NrOfStrippedSequencesMeasured" = "npeptides_quanted",
                        "NrOfStrippedSequencesIdentified" = "npeptides_idented",
                        "NrOfModifiedSequencesMeasured" = "npepmods_quanted",
                        "NrOfModifiedSequencesIdentified" = "npepmods_idented",
                        "NrOfPrecursorsIdentified" = "npepmodstates_idented",
                        "TotalQuantity (Settings)" = "total_intensity"
                        )


pivot_longer.Spectronaut.Metrics <- function(msdata.wide, pkey, colgroup) {
    colgroups <- attr(msdata.wide, "column_groups")
    cols <- colgroups[[colgroup]]

    quantity_prespec_df <- tibble(.name = cols) %>%
      tidyr::extract(.name, c("msrun_ix", "raw_file", ".value"), remove = FALSE,
                     "^\\[(\\d+)\\]\\s(.+)(?:\\.PG|\\.PEP)\\.(.+)$") %>%
      dplyr::mutate(.value = SpectronautMetrics[.value])

    return (tidyr::pivot_longer_spec(
        dplyr::select(msdata.wide, !!pkey, !!cols),
        quantity_prespec_df) %>%
        dplyr::mutate(msrun_ix = parse_integer(msrun_ix)))
}

#' @export
pivot_longer.Spectronaut.ProtgroupIntensities <- function(msdata.wide) {
    return (pivot_longer.Spectronaut.Metrics(msdata.wide, "protgroup_id", "quantity"))
}

#' @export
pivot_longer.Spectronaut.PepmodstateIntensities <- function(msdata.wide) {
    return (pivot_longer.Spectronaut.Metrics(msdata.wide, "pepmodstate_id", "quantity"))
}
