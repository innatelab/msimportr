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
        quantities.df <- proteinGroups.df %>% dplyr::select(tidyselect::matches("\\.PG\\.Quantity"))
        res.df <- dplyr::bind_cols(res.df, quantities.df)
        col_info$quantity <- colnames(quantities.df)
    }
    if ('runstats' %in% import_data) {
      runstats.df <- proteinGroups.df %>% dplyr::select(tidyselect::matches("\\.PG\\.(RunEvidenceCount|NrOfStrippedSequencesMeasured|NrOfStippedSequencesIdentified|NrOfPrecursorsIdentified)"))
      res.df <- dplyr::bind_cols(res.df, runstats.df)
      col_info$runstats <- colnames(runstats.df)
    }
    attr(res.df, "column_groups") <- col_info
    return (res.df)
}

#' @export
read.Spectronaut.PepmodstatesReport <- function(file, nrows = Inf, import_data = c(), guess_max = min(10000L, nrows), delim=',')
{
    report.df <- readr::read_delim(file, delim = delim, n_max = nrows,
                                   #col_types = readr::cols(`Fasta headers` = "c", `id` = "i"),
                                   na = Spectronaut_NAs, guess_max = guess_max)
    col_renames <- c("protein_acs" = "PEP.AllOccurringProteinAccessions",
                     "is_pg_specific" = "PEP.IsProteinGroupSpecific",
                     "peptide_poses" = "PEP.PeptidePosition",
                     "unmod_seq" = "PEP.StrippedSequence",
                     "mod_seq" = "EG.ModifiedSequence",
                     "pepmodstate_seq" = "EG.PrecursorId",
                     "is_decoy" = "EG.IsDecoy")
    res.df <- dplyr::select(report.df, !!!col_renames[col_renames %in% colnames(report.df)]) %>%
        dplyr::mutate(pepmodstate_id = row_number() - 1L) %>%
        tidyr::extract(pepmodstate_seq, c("pepmod_seq", "charge"), "_([^.]+)_\\.(\\d+)", remove=FALSE) %>%
        dplyr::mutate(charge = parse_integer(charge),
                      peptide_seq = str_remove_all(pepmod_seq, "\\[[^]]+\\]"))
    if (!rlang::has_name(res.df, "protein_acs")) {
      res.df <- dplyr::mutate(res.df, protein_acs = majority_protein_acs)
    }
    col_info <- list(pepmodstate = colnames(res.df))
    if ('pg_quantity' %in% import_data) {
      intensities.df <- report.df %>% dplyr::select(matches("\\.PG\\.Quantity"))
      res.df <- dplyr::bind_cols(res.df, intensities.df)
      col_info$pg_quantity <- colnames(intensities.df)
    }
    if ('pg_runstats' %in% import_data) {
        pg_stats.df <- report.df %>% dplyr::select(matches("\\.PG\\.(RunEvidenceCount|NrOfModifiedSequencesMeasured|NrOfModifiedSequencesIdentified|NrOfStrippedSequencesIdentified)"))
        res.df <- dplyr::bind_cols(res.df, pg_stats.df)
        col_info$pg_stats <- colnames(pg_stats.df)
    }
    if ('pep_quantity' %in% import_data) {
        intensities.df <- report.df %>% dplyr::select(matches("\\.PEP\\.(RunEvidenceCount|Quantity)"))
        res.df <- dplyr::bind_cols(res.df, intensities.df)
        col_info$pep_quantity <- colnames(intensities.df)
    }
    if ('eg_quantity' %in% import_data) { # FIXME what's that?
      intensities.df <- report.df %>% dplyr::select(matches("\\.EG\\.(TotalQuantity \\(Settings\\))"))
      res.df <- dplyr::bind_cols(res.df, intensities.df)
      col_info$eg_quantity <- colnames(intensities.df)
    }
    if ('eg_qvalue' %in% import_data) { # FIXME what's that?
      qvalues.df <- report.df %>% dplyr::select(matches("\\.EG\\.(Qvalue)"))
      res.df <- dplyr::bind_cols(res.df, qvalues.df)
      col_info$eg_qvalue <- colnames(qvalues.df)
    }
    if ('eg_snratio' %in% import_data) { # FIXME what's that?
      snratios.df <- report.df %>% dplyr::select(matches("\\.EG\\.(SignalToNoise)"))
      res.df <- dplyr::bind_cols(res.df, snratios.df)
      col_info$eg_snratio <- colnames(snratios.df)
    }
    if ('eg_normfactor' %in% import_data) { # FIXME what's that?
      normfactors.df <- report.df %>% dplyr::select(matches("\\EG\\.(NormalizationFactor)"))
      res.df <- dplyr::bind_cols(res.df, normfactors.df)
      col_info$eg_normfactor <- colnames(normfactors.df)
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
                        "TotalQuantity (Settings)" = "intensity",
                        "NormalizationFactor" = "msrun_normfactor",
                        "SignalToNoise" = "sn_ratio",
                        "Qvalue" = "qvalue"
                        )


pivot_longer.Spectronaut.Metrics <- function(msdata.wide, pkey, colgroup) {
    colgroups <- attr(msdata.wide, "column_groups")
    if (length(colgroup) == 1L) {
      cols <- colgroups[[colgroup]]
    } else {
      cols <- unlist(colgroups[colgroup])
      names(cols) <- NULL
    }
    if (is_empty(cols)) {
      stop("No columns found for column groups: ", str_c(colgroup, collapse=" "))
    }

    quantity_prespec_df <- tibble(.name = cols) %>%
      tidyr::extract(.name, c("rawfile_ix", "rawfile", ".value"), remove = FALSE,
                     "^\\[(\\d+)\\]\\s(.+)\\.(?:PG|PEP|EG)\\.(.+)$") %>%
      dplyr::mutate(.value = SpectronautMetrics[.value],
                    rawfile_ix = parse_integer(rawfile_ix))

    return (tidyr::pivot_longer_spec(
        dplyr::select(msdata.wide, !!pkey, !!cols),
        quantity_prespec_df))
}

#' @export
pivot_longer.Spectronaut.ProtgroupIntensities <- function(msdata.wide) {
    return (pivot_longer.Spectronaut.Metrics(msdata.wide, "protgroup_id", "quantity"))
}

#' @export
pivot_longer.Spectronaut.PepmodstateIntensities <- function(msdata.wide) {
    return (pivot_longer.Spectronaut.Metrics(msdata.wide, "pepmodstate_id", "quantity"))
}
