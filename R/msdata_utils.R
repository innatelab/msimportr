# utilities for preparing msdata list

uniprot_fasta_extract_ac <- function(fasta_headers) {
    str_replace_all(fasta_headers, "(^|;)(?:sp|tr)\\|([a-zA-Z0-9_-]+)\\|\\S+($|;)", "\\1\\2\\3")
}

agg_protgroup_col <- function(data) {
  if (is.character(data)) {
    str_c(unique(data), collapse=';')
  } else if (is.logical(data)) {
    any(data) # FIXME what if all?
  } else {
    error("don't know how to aggregate ", typeof(data), " type")
  }
}

#' @export
append_protgroups_info <- function(msdata, msdata_wide, proteins_info = NULL,
                                   fix_protein_info = TRUE,
                                   import_columns = NULL, verbose = TRUE) {
    if (verbose) message("Preparing msdata$protgroups data frame...")
    msdata_colgroups <- attr(msdata_wide, "column_groups")
    pg_df <- dplyr::select(msdata_wide, !!msdata_colgroups$protgroup) %>%
        dplyr::distinct() %>% dplyr::arrange(protgroup_id)
    # keep only protgroups with data
    intensity_dfnames <- str_subset(names(msdata), "intensities$")
    for (dfname in intensity_dfnames) {
        if (verbose) message("Contstraining msdata$protgroups to protgroups found in ", dfname)
        pg_df <- dplyr::semi_join(pg_df, msdata[[dfname]])
    }
    if (verbose) message(nrow(pg_df), " protgroup(s) found")
    if (fix_protein_info) {
        pg_df <- pg_df %>%
            dplyr::mutate(# fix ACs due to incorrect mqpar parse rules
                          protein_acs = uniprot_fasta_extract_ac(protein_acs),
                          majority_protein_acs = uniprot_fasta_extract_ac(majority_protein_acs))
    }
    protein2pg_df <- dplyr::bind_rows(
        dplyr::mutate(expand_protgroup_acs(pg_df, acs_col="majority_protein_acs"),
                      is_majority = TRUE) %>%
        dplyr::rename(protein_ac = majority_protein_ac),
        dplyr::mutate(expand_protgroup_acs(pg_df, acs_col="protein_acs"),
                      is_majority = FALSE)) %>%
        dplyr::group_by(protgroup_id, protein_ac) %>%
        dplyr::summarise(is_majority = any(is_majority)) %>%
        dplyr::ungroup()
    if (!is.null(proteins_info)) {
        proteins_df <- semi_join(proteins_info, protein2pg_df) %>%
            dplyr::left_join(dplyr::select(dplyr::filter(protein2pg_df, is_majority),
                                           protein_ac, protgroup_id))
        if (fix_protein_info) {
            # fix gene and protein names using FASTA information
            fix_cols <- c("protein_names" = "protein_name", "gene_names" = "gene_name")
        } else {
            fix_cols <- c()
        }
        agg_cols <- fix_cols
        if (!is.null(import_columns)) {
            agg_cols <- c(agg_cols, import_columns)
        }
        fixpg_df <- dplyr::left_join(dplyr::select(pg_df, protgroup_id),
                                     dplyr::group_by(proteins_df, protgroup_id) %>%
                                     dplyr::summarise_at(agg_cols, agg_protgroup_col) %>%
                                     dplyr::ungroup())
        if (any(fixpg_df$protgroup_id != pg_df$protgroup_id)) {
            warning("Incorrect fixed protgroups info, skipping")
        } else {
            for (fix_col in intersect(names(fix_cols), colnames(pg_df))) {
                pg_mask <- !is.na(fixpg_df[[fix_col]]) & is.na(pg_df[[fix_col]])
                if (any(pg_mask)) {
                    if (verbose) message("Fixing ", sum(pg_mask), " ", fix_col, " using provided proteins_info")
                    pg_df[pg_mask, fix_col] <- fixpg_df[pg_mask, fix_col]
                }
            }
            for (imp_col in import_columns) {
                if (verbose) message("Importing ", imp_col, " from proteins_info")
                pg_df[, imp_col] <- fixpg_df[[imp_col]]
            }
        }
    } else {
        proteins_df <- NULL
    }
    msdata$protgroups <- pg_df
    msdata$protein2protgroup <- protein2pg_df
    msdata$proteins <- proteins_df
    return (msdata)
}

#' @export
mschannel_statistics <- function(msdata) {
    res <- dplyr::left_join(tidyr::expand(dplyr::filter(msdata$protgroup_tagintensities, !is.na(protgroup_id)),
                                     protgroup_id, msrun, mstag),
                            dplyr::filter(msdata$protgroup_tagintensities, !is.na(protgroup_id)) %>%
                            dplyr::select(protgroup_id, msrun, mstag, intensity)) %>%
    dplyr::inner_join(dplyr::select(msdata$mschannels, msrun, mstag, mschannel) %>%
                      dplyr::distinct()) %>%
    #dplyr::group_by(protgroup_id, condition) %>% dplyr::filter(any(!is.na(intensity))) %>%
    dplyr::group_by(mschannel, mstag, msrun) %>%
    summarize(log2_intensity.mean = mean(log2(intensity[!is.na(intensity)])),
              log2_intensity.median = median(log2(intensity[!is.na(intensity)])),
              log2_intensity.sd = sd(log2(intensity[!is.na(intensity)])),
              n = n(),
              n_missing = sum(is.na(intensity))) %>%
    dplyr::ungroup()

    pg_idents_df <- msdata$protgroup_idents %||% msdata$protgroup_intensities %||% NULL
    if (!is.null(pg_idents_df)) {
        ident_stats <- dplyr::left_join(dplyr::filter(pg_idents_df, !is.na(protgroup_id)),
            dplyr::select(msdata$mschannels, msrun, msrun_mq) %>% dplyr::distinct()) %>%
            dplyr::group_by(msrun) %>%
        summarize(n_matching = sum(ident_type=="By matching", na.rm = TRUE),
                  n_msms = sum(ident_type=="By MS/MS", na.rm = TRUE)) %>%
        dplyr::ungroup()
        res <- left_join(res, ident_stats)
    } else {
        warning("No protgroup ident_type data found")
    }
    return (res)
}

#' @export
# FIXME more checks/control over the columns of intensities_df/stats_df
impute_intensities <- function(intensities_df, stats_df, log2_mean_offset=-1.8, log2_sd_scale=0.3){
    res <- dplyr::inner_join(intensities_df, stats_df) %>%
        dplyr::mutate(intensity_imputed = if_else(is.na(intensity),
                                                  2^(rnorm(n(), mean=log2_intensity.mean + log2_mean_offset,
                                                                sd=log2_intensity.sd * log2_sd_scale)),
                                                  intensity)) %>%
        dplyr::ungroup()
    # don't take stats_df columns along
    select(res, one_of(colnames(intensities_df)), intensity_imputed)
}
