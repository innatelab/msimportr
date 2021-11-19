# MaxQuant data import utils
#
# Author: astukalov
###############################################################################

require(dplyr)
require(readr)
require(tidyr)
require(stringr)

#' @importFrom rlang %||%
#' @importFrom dplyr tibble
#' @importFrom dplyr na_if starts_with matches row_number n n_distinct if_else case_when
#' @importFrom stringr str_split str_detect str_replace
NULL

#' @export
recombine_dlms <- function(dlms, sep=stringr::fixed(";")) {
    dlms_expanded <- unique(unlist(str_split(dlms, sep)))
    na_if(paste0(dlms_expanded[!is.na(dlms_expanded)], collapse=sep), "")
}

#' @export
expand_protgroups <- function(protgroup_ids, sep=stringr::fixed(";")) {
    protgroups2protgroup.list <- str_split(protgroup_ids, sep)
    tibble(protgroup_ids = rep.int(protgroup_ids, sapply(protgroups2protgroup.list, length)),
           protgroup_id = as.integer(unlist(protgroups2protgroup.list)))
}

#' @export
expand_sites <- function(sites_df, by=c("protgroup_id", "protein_ac"),
                         keep_cols=c("site_id"), sep=stringr::fixed(";"))
{
  exp_by <- match.arg(by)
  collapsed_by <- paste0(exp_by, "s")
  exp_list <- str_split(sites_df[[collapsed_by]], sep)
  exp_lengths <- sapply(exp_list, length)
  res <- sites_df[rep.int(1:nrow(sites_df), exp_lengths), c(collapsed_by, keep_cols)]
  res[[exp_by]] <- unlist(exp_list)
  if ("site_ids" %in% colnames(sites_df)) {
    site_list <- str_split(sites_df$site_ids, sep)
    site_lengths <- sapply(site_list, length)
    if (any(site_lengths != exp_lengths)) stop("site lengths mismatch")
    res$site_id <- unlist(site_list)
  }
  poses_col <- switch(exp_by,
    protgroup_id = "positions",
    protein_ac = "positions_in_proteins"
  )
  if (poses_col %in% colnames(sites_df)) {
    pos_list <- str_split(sites_df[[poses_col]], sep)
    pos_lengths <- sapply(pos_list, length)
    if (any(pos_lengths != exp_lengths)) stop("pos lengths mismatch")
    res$position <- as.integer(unlist(pos_list))
  } else {
    warning("No position information column (", poses_col, ") found")
  }
  if ("protgroup_id" %in% colnames(res)) {
    res$protgroup_id <- as.integer(res$protgroup_id)
  }
  return(res)
}

zero2na <- function(x) if_else(x == 0.0, NA_real_, x)

#' @export
expand_peptides <- function(peptides_df, by=c("protgroup_id", "protein_ac"),
                            keep_cols=c("start_pos", "end_pos"), sep=";")
{
    exp_by <- match.arg(by)
    expand_collapsed(peptides_df, paste0(exp_by,"s"), exp_by,
                     c("peptide_id", keep_cols), sep=sep)
}

#' @export
expand_pepmods <- function(pepmods_df, by=c("protgroup_id", "protein_ac"),
                           keep_cols=c("peptide_id"), sep=";")
{
  exp_by <- match.arg(by)
  expand_collapsed(pepmods_df, paste0(exp_by,"s"), exp_by,
                   c("pepmod_id", keep_cols), sep=sep)
}

selectUniprotACs <- function(acs, valid_acs)
{
    acs_noiso <- stripUniprotIsoform(acs)
    if (!is.null(valid_acs)) {
        acs_noiso <- intersect(acs_noiso, valid_acs)
    }
    if (length(acs_noiso) == 1L) {
        return (acs_noiso)
    } else if (length(acs_noiso) > 1L) {
        # return the first "classical" UniProt AC starting with O, P or Q
        is_classical_uprot <- str_detect(acs_noiso, '^[POQ]')
        return (acs_noiso[[ order(!is_classical_uprot, acs_noiso)[[1]] ]])
    } else {return (NA_character_)}
}

MaxQuant_NAs <- c("", "NA", "NaN", "n. def.", "n.def.")

#' @export
read.MaxQuant <- function(filename, layout = c("wide", "long"),
                          row_id_cols = c('protein_ac_noiso'),
                          protein_ac_cols, protein_info = NULL,
                          measures.regex)
{
    res <- read.table(filename, sep = '\t', header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    #print(colnames(res))
    measures <- gsub('\\\\', '', measures.regex)
    measures.fixed <- measures %>% gsub('\\s*\\[\\%\\]', '', .) %>% gsub('/', '', .) %>% gsub('\\+', 'and', .) %>% gsub('\\s', '_', .)
    names(measures.fixed) <- measures

    # separate measure from SILAC label and from the run label and fix the total/summary names
    colnames(res) <- colnames(res) %>%
        sub('^Sequence coverage (.+) \\[\\%\\]$', 'Sequence coverage [%] \\1', .) %>% # fix the position of [%] to make the naming scheme consitent
        sub(paste0('^(',paste0(measures.regex, collapse='|'),')(\\s([HLM]))?(\\s([^HLM].*))?$'), '\\1.\\3__\\5', .) %>%
        sub('\\.__', '.Sum__', .) %>% sub('\\__$', '\\__Everything', .)
    #print(colnames(res))

    # fix UniProt ACs
    for (i in seq_along(protein_ac_cols)) {
        if (!(protein_ac_cols[[i]] %in% colnames(res))) {
            warning('Protein AC column `', protein_ac_cols[[i]], '` not found')
            next
        }
        fixed_acs <- sapply(strsplit(res[[ protein_ac_cols[[i]] ]], ';', fixed = TRUE), selectUniprotACs,
                            if (is.null(protein_info)) NULL else protein_info$protein_ac_noiso)
        na_mask <- is.na(fixed_acs)
        if ( any( na_mask ) ) {
            warning(sum(na_mask), ' records have no valid protein AC, removed: ',
                    paste0(res[[ protein_ac_cols[[i]] ]][na_mask], collapse = ' '))
            res <- res[!na_mask, , drop = FALSE]
        }
        res[, names(protein_ac_cols)[[i]]] <- fixed_acs[!na_mask]
    }

    quant.columns.list <- lapply(measures.regex, function(mes) grep(paste0('^', mes, '\\.'), colnames(res), value = TRUE))
    names(quant.columns.list) <- measures
    quant.columns <- unlist(quant.columns.list)
    channels <- sort(unique(gsub( '^[^.]+\\.', '', quant.columns)))
    channels.df <- cbind(mschannel = channels, do.call(rbind, strsplit(channels, '__', fixed = TRUE))) %>%
        as.data.frame(stringsAsFactors = FALSE)
    colnames(channels.df) <- c('mschannel', 'mstag', 'msrun')
    channels.df <- dplyr::mutate(channels.df,
        mstag = factor(mstag, levels = intersect(c('H', 'M', 'L', 'Sum'), unique(mstag))),
        mstag_type = if_else(mstag %in% c('H', 'M', 'L'), 'SILAC', 'label_free'))

    quant.columns.df <- expand.grid(measure = measures, mschannel = channels,
                                    stringsAsFactors = FALSE) %>%
       mutate(colname = paste0(measure, '.', channel),
              measure_fixed = measures.fixed[measure],
              colname_fixed = paste0(measure_fixed, '.', mschannel)) %>%
       inner_join(channels.df)

    #print(str(res))
    #print( quant.columns.df )
    res <- res %>% mutate_at(quant.columns, as.numeric) %>%
        group_by(!!row_id_cols) %>%
        summarise_at(quant.columns, max) %>%
        dplyr::ungroup()

    if (match.arg(layout) == 'long') {
        quant.columns.df$exists <- quant.columns.df$colname %in% colnames(res)
        res[, subset(quant.columns.df,!exists)$colname ] <- NA # add missing column names to make it reshapeable
        rename_cols <- quant.columns.df$colname
        names(rename_cols) <- quant.columns.df$colname_fixed
        res <- rename(res, !!rename_cols)
        quant.columns.full_list <- lapply(measures.fixed, function(mes) subset(quant.columns.df,measure_fixed==mes)$colname_fixed)
        names(quant.columns.full_list) <- measures.fixed
        ratio_cols <- str_c("ratio.", ratio_columns_sel.df$suffix)
        res <- tidyr::pivot_longer(select(res, !!row_id_cols, !!quant.columns.full_list),
                                   cols=cols(ratio_cols),
                                   names_to=c('.value', 'mschannel'),
                                   names_pattern="(.+)\\.(.+)")
        # v.names = names(quant.columns.full_list)
        #times = sort(unique(quant.columns.df$mschannel))
        res <- inner_join(res, channels.df) %>% mutate(mschannel = NULL) # add mstags and msrun info
        # fix NA
        for (mes in measures.fixed) {
            res[!is.na(res[[mes]]) & res[[mes]] == 0, mes] <- NA
        }
    }
    if (!is.null(protein_info)) {
        prot_info_cols <- c('protein_label','genename','molweight','contaminant_family')
        join_by <- c(protein_ac_noiso = 'protein_ac_noiso')
        names(join_by) <- names(protein_ac_cols)[[1]]
        res <- left_join(res, dplyr::select(protein_info, protein_ac_noiso, !!prot_info_cols),
                         by = join_by)
    }
    return ( res )
}

#' @export
read.MaxQuant.ProteinGroups <- function(folder_path, file_name = 'proteinGroups.txt',
                                        protein_info = NULL, layout = c("wide", "long"),
                                        nrows = Inf, import_data = c(), guess_max = min(10000L, nrows))
{
    proteinGroups.df <- readr::read_tsv(file.path(folder_path, file_name), n_max = nrows,
                                        col_types = readr::cols(`Fasta headers` = "c", `id` = "i"),
                                        na = MaxQuant_NAs, guess_max = guess_max)
    col_renames <- c("protgroup_id" = "id", "protein_acs" = "Protein IDs",
                     "majority_protein_acs" = "Majority protein IDs",
                     "gene_names" = "Gene names", "protein_names" = "Protein names",
                     "fasta_headers" = "Fasta headers", "n_proteins" = "Number of proteins",
                     "npeptides" = "Peptide counts (all)",
                     "npeptides_unique_razor" = "Peptide counts (razor+unique)",
                     "npeptides_unique" = "Peptide counts (unique)",
                     "npeptides_mutated" = "Mutated peptide count",
                     "score" = "Score", "q_value" = "Q-value",
                     "seqlen" = "Sequence length", "seqlens" = "Sequence lengths",
                     "mol_weight_kDa" = "Mol. weight [kDa]",
                     "seqcov" = "Sequence coverage [%]","unique_razor_seqcov" = "Unique + razor sequence coverage [%]",
                     "mutations" = "Mutation names",
                     "is_contaminant" = "Potential contaminant", "is_reverse" = "Reverse")
    res.df <- dplyr::select(proteinGroups.df, !!col_renames[col_renames %in% colnames(proteinGroups.df)]) %>%
        dplyr::mutate(is_contaminant = replace_na(is_contaminant, "") == '+',
                      is_reverse = replace_na(is_reverse, "") == '+')
    col_info <- list(protgroup = colnames(res.df))
    if ('intensity' %in% import_data) {
        intensities.df <- proteinGroups.df %>% dplyr::select(starts_with("Intensity")) %>%
            dplyr::rename_all(~str_replace(.x, "^Intensity\\s([LMH](\\s|$))", "Intensity.\\1")) %>%
            dplyr::rename_all(~str_replace(.x, "^Intensity(\\s|$)", "Intensity.Sum\\1")) %>%
            dplyr::mutate_all(zero2na)
        res.df <- dplyr::bind_cols(res.df, intensities.df)
        col_info$intensity <- colnames(intensities.df)
    }
    if ('LFQ' %in% import_data) {
        lfq.df <- proteinGroups.df %>% dplyr::select(starts_with("LFQ intensity")) %>%
            dplyr::rename_all(~str_replace(.x, "^LFQ intensity\\s([LMH](\\s|$))", "LFQ_Intensity.\\1")) %>%
            dplyr::rename_all(~str_replace(.x, "^LFQ intensity(\\s|$)", "LFQ_Intensity.Sum\\1")) %>%
            dplyr::mutate_all(zero2na)
        res.df <- dplyr::bind_cols(res.df, lfq.df)
        col_info$LFQ <- colnames(lfq.df)
    }
    if ('iBAQ' %in% import_data) {
        ibaq.df <- proteinGroups.df %>% dplyr::select(starts_with("iBAQ")) %>%
            dplyr::rename_all(~str_replace(.x, "^iBAQ\\s([LMH](\\s|$))", "iBAQ.\\1")) %>%
            dplyr::rename_all(~str_replace(.x, "^iBAQ(\\s|$)", "iBAQ.Sum\\1")) %>%
            dplyr::mutate_all(zero2na)
        res.df <- dplyr::bind_cols(res.df, ibaq.df)
        col_info$iBAQ <- colnames(ibaq.df)
    }
    if ('ratio' %in% import_data) {
        ratios.df <- proteinGroups.df %>% dplyr::select(starts_with("Ratio")) %>%
            dplyr::rename_all(~str_replace(.x, "^Ratio\\s", "Ratio."))
        res.df <- dplyr::bind_cols(res.df, ratios.df)
        col_info$ratio <- colnames(ratios.df)
    }
    if ('ident_type' %in% import_data) {
        ident_types.df <- proteinGroups.df %>% dplyr::select(starts_with("Identification type")) %>%
          dplyr::rename_all(~str_replace(.x, "^Identification\\stype\\s", "ident_type.")) %>%
          dplyr::mutate_all(~factor(.x, levels=c("By matching", "By MS/MS")))
        res.df <- dplyr::bind_cols(res.df, ident_types.df)
        col_info$ident_type <- colnames(ident_types.df)
    }
    if ('ms2_count' %in% import_data) {
      ms2_counts.df <- proteinGroups.df %>% dplyr::select(starts_with("MS/MS count ")) %>%
        dplyr::rename_all(~str_replace(.x, "^MS/MS\\scount\\s", "ms2_count."))
      res.df <- dplyr::bind_cols(res.df, ms2_counts.df)
      col_info$ms2_count <- colnames(ms2_counts.df)
    }
    attr(res.df, "column_groups") <- col_info
    return (res.df)
}

#' @export
read.MaxQuant.Peptides <- function(folder_path, file_name = 'peptides.txt',
                                   import_data = c(), nrows = Inf)
{
    peptides.df <- readr::read_tsv(file.path(folder_path, file_name), n_max = nrows,
                                   col_types = readr::cols(
                                        `Proteins` = "c",
                                        `Protein group IDs` = "c",
                                        `Mod. peptide IDs` = "c",
                                        `Leading razor protein` = "c",
                                        `Charges` = "c",
                                        `Reverse` = "c", `Potential contaminant` = "c",
                                        `Ratio H/L variability [%]` = "n"),
                                   na = MaxQuant_NAs, guess_max = 20000L)
    col_renames <- c(peptide_seq = "Sequence", peptide_len = "Length",
                     n_miscleaved = "Missed cleavages", n_miscleaved = "Missed cleavages (Trypsin/P;LysC/P)", n_miscleaved_lysc = "Missed cleavages (LysC/P)",
                     nterm_window = "N-term cleavage window", cterm_window = "C-term cleavage window",
                     aa_before = "Amino acid before", aa_after = "Amino acid after",
                     aa_first = "First amino acid", aa_last = "Last amino acid",
                     protgroup_ids = "Protein group IDs",
                     pepmod_ids = "Mod. peptide IDs",
                     protein_acs = "Proteins", lead_razor_protein_ac = "Leading razor protein",
                     start_pos = "Start position", end_pos = "End position",
                     mass = "Mass", charges = "Charges", score = "Score", PEP = "PEP",
                     is_mutated = "Mutated", mutations = "Mutation names",
                     is_reverse = "Reverse", is_contaminant = "Potential contaminant",
                     is_group_unique = "Unique (Groups)",
                     is_protein_unique = "Unique (Proteins)")
    res.df <- dplyr::select(peptides.df, !!col_renames[col_renames %in% colnames(peptides.df)]) %>%    #message(paste0(colnames(peptides.df), collapse='\n'))
        mutate(peptide_id = row_number()-1L,
               is_reverse = replace_na(is_reverse, "") == "+",
               is_contaminant = replace_na(is_contaminant, "") == "+",
               is_shared_by_groups = !(replace_na(is_group_unique, "") %in% c('yes', '+')),
               is_shared = is_shared_by_groups,
               is_shared_by_proteins = !(replace_na(is_protein_unique, "") %in% c('yes', '+'))) %>%
        select(-is_group_unique, -is_protein_unique)
    col_info <- list(peptide = colnames(res.df))
    if ('intensity' %in% import_data) {
        intensities.df <- dplyr::select(peptides.df, starts_with("Intensity")) %>%
            dplyr::rename_all(~str_replace(.x, "^Intensity\\s([LMH](\\s|$))", "Intensity.\\1")) %>%
            dplyr::rename_all(~str_replace(.x, "^Intensity(\\s|$)", "Intensity.Sum\\1")) %>%
            dplyr::mutate_all(zero2na)
        res.df <- dplyr::bind_cols(res.df, intensities.df)
        col_info$intensity <- colnames(intensities.df)
    }
    if ('aa_stats' %in% import_data) {
       aa_stats.df <- dplyr::select(peptides.df, matches("^[A-Z] Count$")) %>%
         dplyr::rename_all(~str_replace(.x, "([A-Z]) Count", "count_\\1"))
       res.df <- bind_cols(res.df, aa_stats.df)
       col_info$aa_stats <- colnames(aa_stats.df)
    }
    if ('ident_type' %in% import_data) {
        ident_types.df <- dplyr::select(peptides.df, starts_with("Identification type")) %>%
          dplyr::rename_all(~str_replace(.x, "^Identification\\stype\\s", "ident_type.")) %>%
          dplyr::mutate_all(~factor(.x, levels=c("By matching", "By MS/MS")))
        res.df <- dplyr::bind_cols(res.df, ident_types.df)
        col_info$ident_type <- colnames(ident_types.df)
    }
    attr(res.df, "column_groups") <- col_info
    return (res.df)
}

# read allPeptides.txt table
#' @export
read.MaxQuant.AllPeptides <- function( folder_path, file_name = 'allPeptides.txt', nrows = Inf )
{
    allPeptides.df <- readr::read_tsv(file.path(folder_path, file_name ), n_max = nrows,
                                      col_types = readr::cols(`Type` = "c",
                                                              `Raw file` = "c",
                                                              `Resolution` = "n"
                                      ),
                                      na = MaxQuant_NAs, guess_max = 20000L)
    allPeptides.df <- dplyr::rename(allPeptides.df,
                                    pep_type = Type,
                                    rawfile = `Raw file`,
                                    protein_acs = `Proteins`,
                                    dp_protein_acs = `DP Proteins`,
                                    #protgroup_ids = `Protein group IDs`,
                                    #fasta_headers = `Fasta headers`,
                                    dp_modif = `DP Modification`,
                                    dp_mass_diff = `DP Mass Difference`) %>%
        dplyr::mutate(pep_type = factor(pep_type),
                      rawfile = factor(rawfile),
                      #protgroup_ids = factor(protgroup_ids),
                      #fasta_headers = factor(fasta_headers),
                      is_dp = !is.na(dp_mass_diff),
                      dp_modif = factor(dp_modif),
                      dp_mass_delta = as.integer(dp_mass_diff*100)/100)
    # remove less useful memory-hungry columns
    #allPeptides.df$Intensities <- NULL
    #allPeptides.df$dp_mass_delta <- as.integer(allPeptides.df$`DP Mass Difference`*100)/100
    return(allPeptides.df)
}

#' @export
read.MaxQuant.Sites <- function(folder_path, file_name, nrows = Inf, modif = "Phospho (STY)",
                                import_data = c())
{
    data.df <- readr::read_tsv(file.path(folder_path, file_name),
                        col_names = TRUE, n_max = nrows,
                        col_types = readr::cols(`Protein group IDs` = readr::col_character(),
						                        `Reverse` = "c", `Potential contaminant` = "c",
                                                `id` = "i",
                                                .default = readr::col_guess()),
                        na = MaxQuant_NAs, guess_max = 20000L)
    sites.df <- data.df %>%
        dplyr::mutate(is_contaminant = replace_na(`Potential contaminant`, "") == '+',
                      is_reverse = replace_na(`Reverse`, "") == '+') %>%
        dplyr::select(site_id = id,
                      protgroup_ids = `Protein group IDs`,
                      leading_protein_acs = `Leading proteins`,
                      protein_ac = `Protein`,
                      protein_acs = `Proteins`,
                      loc_prob = `Localization prob`,
                      delta_score = `Delta score`,
                      score = `Score`,
                      loc_score = `Score for localization`,
                      peptide_ids = `Peptide IDs`,
                      pepmod_ids = `Mod. peptide IDs`,
                      positions = `Positions`,
                      position = `Position`,
                      positions_in_proteins = `Positions within proteins`,
                      aa = `Amino acid`,
                      seq_context = `Sequence window`,
                      modif_context = `Modification window`,
                      gene_names = `Gene names`,
                      protein_names = `Protein names`,
                      fasta_headers = `Fasta headers`,
                      charge = `Charge`,
                      mass_error_ppm = `Mass error [ppm]`,
                      is_contaminant, is_reverse)
    sites.df["n_of_modifs"] <- data.df[[paste0("Number of ", modif)]]
    col_info <- list( site = colnames(sites.df) )
    res.df <- sites.df
    if ('intensity' %in% import_data) {
        intensities.df <- data.df %>% dplyr::select(starts_with("Intensity")) %>%
            dplyr::rename_all(~str_replace(.x, "^Intensity\\s([LMH](\\s|_|$))", "Intensity.\\1")) %>%
            dplyr::rename_all(~str_replace(.x, "^Intensity(\\s|_|$)", "Intensity.Sum\\1")) %>%
            dplyr::rename_all(~str_replace(.x, "^Intensity.Sum(.+)(___\\d+)$", "Intensity.Sum\\2\\1")) %>%
            dplyr::mutate_all(~zero2na(as.numeric(.x)))
        res.df <- dplyr::bind_cols(res.df, intensities.df)
        col_info$intensity <- colnames(intensities.df)
    }
    if ('occupancy' %in% import_data) {
        occupancies.df <- data.df %>% dplyr::select(starts_with("Occupancy")) %>%
            dplyr::rename_all(~str_replace(.x, "^(Occupancy\\s|Occupancy ratio|Occupancy error scale\\s)([LMH](\\s|$))", "\\1.\\2")) %>%
            dplyr::rename_all(~str_replace(.x, "^Occupancy ratio", "occupancy_ratio")) %>%
            dplyr::rename_all(~str_replace(.x, "^Occupancy error scale", "occupancy_error_scale")) %>%
            dplyr::mutate_all(as.numeric)
        res.df <- dplyr::bind_cols(res.df, occupancies.df)
        col_info$occupancy <- colnames(occupancies.df)
    }
    if ('ratio' %in% import_data) {
        ratios.df <- data.df %>% dplyr::select(starts_with("Ratio mode/base")) %>%
            dplyr::rename_all(~str_replace(.x, "^Ratio mod/base\\s", "ratio_mod2base.")) %>%
            dplyr::mutate_all(as.numeric)
        res.df <- dplyr::bind_cols(res.df, ratios.df)
        col_info$ratio <- colnames(ratios.df)
    }
    attr(res.df, "column_groups") <- col_info
    return (res.df)
}

backquote <- function( vars ) {
    res <- paste0( "`", vars, "`" )
    names(res) <- names(vars)
    return (res)
}

MaxQuant_IdentTypes <- c("ISO-MSMS", "MULTI-MSMS", "MSMS", "MULTI-SECPEP", "MULTI-MATCH", "MULTI-MATCH-MSMS")

read.MaxQuant.Evidence_internal <- function(folder_path, file_name = 'evidence.txt',
                                            nrows = Inf, guess_max = min(20000L, nrows)) {
    message('Reading evidence table...')
    evidence.df <- readr::read_tsv(file.path(folder_path, file_name), col_names = TRUE, n_max = nrows,
                            col_types = readr::cols(Resolution = 'n', Fraction = 'i',
                                            `Protein group IDs` = 'c',
                                            `Oxidation (M) site IDs` = 'c',
                                            `Met->AHA site IDs` = 'c',
                                            `Peptide ID` = 'i', `Mod. peptide ID` = 'i', `Charge` = 'i',
                                             Type = 'c', `Raw file` = 'c', Experiment = 'c',
                                             Modifications = 'c', `Labeling State` = 'c',
                                            `Reverse` = "c", `Potential contaminant` = "c",
                                             .default = readr::col_guess()),
                            na = MaxQuant_NAs, guess_max = guess_max)

    message( 'Renaming and converting evidence table columns...' )
    # fix column names from different MQ versions
    col_renames <- c(evidence_id = 'id',
                     pepmod_id = 'Mod. peptide ID',
                     rawfile = "Raw file", msexperiment_mq = "Experiment", msfraction_mq = "Fraction",
                     ident_type = "Type", label_state = 'Labeling State',
                     mass_error_ppm = "Mass Error [ppm]", mass_error_da = "Mass Error [Da]",
                     mass_error_ppm = "Mass error [ppm]", mass_error_da = "Mass error [Da]",
                     uncalib_mass_error_ppm = "Uncalibrated Mass Error [ppm]",
                     uncalib_mass_error_da = "Uncalibrated Mass Error [Da]",
                     peptide_seq = 'Sequence', peptide_len = 'Length', count_K = 'K Count', count_R = 'R Count', modifs = 'Modifications',
                     n_miscleaved = "Missed cleavages", n_miscleaved = "Missed cleavages (Trypsin/P;LysC/P)", n_miscleaved_lysc = "Missed cleavages (LysC/P)",
                     pepmod_seq = 'Modified sequence', charge = 'Charge',
                     protein_acs = 'Proteins', lead_protein_acs = 'Leading proteins', lead_razor_protein_ac = 'Leading razor protein',
                     gene_names = 'Gene Names', gene_names = 'Gene names',
                     protein_names = 'Protein Names', protein_names = 'Protein names',
                     fasta_headers = 'Fasta headers',
                     is_reverse = 'Reverse',
                     is_contaminant = "Potential contaminant", is_contaminant = 'Contaminant',
                     rt_len = 'Retention length', rt_avg = 'Calibrated retention time', rt_start = 'Calibrated retention time start', rt_finish = 'Calibrated retention time finish',
                     protgroup_ids = 'Protein group IDs', # IDs between different folders do not match
                     # Carbamidomethyl (C) site IDs', Oxidation (M) site IDs
                     peptide_id = 'Peptide ID',
                     mz = 'm/z', ms2_mz = 'MS/MS m/z', mass_da = 'Mass', resolution = 'Resolution',
                     delta_mz_da = "Uncalibrated - Calibrated m/z [Da]",
                     delta_mz_ppm = "Uncalibrated - Calibrated m/z [ppm]"
    )
    if (any(str_detect(colnames(evidence.df), "^Intensity \\S+"))) {
        # intensity is an aggregation column of labeled intensities
        col_renames = c(col_renames, "Intensity Sum" = "Intensity")
    } else {
        # label-free data, rename column so it gets the tag
        col_renames = c(col_renames, "Intensity F" = "Intensity")
    }
    col_renames <- col_renames[col_renames %in% colnames(evidence.df)]
    if (length(col_renames) > 0) {
        evidence.df <- dplyr::rename(evidence.df, !!col_renames)
    }
    evidence.df <- dplyr::mutate(evidence.df,
                                 msexperiment_mq = factor(msexperiment_mq),
                                 rawfile = factor(rawfile),
                                 ident_type = factor(ident_type, levels=MaxQuant_IdentTypes),
                                 is_contaminant = replace_na(is_contaminant, "") == '+',
                                 is_reverse = replace_na(is_reverse, "") == '+')
    return ( evidence.df )
}

process.MaxQuant.Evidence <- function( evidence.df,
                                       import_data = c("intensity"),
                                       mschannel_annotate.f = NULL,
                                       na_weight = 1E-5, min_intensity = 1E+3,
                                       multidplyr_cluster = NULL )
{
    complete_names <- function( vars ) {
        res <- vars
        res_names <- names(vars)
        res_names[res_names==""] <- vars[res_names==""]
        names(res) <- res_names
        return (res)
    }
    message("Extracting MS runs and MS channels info...")
    intensity_columns.df <- tidyr::expand_grid(measure = 'intensity',
            mstag = factor(c("H", "M", "L", "F", 'Sum'),
                           levels = c("H", "M", "L", "F", 'Sum'))) %>%
        mutate(old_name = paste0('Intensity ', mstag),
               new_name = paste0(measure, '.', mstag),
               quant_type = case_when(mstag %in% c('H','M','L') ~ 'SILAC',
                                      mstag == 'F' ~ 'label_free',
                                      mstag == 'Sum' ~ 'aggregate',
                                      TRUE ~ NA_character_)) %>%
        dplyr::filter(old_name %in% colnames(evidence.df)) %>%
        dplyr::mutate(quant_type = factor(quant_type, levels=c('label_free', 'SILAC', 'aggregate'))) %>%
        dplyr::arrange(mstag) %>%
        # restrict to the labels actually used
        dplyr::mutate(mstag = factor(mstag, levels=as.character(unique(mstag))))
    rawfiles.df <- dplyr::select(evidence.df, rawfile, msexperiment_mq, msfraction_mq) %>% dplyr::distinct()
    mschannels.df <- tidyr::expand_grid(rawfile = unique(evidence.df$rawfile),
                                        mstag = unique(intensity_columns.df$mstag)) %>%
        dplyr::inner_join(rawfiles.df, by="rawfile") %>%
        dplyr::inner_join(dplyr::select(intensity_columns.df, mstag, quant_type) %>% dplyr::distinct(),
                          by="mstag")
    message("Data contains ", nrow(mschannels.df), " mschannel(s)")
    if (!is.null(mschannel_annotate.f)) {
        message("Requesting user-defined MS channel annotations...")
        # get user-defined mschannel annotations
        annot.df <- mschannel_annotate.f(mschannels.df)
        checkmate::assert_data_frame(annot.df)
        checkmate::assert_names(colnames(annot.df), must.include = "rawfile")
        checkmate::assert_set_equal(as.character(annot.df$rawfile),
                                    as.character(mschannels.df$rawfile))
        if (rlang::has_name(attributes(annot.df), "column_scopes")) {
            message("  * detected user-provided MS channel column scopes")
            annot_col_scopes <- attr(annot.df, "column_scopes", exact=TRUE)
            checkmate::assert_character(annot_col_scopes, any.missing=FALSE, names="unique")
            checkmate::assert_subset(annot_col_scopes, choices=c("msexperiment", "msrun", "mschannel"))
            checkmate::assert_subset(names(annot_col_scopes), colnames(annot.df))
        } else {
            annot_col_scopes <- c()
        }
        if (rlang::has_name(colnames(annot.df), "msexperiment_mq")) {
            checkmate::assert_set_equal(as.character(annot.df$msexperiment_mq),
                                        as.character(mschannels.df$msexperiment_mq))
        }
        if (nrow(annot.df) != nrow(mschannels.df)) {
            stop("User-specified mschannel annotations has incorrect number of rows")
        }
        mschannels_annot.df <- dplyr::left_join(mschannels.df, annot.df)
        if (nrow(mschannels_annot.df) != nrow(mschannels.df)) {
            stop("User-specified mschannel annotations do not correctly match mschannels")
        }
        checkmate::assert_set_equal(as.character(mschannels_annot.df$rawfile),
                                    as.character(mschannels.df$rawfile))
        mschannels.df <- mschannels_annot.df
        if (rlang::has_name(mschannels.df, 'is_skipped')) {
            message('  * skipping ', sum(mschannels.df$is_skipped, na.rm=TRUE),
                    ' of ', nrow(mschannels.df), ' mschannels(s) excluded by user')
            mschannels.df <- dplyr::filter(mschannels.df, !dplyr::coalesce(is_skipped, FALSE))
            # drop unused levels
            mschannels.df <- dplyr::mutate_at(mschannels.df,
                    vars(any_of(c("msexperiment_mq", "msexperiment", "msrun", "mschannel", "rawfile"))),
                    ~factor(., levels=unique(as.character(.))))
            evidence.df <- dplyr::filter(evidence.df, rawfile %in% mschannels.df$rawfile)
            evidence.df <- dplyr::mutate(evidence.df,
                    rawfile = factor(rawfile, levels=levels(mschannels.df$rawfile)),
                    msexperiment_mq = factor(msexperiment_mq, levels=levels(mschannels.df$msexperiment_mq)))
            if (nrow(evidence.df) == 0) {
                stop("Empty evidence frame after unused MS filtering")
            }
        }
        if (rlang::has_name(mschannels.df, "msrun")) {
            if (n_distinct(mschannels.df$msrun) != n_distinct(mschannels.df$rawfile)) {
                stop("User-specified msrun IDs in mschannel annotations do not correctly match rawfiles")
            }
        }
    } else {
        annot_col_scopes <- c()
    }
    if (rlang::has_name(mschannels.df, "msexperiment")) {
        message('Overriding MaxQuant experiment IDs (msexperiment)')
    } else {
        mschannels.df$msexperiment <- mschannels.df$msexperiment_mq
    }
    mschannel_cols <- "msexperiment"
    if (rlang::has_name(mschannels.df, "msfraction")) {
        message('Overriding MaxQuant fraction IDs (msfraction)')
        mschannel_cols <- c(mschannel_cols, "msfraction")
    } else {
        mschannels.df <- dplyr::mutate(mschannels.df, msfraction = msfraction_mq)
    }
    if (any(intensity_columns.df$quant_type == 'SILAC')) {
        message('Evidence table contains SILAC-labeled data')
        mschannel_cols <- c(mschannel_cols, "mstag")
    }
    if (!rlang::has_name(mschannels.df, "msprotocol")) {
        message("No MS protocols specified")
        mschannels.df <- dplyr::mutate(mschannels.df, msprotocol = NA_integer_)
    }
    if (!rlang::has_name(mschannels.df, "msrun")) {
        message("No MS runs specified, Generating IDs")
        if (!rlang::has_name(mschannels.df, "msfraction")) {
            message("Assuming msrun=msexperiment")
            mschannels.df <- dplyr::mutate(mschannels.df, msrun = msexperiment)
        } else {
            message("Assuming msrun=msexperiment X msfraction")
            mschannels.df <- dplyr::mutate(mschannels.df,
                `__msfraction_sym__` = paste0("F", msfraction),
                msrun = interaction(msexperiment, `__msfraction_sym__`, drop=TRUE, lex.order=TRUE, sep='_'),
                `__msfraction_sym__` = NULL)
        }
    }
    if (!rlang::has_name(mschannels.df, "mschannel")) {
        message("Generating MS channel IDs from ", paste0(mschannel_cols, collapse=", "))
        if ("mstag" %in% mschannel_cols) {
            mschannels.df <- dplyr::mutate(mschannels.df,
                mschannel = interaction(msrun, mstag, drop=TRUE, lex.order=TRUE, sep='_'))
        } else {
            mschannels.df <- dplyr::mutate(mschannels.df, mschannel = msrun)
        }
    }
    if (n_distinct(mschannels.df$mschannel) != nrow(mschannels.df)) {
        stop('mschannel ids are not unique')
    }
    annot_col_scopes[['msexperiment']] <- "msexperiment"
    annot_col_scopes[['msexperiment_mq']] <- "msexperiment"
    annot_col_scopes[['msprotocol']] <- "msexperiment"
    annot_col_scopes[['rawfile']] <- "msrun"
    annot_col_scopes[['msrun']] <- "msrun"
    annot_col_scopes[['msfraction']] <- "msrun"
    annot_col_scopes[['msfraction_mq']] <- "msrun"
    annot_col_scopes[['mschannel']] <- "mschannel"
    annot_col_scopes[['mstag']] <- "mschannel"
    annot_col_scopes[['quant_type']] <- "mschannel"
    mschannels.df <- dplyr::mutate(mschannels.df, quant_type = factor(quant_type))
    msruns.df <- dplyr::select_at(mschannels.df, names(annot_col_scopes)[annot_col_scopes %in% c('msrun', 'msexperiment')]) %>%
                 dplyr::distinct()
    msexperiments.df <- dplyr::select_at(msruns.df, names(annot_col_scopes)[annot_col_scopes == 'msexperiment']) %>% dplyr::distinct()
    # NOTE?: the same mod. peptide Id can have multiple mod. sequences (mod at different poses)
    message('Enumerating pepmod states and summarizing channel intensities...')
    intensities.mtx <- as.matrix(evidence.df[,dplyr::filter(intensity_columns.df, quant_type != 'aggregate')$old_name])
    intens_cols <- rlang::set_names(intensity_columns.df$old_name, intensity_columns.df$new_name)
    evidence.df <- dplyr::mutate(evidence.df,
                          pepmod_id_mq = pepmod_id,
                          # MaxQuant has strange way of indexing modified peptides with PTM
                          # that don't take into account the position of the modification
                          # redefine the pepmod_id based on the pepmod_seq
                          pepmod_id = match(pepmod_seq, unique(pepmod_seq)) - 1L,
                          n_quants = Matrix::rowSums(!is.na(intensities.mtx) & intensities.mtx > 0.0),
                          is_full_quant = !is.na(Matrix::rowSums(intensities.mtx > 0.0))) %>%
        dplyr::rename(!!intens_cols) %>%
        dplyr::inner_join(dplyr::select(msruns.df, msexperiment, msrun, rawfile, msfraction, msprotocol),
                          by=c("rawfile"))

    pepmodstate_cols <- c("pepmod_id", "charge")
    if (any(!is.na(mschannels.df$msfraction))) {
        message("MS data contains MS fractions, pepmodstate in different fractions get unique IDs")
        pepmodstate_cols <- c(pepmodstate_cols, "msfraction")
    }
    message('Extracting pepmod states (', paste0(pepmodstate_cols, collapse=' X '), ')...')
    pepmodstates.df <- dplyr::select(evidence.df, !!!syms(pepmodstate_cols), pepmod_id_mq, protgroup_ids) %>%
      dplyr::distinct() %>% dplyr::arrange_at(pepmodstate_cols) %>%
      dplyr::mutate(pepmodstate_id = row_number())
    evidence.df <- dplyr::left_join(evidence.df,
                    dplyr::select(pepmodstates.df, pepmodstate_id, !!!syms(pepmodstate_cols)),
                    by=pepmodstate_cols)

    # summarize intensities for pepmodstate_id X msexperiment pair (there could be multiple ones)
    message('Reshaping, summarizing & weighting pepmodstate intensities...')
    pms_intensities_long.df <- tidyr::pivot_longer(evidence.df, starts_with("intensity"), names_to="mstag", values_to="intensity", names_prefix="intensity.") %>%
        dplyr::mutate(mstag = factor(mstag, levels = levels(mschannels.df$mstag)),
                      mass_error_w = pmax(replace_na(abs(mass_error_ppm), 1000), 0.01)^(-0.5)) %>%
        dplyr::inner_join(dplyr::select(mschannels.df, rawfile, mschannel, mstag), by=c("rawfile", "mstag")) %>%
        dplyr::group_by(pepmod_id, pepmodstate_id, msexperiment, msfraction, mstag, msrun, mschannel, rawfile) %>%
        dplyr::summarise(weight = weighted.mean(intensity, mass_error_w),
                         intensity = sum(intensity, na.rm=TRUE)) %>%
        dplyr::group_by(pepmod_id, mschannel) %>%
        dplyr::mutate(weight = pmax(weight/sum(weight), na_weight)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(intensity = if_else(!is.na(intensity) & intensity>0, intensity, NA_real_))

    message('Extracting pepmod information...')
    pepmodXmsrun_stats.df <- dplyr::mutate(evidence.df, has_quants = n_quants > 0) %>%
        dplyr::group_by(pepmod_id, msrun) %>%
        dplyr::summarize(n_charges = n_distinct(charge),
                         n_evidences = n_distinct(evidence_id),
                         n_quants = sum(has_quants),
                         n_full_quants = sum(is_full_quant)) %>%
        dplyr::ungroup()
    pepmod_stats.df <- dplyr::mutate(pepmodXmsrun_stats.df,
                                     has_quants = n_quants > 0,
                                     has_full_quants = n_full_quants > 0) %>%
        dplyr::group_by(pepmod_id) %>%
        dplyr::summarise(n_msruns = n(), # actually, rawfiles
                         n_quant_msruns = sum(has_quants),
                         n_full_quant_msruns = sum(has_quants),
                         n_max_charges=max(n_charges),
                         n_max_evidences=max(n_evidences)) %>%
        dplyr::ungroup()

    pepmods.df <- dplyr::select(evidence.df,
                                any_of(c("pepmod_id", "pepmod_id_mq", "peptide_id", "protgroup_ids",
                                       "protein_acs", "lead_protein_acs", "lead_razor_protein_ac",
                                       "gene_names", "protein_names", "is_reverse", "is_contaminant",
                                       "pepmod_seq", "peptide_seq", "peptide_len", "count_K", "count_R", "modifs", "n_miscleaved",
                                       "charge")),
                                tidyselect::ends_with(" site IDs")) %>%
        dplyr::distinct() %>%
        dplyr::group_by(pepmod_id) %>%
        dplyr::mutate(charges = paste0(sort(charge), collapse=' ')) %>%
        dplyr::filter(row_number() == 1L) %>% dplyr::ungroup() %>%
        dplyr::select(-charge) %>%
        dplyr::left_join(pepmod_stats.df) %>%
        dplyr::mutate(is_shared_by_groups = str_detect(protgroup_ids, stringr::fixed(';')),
                      is_shared = is_shared_by_groups,
                      is_shared_by_proteins = str_detect(protein_acs, stringr::fixed(';')))

    message( 'Extracting peaks information...' )
    peak_columns <- c('pepmodstate_id', 'pepmod_id', 'charge', 'msexperiment', 'msrun', 'msfraction', 'rawfile', 'evidence_id', 'n_quants', 'is_full_quant',
                      # Carbamidomethyl (C) Probabilities, Oxidation (M) Probabilities, Carbamidomethyl (C) Score Diffs, Oxidation (M) Score Diffs, Acetyl (Protein N-term), Carbamidomethyl (C), Oxidation (M),
                      'ident_type', 'label_state', 'ms2_mz', 'mz', 'mass_da', 'resolution',
                      'delta_mz_ppm', 'delta_mz_da', 'mass_error_ppm', 'Mass Error [mDa]', 'Mass Error [Da]',
                      'Uncalibrated Mass Error [ppm]', 'Uncalibrated Mass Error [mDa]', 'Uncalibrated Mass Error [Da]',
                      'Max intensity m/z 0', 'Max intensity m/z 1',
                      'Retention time', 'rt_len', 'rt_avg', 'rt_start', 'rt_finish',
                      'Retention time calibration', 'Match time difference', 'Match q-value', 'Match score',
                      'Number of data points', 'Number of scans', 'Number of isotopic peaks', 'PIF', 'Fraction of total spectrum',
                      'Base peak fraction', 'PEP', 'MS/MS Count', 'MS/MS Scan Number', 'Score', 'Delta score', 'Combinatorics',
                      'MS/MS IDs', 'Best MS/MS', 'AIF MS/MS IDs' ) %>%
      .[ . %in% colnames(evidence.df) ]
    peaks.df <- evidence.df %>% dplyr::select(!!peak_columns) %>% dplyr::distinct()

    ratio_columns.df <- tidyr::expand_grid(measure = 'ratio',
                                 mstag_nom = unique(dplyr::filter(intensity_columns.df, quant_type != 'aggregate')$mstag),
                                 mstag_denom = unique(dplyr::filter(intensity_columns.df, quant_type != 'aggregate')$mstag),
                                 type = c('', 'normalized', 'shift')) %>%
        dplyr::filter(("ratio" %in% import_data) & mstag_nom != mstag_denom) %>%
        dplyr::mutate(old_name = str_replace(paste0('Ratio ', mstag_nom, '/', mstag_denom, ' ', type), '\\s$', ''),
                      inverted_old_name = str_replace(paste0('Ratio ', mstag_denom, '/', mstag_nom, ' ', type), '\\s$', ''),
                      suffix = paste0(type, if_else(type != "", '_', ""), mstag_nom, mstag_denom),
                      new_name = paste0(measure, '.', suffix),
                      exists = old_name %in% colnames(evidence.df),
                      inverted_exists = inverted_old_name %in% colnames(evidence.df))

    message('Converting intensities to row format...')
    intensities.df <- pms_intensities_long.df %>%
        dplyr::mutate(mstag = factor(mstag, levels = as.character(intensity_columns.df$mstag)))
    ident_types.df <- dplyr::distinct(dplyr::select(peaks.df, pepmodstate_id, msexperiment, msrun, rawfile, ident_type)) %>%
        dplyr::mutate(ident_type_i = as_integer(ident_type)) %>%
        dplyr::group_by(pepmodstate_id, msexperiment, msrun, rawfile) %>%
        dplyr::summarise(ident_type_i = min(ident_type_i, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(ident_type = factor(levels(peaks.df$ident_type)[ident_type_i], levels=levels(peaks.df$ident_type)),
                      ident_type_i = NULL)
    intensities.df <- dplyr::left_join(intensities.df, ident_types.df)

    ratio_columns_sel.df <- dplyr::filter(ratio_columns.df,
                                          mstag_nom < mstag_denom & mstag_denom != 'Sum')
    if (any(ratio_columns_sel.df$exists)) {
        message('Extracting ratio data...')
        if (!all(ratio_columns_sel.df$exists)) {
            missing_cols <- dplyr::filter(ratio_columns_sel.df, !exists) %>% .$old_name
            warning(nrow(missing_cols), " ratio column(s) are missing: ", paste(missing_cols, collapse=' '))
        }
        sel_cols <- ratio_columns_sel.df$old_name
        names(sel_cols) <- ratio_columns_sel.df$new_name
        ratios.df <- intensities.df %>%
            dplyr::summarise_at(sel_cols, ~mean(.x, na.rm=TRUE)) %>%
            dplyr::ungroup()
        ratios.df[,dplyr::filter(ratio_columns_sel.df,!exists)$new_name] <- NA_real_ # add missing column names to make it evidence.df hapeable
        # fix NA
        ratios.df <- ratios.df %>% dplyr::mutate_at(ratio_columns_sel.df$new_name, zero2na)
        message('Converting ratios to long format...')
        ratio_cols <- str_c("ratio.", ratio_columns_sel.df$suffix)
        ratios.df <- tidyr::pivot_longer(select(ratios.df, !!id_cols, !!ratio_cols), cols=cols(ratio_cols),
                                            names_to='ratio_type', names_prefix="ratio.") %>%
            dplyr::inner_join(dplyr::select(ratio_columns_sel.df, suffix, type, mstag_nom, mstag_denom),
                                by = c("ratio_type" = "suffix"))
    } else {
        message('No channel ratios')
        ratios.df <- NULL
    }
    protgroup_ids <- unique(pepmods.df$protgroup_ids)
    protgroups2protgroup.list <- str_split(protgroup_ids, stringr::fixed(';'))
    protgroups2protgroup.df <- tibble(protgroup_ids = rep.int(protgroup_ids, sapply(protgroups2protgroup.list, length)),
                                      protgroup_id = as.integer(unlist(protgroups2protgroup.list)))
    res <- list(pepmods = pepmods.df,
                protgroups2protgroup = protgroups2protgroup.df,
                rawfiles = rawfiles.df,
                msexperiments = msexperiments.df,
                msruns = msruns.df,
                mschannels = mschannels.df,
                peaks = peaks.df,
                pepmodstate_intensities = intensities.df,
                pepmodstate_ratios = ratios.df,
                pepmodstates = pepmodstates.df)

    return (res)
}

#' @export
read.MaxQuant.Evidence <- function(folder_path, file_name = 'evidence.txt',
                                   nrows = Inf, guess_max = min(10000L, nrows),
                                   mschannel_annotate.f = NULL)
{
    process.MaxQuant.Evidence(read.MaxQuant.Evidence_internal(
        folder_path = folder_path, file_name = file_name, nrows = nrows, guess_max=guess_max),
        mschannel_annotate.f = mschannel_annotate.f,
        multidplyr_cluster = multidplyr_cluster)
}

msrun_code_parser.f = function(msruns, chunk_names = c('dataset', 'batch', 'frac_protocol', 'fraction', 'tech_replicate')) {
    msrun_chunks <- str_split(as.character(msruns), stringr::fixed('_'))
    n_missing_chunks <- length(chunk_names) - sapply(msrun_chunks, length)
    msrun_chunks <- lapply(seq_along(msrun_chunks), function(i) c(msrun_chunks[[i]], rep.int(NA, n_missing_chunks[[i]])))
    res <- cbind(msrun = msruns, as_tibble(do.call(rbind, msrun_chunks)))
    colnames(res) <- c('msrun', chunk_names)
    canonical_order <- c('msrun', 'dataset', 'experiment', 'batch', 'frac_protocol', 'replicate', 'fraction', 'tech_replicate')
    res[, setdiff(canonical_order, colnames(res))] <- NA
    res[, canonical_order]
}

#' @export
expand_protgroup_acs <- function(protgroups, acs_col,
                                 ac_col = str_replace(acs_col, "_(ac|id)s", "_\\1"),
                                 id_col="protgroup_id", sep=";") {
  acs <- str_split(protgroups[[acs_col]], sep, simplify = FALSE)
  res <- tibble(row_ix = rep.int(seq_along(acs), sapply(acs, length)),
                prot_ix = unlist(lapply(acs, function(acs) seq_along(acs))))
  res[[id_col]] <- protgroups[[id_col]][res$row_ix]
  res[[ac_col]] <- unlist(acs)
  return (res)
}

#' @export
match_protgroups_by_acs <- function(pgs1, pgs2, acs_col, suffix = c(".x", ".y")) {
    ac_col <- str_replace(acs_col, "_acs", "_ac")
    expd_pgs1 <- expand_protgroup_acs(pgs1, acs_col, ac_col)
    expd_pgs2 <- expand_protgroup_acs(pgs2, acs_col, ac_col)
    res <- dplyr::full_join(expd_pgs1, expd_pgs2, by = ac_col, suffix = suffix)
    res$is_matching <- !is.na(res[[paste0("protgroup_id", suffix[[1]])]]) &
        !is.na(res[[paste0("protgroup_id", suffix[[2]])]])
    return ( res )
}

#' @export
match_protgroups <- function(pgs1, pgs2, suffix = c(".x", ".y")) {
  pg2pg_by_majority_acs.df <- match_protgroups_by_acs(pgs1, pgs2, "majority_protein_acs") %>%
    dplyr::mutate(ac_match_rank = if_else(is_matching, 1L, NA_integer_))
  pg2pg_by_simple_acs.df <- match_protgroups_by_acs(pgs1, pgs2, "protein_acs") %>%
    dplyr::mutate(ac_match_rank = if_else(is_matching, 2L, NA_integer_))
  protgroup_vars <- c("protgroup_id.x", "protgroup_id.y")
  names(protgroup_vars) <- paste0("protgroup_id", suffix)
  pg2pg.df <- dplyr::bind_rows(pg2pg_by_majority_acs.df, pg2pg_by_simple_acs.df) %>%
    dplyr::group_by(protgroup_id.x, protgroup_id.y) %>%
    dplyr::summarize(min_ac_match_rank = min(ac_match_rank, na.rm = TRUE),
                     n_ac_matches = sum(!is.na(ac_match_rank) & ac_match_rank == min_ac_match_rank)) %>%
    # remove non-matches if there are matches
    dplyr::group_by(protgroup_id.x) %>%
    dplyr::filter(is.na(protgroup_id.x) | (is.na(min_ac_match_rank) == all(is.na(min_ac_match_rank)))) %>%
    dplyr::group_by(protgroup_id.y) %>%
    dplyr::filter(is.na(protgroup_id.y) | (is.na(min_ac_match_rank) == all(is.na(min_ac_match_rank)))) %>%
    dplyr::ungroup() %>%
    dplyr::rename(!!protgroup_vars)
}

#' @export
split_protgroups <- function(protgroups.df) {
    majority_acs.df <- expand_protgroup_acs(protgroups.df, acs_col = 'majority_protein_acs') %>%
        dplyr::mutate(is_majority = TRUE)
    expand_protgroup_acs(protgroups.df, acs_col = 'protein_acs') %>%
        mutate(is_contaminant = str_detect(protein_ac, "^CON_"),
               is_reverse = str_detect(protein_ac, "^REV_")) %>%
        dplyr::left_join(majority_acs.df, by=c("protgroup_id"="protgroup_id", "protein_ac" = "majority_protein_ac")) %>%
        dplyr::mutate(is_majority = if_else(is.na(is_majority), FALSE, TRUE))
}
