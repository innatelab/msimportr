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

gsub_columns <- function(df, from, to) {
  colnames(df) <- str_replace(colnames(df), from, to)
  df
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
            gsub_columns("^Intensity\\s([LMH](\\s|$))", "Intensity.\\1") %>%
            gsub_columns("^Intensity(\\s|$)", "Intensity.Sum\\1") %>%
            dplyr::mutate_all(zero2na)
        res.df <- dplyr::bind_cols(res.df, intensities.df)
        col_info$intensity <- colnames(intensities.df)
    }
    if ('LFQ' %in% import_data) {
        lfq.df <- proteinGroups.df %>% dplyr::select(starts_with("LFQ intensity")) %>%
            gsub_columns("^LFQ intensity\\s([LMH](\\s|$))", "LFQ_Intensity.\\1") %>%
            gsub_columns("^LFQ intensity(\\s|$)", "LFQ_Intensity.Sum\\1") %>%
            dplyr::mutate_all(zero2na)
        res.df <- dplyr::bind_cols(res.df, lfq.df)
        col_info$LFQ <- colnames(lfq.df)
    }
    if ('iBAQ' %in% import_data) {
        ibaq.df <- proteinGroups.df %>% dplyr::select(starts_with("iBAQ")) %>%
            gsub_columns("^iBAQ\\s([LMH](\\s|$))", "iBAQ.\\1") %>%
            gsub_columns("^iBAQ(\\s|$)", "iBAQ.Sum\\1") %>%
            dplyr::mutate_all(zero2na)
        res.df <- dplyr::bind_cols(res.df, ibaq.df)
        col_info$iBAQ <- colnames(ibaq.df)
    }
    if ('ratio' %in% import_data) {
        ratios.df <- proteinGroups.df %>% dplyr::select(starts_with("Ratio")) %>%
            gsub_columns("^Ratio\\s", "Ratio.")
        res.df <- dplyr::bind_cols(res.df, ratios.df)
        col_info$ratio <- colnames(ratios.df)
    }
    if ('ident_type' %in% import_data) {
        ident_types.df <- proteinGroups.df %>% dplyr::select(starts_with("Identification type")) %>%
          gsub_columns("^Identification\\stype\\s", "ident_type.") %>%
          dplyr::mutate_all(~factor(.x, levels=c("By matching", "By MS/MS")))
        res.df <- dplyr::bind_cols(res.df, ident_types.df)
        col_info$ident_type <- colnames(ident_types.df)
    }
    if ('ms2_count' %in% import_data) {
      ms2_counts.df <- proteinGroups.df %>% dplyr::select(starts_with("MS/MS count ")) %>%
        gsub_columns("^MS/MS\\scount\\s", "ms2_count.")
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
    col_renames <- c(seq = "Sequence", seq_len = "Length",
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
            gsub_columns("^Intensity\\s([LMH](\\s|$))", "Intensity.\\1") %>%
            gsub_columns("^Intensity(\\s|$)", "Intensity.Sum\\1") %>%
            dplyr::mutate_all(zero2na)
        res.df <- dplyr::bind_cols(res.df, intensities.df)
        col_info$intensity <- colnames(intensities.df)
    }
    if ('aa_stats' %in% import_data) {
       aa_stats.df <- dplyr::select(peptides.df, matches("^[A-Z] Count$")) %>%
        gsub_columns("([A-Z]) Count", "count_\\1")
       res.df <- bind_cols(res.df, aa_stats.df)
       col_info$aa_stats <- colnames(aa_stats.df)
    }
    if ('ident_type' %in% import_data) {
        ident_types.df <- dplyr::select(peptides.df, starts_with("Identification type")) %>%
          gsub_columns("^Identification\\stype\\s", "ident_type.") %>%
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
                                    raw_file = `Raw file`,
                                    protein_acs = `Proteins`,
                                    dp_protein_acs = `DP Proteins`,
                                    #protgroup_ids = `Protein group IDs`,
                                    #fasta_headers = `Fasta headers`,
                                    dp_modif = `DP Modification`,
                                    dp_mass_diff = `DP Mass Difference`) %>%
        dplyr::mutate(pep_type = factor(pep_type),
                      raw_file = factor(raw_file),
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
            gsub_columns("^Intensity\\s([LMH](\\s|_|$))", "Intensity.\\1") %>%
            gsub_columns("^Intensity(\\s|_|$)", "Intensity.Sum\\1") %>%
            gsub_columns("^Intensity.Sum(.+)(___\\d+)$", "Intensity.Sum\\2\\1") %>%
            dplyr::mutate_all(~zero2na(as.numeric(.x)))
        res.df <- dplyr::bind_cols(res.df, intensities.df)
        col_info$intensity <- colnames(intensities.df)
    }
    if ('occupancy' %in% import_data) {
        occupancies.df <- data.df %>% dplyr::select(starts_with("Occupancy")) %>%
            gsub_columns("^(Occupancy\\s|Occupancy ratio|Occupancy error scale\\s)([LMH](\\s|$))", "\\1.\\2") %>%
            gsub_columns("^Occupancy ratio", "occupancy_ratio") %>%
            gsub_columns("^Occupancy error scale", "occupancy_error_scale") %>%
            dplyr::mutate_all(as.numeric)
        res.df <- dplyr::bind_cols(res.df, occupancies.df)
        col_info$occupancy <- colnames(occupancies.df)
    }
    if ('ratio' %in% import_data) {
        ratios.df <- data.df %>% dplyr::select(starts_with("Ratio mode/base")) %>%
            gsub_columns("^Ratio mod/base\\s", "ratio_mod2base.") %>%
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
                            col_types = readr::cols(Resolution = 'n',
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
                     raw_file = "Raw file", msrun = "Experiment", ident_type = "Type", label_state = 'Labeling State',
                     mass_error_ppm = "Mass Error [ppm]", mass_error_da = "Mass Error [Da]",
                     mass_error_ppm = "Mass error [ppm]", mass_error_da = "Mass error [Da]",
                     uncalib_mass_error_ppm = "Uncalibrated Mass Error [ppm]",
                     uncalib_mass_error_da = "Uncalibrated Mass Error [Da]",
                     seq = 'Sequence', seq_len = 'Length', count_K = 'K Count', count_R = 'R Count', modifs = 'Modifications',
                     n_miscleaved = "Missed cleavages", n_miscleaved = "Missed cleavages (Trypsin/P;LysC/P)", n_miscleaved_lysc = "Missed cleavages (LysC/P)",
                     mod_seq = 'Modified sequence', charge = 'Charge',
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
    if (any(str_detect(colnames(evidence.df), "Intensity "))) {
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
                                 msrun = factor(msrun),
                                 raw_file = factor(raw_file),
                                 ident_type = factor(ident_type, levels=MaxQuant_IdentTypes),
                                 is_contaminant = replace_na(is_contaminant, "") == '+',
                                 is_reverse = replace_na(is_reverse, "") == '+')
    return ( evidence.df )
}

glm_corrected_intensities <- function(intensities.df,
                                      use_mstags=TRUE, pepmodXmsrun_iaction = FALSE,
                                      glm_max_factor = 50.0, max_npepmods = 50L,
                                      glm_reldelta_range = c(-1+1E-5, 9.0),
                                      min_intensity = 1E+3) {
    npepmods <- n_distinct(intensities.df$pepmod_id)
    if (max_npepmods > 0 && npepmods > max_npepmods) {
        pepmod_stats.df <- dplyr::group_by(intensities.df, pepmod_id) %>%
            dplyr::summarise(n_observed = sum(observed),
                             n_charges = n_distinct(charge),
                             max_intensity = max(intensity, na.rm=TRUE)) %>%
            dplyr::ungroup() %>% dplyr::arrange(desc(n_observed), desc(n_charges), desc(max_intensity)) %>%
            dplyr::filter(row_number() <= max_npepmods)
        intensities.df <- dplyr::semi_join(intensities.df, pepmod_stats.df)
    }
    intensities.df <- dplyr::group_by(intensities.df, msrun) %>%
        dplyr::mutate(glm_used = !all(is.na(intensity))) %>% dplyr::ungroup()
    glm_mask <- intensities.df$glm_used
    glm_intensities.df <- intensities.df[glm_mask,] %>%
        dplyr::mutate(intensity_fixed = if_else(!is.na(intensity) & intensity > min_intensity,
                                                intensity, min_intensity)+
                                         0.05*min_intensity*runif(n(), -1, 1))
    ndatapoints <- nrow(glm_intensities.df)
    nmsprotocols <- if ('msprotocol' %in% colnames(glm_intensities.df)) n_distinct(glm_intensities.df$msprotocol) else 1L
    msrunXcharge <- dplyr::distinct(dplyr::select(dplyr::filter(glm_intensities.df, observed), msrun, charge))
    # GLM could be run also with ncharges = n_distinct(charge) > 1, but often generates artifacts
    ncharges <- max(table(msrunXcharge$msrun))
    nmsruns <- n_distinct(glm_intensities.df$msrun)
    msrunXtag <- dplyr::distinct(dplyr::select(dplyr::filter(glm_intensities.df, observed), msrun, mstag))
    nmstags <- max(table(msrunXtag$msrun))
    rq.method <- ifelse(nmsruns*npepmods > 100L, "lasso", "lasso")
    rhs_str <- "0"
    if (nmsruns >= 1) {
        if (nmstags > 1) {
            rhs_str <- paste0(rhs_str, " + msrun + msrun:mstag")
        } else {
            rhs_str <- paste0(rhs_str, " + msrun")
        }
        #if (nmsprotocols > 1L) {
        #  rhs_str <- str_replace_all(rhs_str, " + (msrun|mschannel)", " + msprotocol + msprotocol:\\1")
        #}
    }
    if (npepmods > 1) {
        rhs_str <- paste0(rhs_str, " + pepmod_id")
        if (nmsruns > 1) {
            # model both individual pepmod regulation and LC variation
            # (the latter affects all charges of the pepmod in the same way)
            rhs_str <- paste0(rhs_str, " + msrun:pepmod_id")
            if (nmstags > 1 && pepmodXmsrun_iaction) {
                # model both individual pepmod regulation in mstags
                rhs_str <- paste0(rhs_str, " + msrun:mstag:pepmod_id")
            }
        }
    }
    if (ncharges > 1) {
        if (npepmods > 1) {
            #glm_intensities.df <- dplyr::mutate(glm_intensities.df,
            #                                    pepmod_charge = factor(paste0(pepmod_id, '_', charge), ordered = FALSE))
            # FIXME should be pepmod_id + pepmod_id: but results in singular matrix
            rhs_str <- paste0(rhs_str, " + pepmod_id:charge")
        } else {
            rhs_str <- paste0(rhs_str, " + charge")
        }
    }
    max_intensity_fixed = max(glm_intensities.df$intensity_fixed, na.rm = TRUE)
    intensities.df$intensity_glm <- NA_real_
    intensities.df$glm_rhs <- rhs_str
    intensities.df$glm_method <- rq.method
    intensities.df$glm_ndup_effects <- 0L
    if (ndatapoints > 1L && rhs_str != "0 + msrun" && rhs_str != "0 + mschannel") {
        intensities.df$glm_status <- "failed"

        # convert GLM terms into factors to generate proper model and
        # use only the factors present in the given intensities chunk
        glm_intensities.df <- dplyr::mutate(glm_intensities.df,
            intensity_fixed_norm = intensity_fixed/max_intensity_fixed,
            pepmod_id = factor(pepmod_id, ordered = FALSE),
            charge = factor(charge, ordered = FALSE),
            mstag = factor(mstag, ordered = FALSE),
            msprotocol = factor(msprotocol, ordered = FALSE)
        )
        message("Pepmod_id=", glm_intensities.df$pepmod_id[1], " protgroup_ids=", glm_intensities.df$protgroup_ids[1],
                ' glm_rhs=', rhs_str, ' rq.method=', rq.method)
        fla <- as.formula(paste0("log(intensity_fixed_norm) ~ ", rhs_str))
        mod_mtx <- model.matrix(fla, data=glm_intensities.df)
        glm_intensities.df$glm_ndup_effects <- sum(duplicated(mod_mtx, MARGIN=2L))
        mod_mtx_attrs <- attributes(mod_mtx)
        mod_mtx <- mod_mtx[,!duplicated(mod_mtx, MARGIN=2L) & colSums(mod_mtx) > 0]
        mod_mtx_attrs[names(attributes(mod_mtx))] <- attributes(mod_mtx)
        attributes(mod_mtx) <- mod_mtx_attrs
        #for (ltl in -6:1) { # try different scales
        if (requireNamespace("quantreg")) {
        tryCatch({
            # FIXME GLM assumes one pepmod_id
            #glm_res <- glm(fla,
            #               data = intensities.df,
            #               weights = intensities.df$weight)
            glm_res <- rq.wfit(mod_mtx, log(glm_intensities.df$intensity_fixed_norm),
                               weights = glm_intensities.df$weight, method=rq.method)
            class(glm_res) <- "rq"
            #glm_res <- rq(fla, data=glm_intensities.df, weights = glm_intensities.df$weight, method = rq.method)
            #     control = lmRob.control(tl=exp(ltl), tlo=1E-4, tua=1.5E-4, mxs=100, mxr=100, mxf=100,
            #                            initial.alg = "random"))#, final.alg = "adaptive"))
            #glm_res <- lmRob(fla,
            #                data = glm_intensities.df,
            #                weights = glm_intensities.df$weight,
            #                control = control.params)
            #glm_res <- rq(log(intensity.Sum) ~ msrun + pepmod_id,
            #              data = glm_intensities.df %>% dplyr::mutate(pepmod_id = factor(pepmod_id)),
            #              weights = glm_intensities.df$weight)
            intensities.df$intensity_glm <- NA_real_
            intensities.df$intensity_fixed <- NA_real_
            intensities.df[glm_mask, 'intensity_glm'] <- exp(predict(glm_res))*max_intensity_fixed
            intensities.df[glm_mask, 'intensity_fixed'] <- glm_intensities.df$intensity_fixed
            intensities.df$glm_status <- "success"
            #    break
        }, error = function(e) warning(e), finally = invisible())
        }
        #}
    } else {
        intensities.df$glm_status <- "skipped"
    }
    if (intensities.df$glm_status[1] == "success") {
        # GLM succeeded, now check if its results are valid and correct
        intensities.df <- dplyr::mutate(intensities.df,
                                        intensity_glm_reldelta = intensity_glm/intensity_fixed - 1,
                                        is_valid = if_else(observed, dplyr::between(intensity_glm_reldelta, glm_reldelta_range[[1]], glm_reldelta_range[[2]]), TRUE))
        if (all(intensities.df$is_valid)) {
            # use GLM predictions correcting for
            # anomalously big GLM predictions (e.g. when charge of higest intensity was missing and got predicted)
            intensities.df <- dplyr::mutate(intensities.df,
                                            intensity_corr = pmin(intensity_glm, glm_max_factor*max_intensity_fixed) )
        } else {
            # GLM generates outliers that do not make sense, use the original intensities
            intensities.df$intensity_corr <-intensities.df$intensity
            intensities.df$glm_status <- "reverted"
        }
        intensities.df <- dplyr::mutate(intensities.df,
                                        intensity_corr_reldelta = (intensity_corr - intensity_fixed) / pmax(intensity_corr, intensity_fixed) )
    } else {
        # GLM failed, use the original intensities
        intensities.df <- dplyr::mutate(intensities.df,
                                        intensity_glm_reldelta = NA_real_,
                                        is_valid = NA,
                                        intensity_corr = intensity,
                                        intensity_corr_reldelta = if_else(!is.na(intensity), 0.0, NA_real_))
    }
    return (intensities.df)
}

process.MaxQuant.Evidence <- function( evidence.df, evidence.pepobj = c("pepmod", "pepmodstate"),
                                       evidence.msobj = c("msrun", "mschannel"),
                                       min_pepmodstate_freq = 0.9, min_essential_freq = 0.0,
                                       import_data = c("intensity"),
                                       correct_ratios = TRUE,
                                       correct_by_ratio.ref_label = NA,
                                       mode = c("labeled", "label-free"),
                                       mschannel_annotate.f = NULL,
                                       na_weight = 1E-5, min_intensity = 1E+3,
                                       glm.min_intensity = 50*min_intensity,
                                       glm.max_factor = 50, glm.max_npepmods=25L,
                                       glm.context = c("protgroup", "pepmod"),
                                       glm.reldelta_range = c(-1+1E-5, 99.0),
                                       multidplyr_cluster = NULL )
{
    quant_mode <- match.arg(mode)
    evidence_pepobj <- match.arg(evidence.pepobj)
    evidence_msobj <- match.arg(evidence.msobj)
    glm_context <- match.arg(glm.context)
    complete_names <- function( vars ) {
        res <- vars
        res_names <- names(vars)
        res_names[res_names==""] <- vars[res_names==""]
        names(res) <- res_names
        return (res)
    }
    message("Extracting MS runs and MS channels info...")
    ilabels <- factor(c("H", "M", "L", "F", 'Sum'),
                      levels = c("H", "M", "L", "F", 'Sum'))
    intensity_columns.df <- tidyr::expand_grid(measure = 'intensity', mstag = ilabels) %>%
        mutate(old_name = paste0('Intensity ', mstag),
               new_name = paste0(measure, '.', mstag),
               type = case_when(mstag == 'Sum' ~ 'aggregate',
                                old_name %in% colnames(evidence.df) ~ 'measured',
                                TRUE ~ 'missing')) %>%
        dplyr::filter(type != 'missing') %>%
        dplyr::mutate(type = factor(type)) %>%
        dplyr::arrange(mstag)
    ilabels <- intensity_columns.df$mstag # restrict to the labels actually used
    msruns.df <- dplyr::select(evidence.df, msrun, raw_file) %>% dplyr::distinct()
    mschannels.df <- tidyr::expand_grid(raw_file = msruns.df$raw_file,
                                        mstag = ilabels) %>%
        dplyr::inner_join(msruns.df) %>%
        dplyr::mutate(
            mschannel = interaction(msrun, mstag, drop=TRUE, lex.order=TRUE, sep='_'),
            quant_type = if_else(mstag %in% c('H','M','L'), 'SILAC',
                                 if_else(mstag == 'Sum', 'aggregate', 'label_free')))
    if (any(mschannels.df$quant_type == 'SILAC')) {
        message('Evidence table contains SILAC-labeled data')
        mode <- 'labeled'
    } else {
        message('Evidence table contains label-free data')
        model <- 'label-free'
    }
    if (!is.null(mschannel_annotate.f)) {
        # get user-defined mschannel annotations
        annot.df <- mschannel_annotate.f(mschannels.df)
        if (nrow(annot.df) != nrow(mschannels.df)) {
            stop("mschannel annotations has incorrect number of rows")
        }
        mschannels_annot.df <- dplyr::left_join(mschannels.df, annot.df)
        if (nrow(mschannels_annot.df) != nrow(mschannels.df)) {
            stop("mschannel annotations do not correctly match mschannels")
        }
        mschannels.df <- mschannels_annot.df
        if ('is_msrun_used' %in% colnames(mschannels.df)) {
            message('Removing unused user-specified msruns')
            mschannels.df <- dplyr::filter(mschannels.df, !is.na(is_msrun_used) & is_msrun_used) %>%
                # drop unused levels
                dplyr::mutate(msrun = factor(msrun),
                              raw_file = factor(raw_file),
                              mschannel = factor(mschannel))
            evidence.df <- dplyr::filter(evidence.df, as.character(raw_file) %in% mschannels.df$raw_file) %>%
                dplyr::mutate(raw_file = factor(raw_file, levels=levels(mschannels.df$raw_file)),
                              msrun = factor(msrun, levels=levels(mschannels.df$msrun)))
            if (nrow(evidence.df) == 0) {
                stop("Empty evidence frame after is_msrun_used filtering")
            }
        }
        if (n_distinct(mschannels.df$mschannel) != nrow(mschannels.df)) {
            stop('mschannel ids are not unique')
        }
    }
    mschannels.df <- dplyr::mutate(mschannels.df, quant_type = factor(quant_type))
    # NOTE?: the same mod. peptide Id can have multiple mod. sequences (mod at different poses)
    message('Enumerating pepmod states and summarizing channel intensities...')
    intensities.mtx <- as.matrix(evidence.df[,dplyr::filter(intensity_columns.df, mstag != 'Sum')$old_name])
    intens_cols <- intensity_columns.df$old_name
    names(intens_cols) <- intensity_columns.df$new_name
    evidence.df <- dplyr::mutate(evidence.df,
                          pepmodstate_id = as.integer(interaction(pepmod_id, charge, drop = TRUE, lex.order = TRUE, sep = '_')),
                          `Intensity Sum` = Matrix::rowSums(intensities.mtx, na.rm = TRUE),
                          n_quants = Matrix::rowSums(!is.na(intensities.mtx) & intensities.mtx > 0.0),
                          is_full_quant = !is.na(Matrix::rowSums(intensities.mtx > 0.0))) %>%
        dplyr::rename(!!intens_cols)

    message('Extracting pepmod states...')
    pepmodstates.df <- dplyr::select(evidence.df, protgroup_ids, pepmodstate_id, pepmod_id, charge) %>% dplyr::distinct() %>%
      dplyr::arrange(pepmodstate_id)

    # summarize intensities for pepmodstate_id X msrun pair (there could be multiple ones)
    message('Reshaping, summarizing & weighting pepmodstate intensities...')
    pms_intensities_long.df <- tidyr::pivot_longer(evidence.df, starts_with("intensity"), names_to="mstag", values_to="intensity", names_prefix="intensity.") %>%
        dplyr::mutate(mstag = factor(mstag, levels = levels(mschannels.df$mstag))) %>%
        dplyr::inner_join(dplyr::select(mschannels.df, msrun, raw_file, mstag, mschannel)) %>%
        dplyr::group_by(pepmod_id, pepmodstate_id, msrun, raw_file, mstag, mschannel) %>%
        dplyr::summarise(weight = weighted.mean(intensity, abs(1/mass_error_ppm)^0.5),
                         intensity = sum(intensity, na.rm=TRUE)) %>%
        dplyr::group_by(pepmod_id, mschannel) %>%
        dplyr::mutate(weight = pmax(weight/sum(weight), na_weight)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(intensity = if_else(!is.na(intensity) & intensity>0, intensity, NA_real_))

    message('Expanding intensities data.frame for every pepmod state and mschannel...')
    pms_full_intensities_long.df <- tidyr::expand(pms_intensities_long.df,
                                                  nesting(raw_file, mschannel), pepmodstate_id) %>%
        dplyr::left_join(dplyr::select(mschannels.df, mschannel, msrun, raw_file, any_of(c("mstag", "msprotocol")))) %>%
        dplyr::left_join(pepmodstates.df) %>%
        dplyr::left_join(pms_intensities_long.df) %>%
        dplyr::mutate(observed = !is.na(intensity) & intensity > 0,
                      weight = if_else(observed & is.finite(weight), weight, na_weight))
    if (!("msprotocol" %in% colnames(pms_full_intensities_long.df))) {
        pms_full_intensities_long.df$msprotocol <- "default"
    }

    message('Correcting & predicting missing intensities...')
    glm_group_col <- switch(glm_context,
                            pepmod = 'pepmod_id',
                            protgroup = 'protgroup_ids',
                            stop('Unsupported glm.context ', glm_context))
    pms_full_intensities_long_glm.df <- dplyr::filter(pms_full_intensities_long.df, mstag != 'Sum')
    pms_full_intensities_long_glm.df <- if (!correct_ratios) {
      dplyr::mutate(pms_full_intensities_long_glm.df,
                    intensity_fixed = intensity,
                    intensity_corr = NA_real_,
                    intensity_glm = NA_real_,
                    glm_rhs = NA_character_,
                    glm_status = "skipped",
                    glm_method = NA_character_,
                    glm_ndup_effects = 0L)
    } else if (!is.null(multidplyr_cluster)) {
        message('Distributing GLM tasks over the cluster')
        parallel::clusterEvalQ(cl=multidplyr_cluster, require(quantreg))
        parallel::clusterExport(cl=multidplyr_cluster, "glm_corrected_intensities")
        parallel::clusterExport(cl=multidplyr_cluster, c("glm.max_factor", "glm.reldelta_range", "glm.max_npepmods", "min_intensity"),
                                envir=environment())
        distr.df <- switch(glm_context,
                           pepmod = multidplyr::partition(pms_full_intensities_long_glm.df, pepmod_id, msprotocol, cluster=multidplyr_cluster),
                           protgroup = multidplyr::partition(pms_full_intensities_long_glm.df, protgroup_ids, msprotocol, cluster=multidplyr_cluster)) %>%
        #multidplyr::partition_(pms_full_intensities_long.df, c(glm_group_col, "msprotocol"), cluster=multidplyr_cluster) %>%
        dplyr::do({glm_corrected_intensities(., glm_max_factor = glm.max_factor,
                                             glm_reldelta_range = glm.reldelta_range,
                                             max_npepmods = glm.max_npepmods,
                                             min_intensity = min_intensity)})
        message("Collecting GLM results")
        dplyr::collect(distr.df)
    } else {
        dplyr::group_by(pms_full_intensities_long_glm.df, !!glm_group_col, msprotocol) %>%
            dplyr::do({glm_corrected_intensities(., glm_max_factor = glm.max_factor, glm_reldelta_range = glm.reldelta_range,
                                                 max_npepmods = glm.max_npepmods, min_intensity=min_intensity)})
    }
    pms_full_intensities_long_glm.df <- dplyr::ungroup(pms_full_intensities_long_glm.df) %>%
        dplyr::select(-intensity_fixed) %>% # remove temporary non-NA column
        dplyr::bind_rows(dplyr::filter(pms_full_intensities_long.df, mstag == 'Sum'))

    if (evidence_pepobj == "pepmod") {
      message('Summing intensities of different charges...')
      full_intensities_long.df <- pms_full_intensities_long_glm.df %>%
          dplyr::group_by(pepmod_id, msrun, raw_file, mschannel, mstag) %>%
          dplyr::summarise(intensity = sum(intensity, na.rm=TRUE),
                           intensity_glm = sum(intensity_glm, na.rm=TRUE),
                           intensity_corr = sum(intensity_corr, na.rm=TRUE))
    } else {
      full_intensities_long.df <- pms_full_intensities_long_glm.df
    }
    full_intensities_long.df <- dplyr::ungroup(full_intensities_long.df) %>%
          dplyr::mutate(intensity = if_else(intensity > 0, intensity, NA_real_),
                        intensity_glm = if_else(mstag == 'Sum', NA_real_, intensity_glm),
                        intensity_corr = if_else(intensity_corr > glm.min_intensity, intensity_corr, NA_real_),
                        intensity_glm_reldelta = intensity_glm/intensity - 1.0,
                        intensity_corr_reldelta = intensity_corr/intensity - 1.0)

    message('Extracting pepmod information...')
    pepmodXmsrun_stats.df <- evidence.df %>% dplyr::group_by(protgroup_ids, pepmod_id, msrun, raw_file) %>%
        dplyr::summarize(n_charges = n_distinct(charge),
                         n_evidences = n_distinct(evidence_id),
                         n_quants = sum(n_quants > 0),
                         n_full_quants = sum(is_full_quant)) %>%
        dplyr::ungroup()
    prediction.df <- dplyr::group_by_at(pms_full_intensities_long_glm.df, glm_group_col) %>% # FIXME support multiple protocols
        # FIXME for different msprotocols GLM could be different
        dplyr::summarise(glm_rhs = glm_rhs[1], glm_status = glm_status[1],
                         glm_method = glm_method[1], glm_ndup_effects = glm_ndup_effects[1]) %>% dplyr::ungroup() %>%
        dplyr::mutate(glm_rhs = factor(glm_rhs),
                      glm_status = factor(glm_status, levels=c("success", "skipped", "reverted", "failed")),
                      glm_method = factor(glm_method))
    pepmod_stats.df <- dplyr::group_by(pepmodXmsrun_stats.df, protgroup_ids, pepmod_id) %>%
        dplyr::summarise(n_msruns = n(), # actually, raw_files
                         n_quant_msruns = sum(n_quants > 0),
                         n_full_quant_msruns = sum(n_full_quants > 0),
                         n_max_charges=max(n_charges),
                         n_max_evidences=max(n_evidences)) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(prediction.df)

    pepmods.df <- dplyr::select(evidence.df,
                                any_of(c("pepmod_id", "peptide_id", "protgroup_ids",
                                       "protein_acs", "lead_protein_acs", "lead_razor_protein_ac",
                                       "gene_names", "protein_names", "is_reverse", "is_contaminant",
                                       "mod_seq", "seq", "seq_len", "count_K", "count_R", "modifs", "n_miscleaved",
                                       # Carbamidomethyl (C) site IDs', Oxidation (M) site IDs
                                       "charge"))) %>%
        dplyr::distinct() %>%
        dplyr::group_by(pepmod_id) %>%
        dplyr::mutate(charges = paste0(sort(charge), collapse=' ')) %>%
        dplyr::filter(row_number() == 1L) %>%
        dplyr::ungroup() %>%
        dplyr::select(-charge) %>%
        dplyr::left_join(pepmod_stats.df) %>%
        dplyr::mutate(is_shared_by_groups = str_detect(protgroup_ids, stringr::fixed(';')),
                      is_shared = is_shared_by_groups,
                      is_shared_by_proteins = str_detect(protein_acs, stringr::fixed(';')))

    message( 'Extracting peaks information...' )
    peak_columns <- c('pepmodstate_id', 'pepmod_id', 'charge', 'msrun', 'raw_file', 'evidence_id', 'n_quants', 'is_full_quant',
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
                                 mstag_nom = unique(dplyr::filter(intensity_columns.df, type == 'measured')$mstag),
                                 mstag_denom = unique(dplyr::filter(intensity_columns.df, type == 'measured')$mstag),
                                 type = c('', 'normalized', 'shift')) %>%
        dplyr::filter(("ratio" %in% import_data) & mstag_nom != mstag_denom) %>%
        dplyr::mutate(old_name = str_replace(paste0('Ratio ', mstag_nom, '/', mstag_denom, ' ', type), '\\s$', ''),
                      inverted_old_name = str_replace(paste0('Ratio ', mstag_denom, '/', mstag_nom, ' ', type), '\\s$', ''),
                      suffix = paste0(type, if_else(type != "", '_', ""), mstag_nom, mstag_denom),
                      new_name = paste0(measure, '.', suffix),
                      exists = old_name %in% colnames(evidence.df),
                      inverted_exists = inverted_old_name %in% colnames(evidence.df))

    # quant_mode defines which features are used for quantitation
    # label-free mode favours features that are observed in multiple runs
    # labeled mode favours features that have most quantitations per label

    if ("ratio" %in% import_data) {
        if (correct_ratios) {
            # channel/msrun ratios are more precisely measured than intensities even within the features of same pepmod
            # so calculate average ratios per pepmod and use these ratios to correct the absolute intensities:
            # the ratios of the corrected intensities should match
            message( "Correcting intensities by averaging ratios..." )
            if (quant_mode == "label-free") {
                if (!("condition" %in% colnames(mschannels.df))) {
                    stop("No condition annotations, cannot average ratios, specify mschannel_annotate.f function")
                }
                rlm.control <- lmrob.control("KS2014")
                intensities.df <- dplyr::left_join(intensities.df, mschannels.df) %>%
                    dplyr::left_join(dplyr::select(pepmods.df, pepmod_id, protgroup_ids)) %>%
                    dplyr::group_by(condition, protgroup_ids)
                intensities.df <- dplyr::do(intensities.df, {glm_corrected_intensities(., rlm.control)})
                intensities.df <- intensities.df %>%
                    dplyr::ungroup() %>% dplyr::select(-protgroup_ids, -condition)
            } else if (quant_mode == "labeled") {
                ref_ratio_cols.df <- dplyr::bind_rows(
                    dplyr::filter(ratio_columns.df, exists & mstag_denom == correct_by_ratio.ref_label) %>%
                    dplyr::mutate(is_inverted = FALSE,
                                  trf_old_name = old_name),
                    dplyr::filter(ratio_columns.df, exists & mstag_nom == correct_by_ratio.ref_label) %>%
                    dplyr::mutate(mstag_nom = mstag_denom,
                                  mstag_denom = correct_by_ratio.ref_label,
                                  is_inverted = TRUE,
                                  trf_old_name = gsub('\\s$', '', paste0(measure, ' ', mstag_nom, '/', mstag_denom, ' ', type)))
                    ) %>%
                    dplyr::filter(type == "")
                inverted_cols <- ref_ratio_cols.df$old_name[ref_ratio_cols.df$is_inverted]
                names(inverted_cols) <- ref_ratio_cols.df$trf_old_name[ref_ratio_cols.df$is_inverted]
                pre_intensities.df$ref_intensity <- pre_intensities.df[[paste0("Intensity ", correct_by_ratio.ref_label)]]

                agg_ratios.df <- pre_intensities.df %>% dplyr::select(pepmod_id, msrun, raw_file, ref_intensity, starts_with("Ratio"), mass_error_ppm) %>%
                    dplyr::mutate(ratio_weight = abs(1/mass_error_ppm)/sum(abs(1/mass_error_ppm), na.rm=TRUE)) %>%
                    dplyr::mutate_at(inverted_cols, ~ 1/.) %>%
                    dplyr::summarise_at(ref_ratio_cols.df$trf_old_name,
                                        ~ if_else(all(is.na(ratio_weight)), NA_real_,
                                                 weighted.mean(.[!is.na(ratio_weight)],
                                                               ratio_weight[!is.na(ratio_weight)], na.rm=TRUE)))
                ilabels_ordered <- c(as.character(ref_ratio_cols.df$mstag_nom), correct_by_ratio.ref_label)
                intens_mtx <- as.matrix(intensities.df[paste0("intensity.", ilabels_ordered)])
                ratios_mtx <- cbind(as.matrix(agg_ratios.df[ref_ratio_cols.df$trf_old_name]),
                                    if_else(is.na(intens_mtx[[ncol(intens_mtx)]]), 0, 1))
                w_avg_intens <- Matrix::rowSums(intens_mtx * ratios_mtx, na.rm = TRUE)
                ratio_sqr_sum <- Matrix::rowSums(ratios_mtx * ratios_mtx, na.rm = TRUE)
                # corrected reference intensity that would minimize
                # the total sqr deviation of corrected intensities vs original intensities
                ref_intensity_corr <- w_avg_intens/ratio_sqr_sum
                intens_corr_mtx <- ratios_mtx * replicate(ncol(ratios_mtx), ref_intensity_corr)
                intensities.df[paste0("intensity_corr.", ilabels_ordered)] <- intens_corr_mtx
                intensities.df$intensity_corr.Sum <- Matrix::rowSums(intens_corr_mtx, na.rm = TRUE)
                # restore intensities that are NA because no ratio available
                for (lbl in ilabels) {
                    col_corr <- paste0("intensity_corr.", lbl)
                    col_orig <- paste0("intensity.", lbl)
                    intensities.df[col_corr] <- if_else(is.na(intensities.df[[col_corr]]) | intensities.df[[col_corr]]==0.0,
                                                        intensities.df[[col_orig]],
                                                        intensities.df[[col_corr]])
                }
            }
        }
    }
    message('Converting intensities to ', evidence_msobj, ' row format...')
    if (evidence_msobj == 'mschannel') {
        #print(str(intensities.df))
        intensities.df <- full_intensities_long.df %>% dplyr::filter(observed) %>%
            dplyr::mutate(mstag = factor(mstag, levels = as.character(intensity_columns.df$mstag)))
    } else if (evidence_msobj == 'msrun') {
        intensities.df <- full_intensities_long.df %>% ungroup() %>% select(-mschannel) %>%
          tidyr::pivot_wider(id_cols=c("pepmod_id", "pepmodstate_id", "msrun", "raw_file"),
                             names_from=mstag,
                             values_from=starts_with("intensity"), names_sep=".")
    }
    ident_types.df <- dplyr::distinct(dplyr::select(peaks.df, pepmodstate_id, msrun, raw_file, ident_type)) %>%
        dplyr::mutate(ident_type_i = as_integer(ident_type)) %>%
        dplyr::group_by(pepmodstate_id, msrun, raw_file) %>%
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
        if (evidence_msobj == 'mschannel') {
            message('Converting ratios to long format...')
            ratio_cols <- str_c("ratio.", ratio_columns_sel.df$suffix)
            ratios.df <- tidyr::pivot_longer(select(ratios.df, !!id_cols, !!ratio_cols), cols=cols(ratio_cols),
                                             names_to='ratio_type', names_prefix="ratio.") %>%
                dplyr::inner_join(dplyr::select(ratio_columns_sel.df, suffix, type, mstag_nom, mstag_denom),
                                  by = c("ratio_type" = "suffix"))
        }
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
                raw_files = msruns.df,
                mschannels = mschannels.df,
                peaks = peaks.df)
    res[[paste0(evidence_pepobj, "_intensities")]] <- intensities.df
    res[[paste0(evidence_pepobj, "_ratios")]] <- ratios.df
    if (evidence_pepobj == "pepmodstate") {
        res$pepmodstates <- pepmodstates.df
    }
    return (res)
}

#' @export
read.MaxQuant.Evidence <- function(folder_path, file_name = 'evidence.txt',
                                   evidence.pepobj = c("pepmod", "pepmodstate"),
                                   evidence.msobj = c("msrun", "mschannel"),
                                   nrows = Inf, guess_max = min(10000L, nrows),
                                   min_pepmodstate_freq = 0.9, min_essential_freq = 0.0,
                                   correct_ratios = TRUE, correct_by_ratio.ref_label = NA,
                                   mode = c("labeled", "label-free"), mschannel_annotate.f = NULL,
                                   glm.context = c("protgroup", "pepmod"), glm.max_npepmods=50L,
                                   multidplyr_cluster = NULL)
{
    process.MaxQuant.Evidence(read.MaxQuant.Evidence_internal(
        folder_path = folder_path, file_name = file_name, nrows = nrows, guess_max=guess_max),
        evidence.pepobj = evidence.pepobj, evidence.msobj = evidence.msobj,
        min_pepmodstate_freq = min_pepmodstate_freq,
        min_essential_freq = min_essential_freq,
        correct_ratios = correct_ratios,
        correct_by_ratio.ref_label = correct_by_ratio.ref_label,
        mode = mode, glm.context = glm.context, glm.max_npepmods = glm.max_npepmods,
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
