# Data helpers: wide <-> long conversion and parsing of pasted spreadsheet
# text. Pure functions (no Shiny) so they can be unit-tested.

META_COLS <- c("Sample", "Group", "Subgroup", "Ct_reference")

# Header spellings accepted for the optional second factor column.
SUBGROUP_ALIASES <- c("subgroup", "sub group", "sub-group", "factor2", "factor 2",
                      "factor_2", "treatment", "condition")
# Header spellings accepted for the reference (housekeeping) gene column.
REFERENCE_ALIASES <- c("ct_ref", "ct_reference", "reference", "reference gene",
                       "ref", "ref gene", "housekeeping", "hk", "hkg",
                       "ct reference", "ct ref")
# Cell contents that mean "no Ct value" (case-insensitive).
UNDETERMINED_TOKENS <- c("", "na", "n/a", "nan", "null", "undetermined", "undet",
                         "undeter", "no ct", "-", "--")

empty_wide <- function(with_subgroup = FALSE) {
  if (with_subgroup) {
    data.frame(Sample = character(0), Group = character(0),
               Subgroup = character(0), Ct_reference = numeric(0),
               stringsAsFactors = FALSE)
  } else {
    data.frame(Sample = character(0), Group = character(0),
               Ct_reference = numeric(0), stringsAsFactors = FALSE)
  }
}

empty_long <- function() {
  data.frame(Sample = character(0), Group = character(0),
             Gene = character(0), Ct_target = numeric(0),
             Ct_reference = numeric(0), Subgroup = character(0),
             stringsAsFactors = FALSE)
}

meta_cols <- function(wide) intersect(META_COLS, names(wide))

gene_cols <- function(wide) setdiff(names(wide), META_COLS)

# Convert a single spreadsheet cell to a Ct value. Accepts decimal commas
# ("25,3") and maps the usual "Undetermined"-style tokens to NA.
parse_ct <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  out <- rep(NA_real_, length(x))
  is_undet <- tolower(x) %in% UNDETERMINED_TOKENS
  y <- x[!is_undet]
  # "25,3" -> "25.3" but leave "1,234.5" alone (thousands separators are
  # not expected in Ct values, so a single comma is treated as decimal).
  y <- ifelse(grepl("^-?[0-9]+,[0-9]+$", y), sub(",", ".", y, fixed = TRUE), y)
  out[!is_undet] <- suppressWarnings(as.numeric(y))
  out
}

to_long_df <- function(wide) {
  gcols <- gene_cols(wide)
  if (length(gcols) == 0 || nrow(wide) == 0) return(empty_long())
  has_sub <- "Subgroup" %in% names(wide)
  n <- nrow(wide)
  blocks <- lapply(gcols, function(g) {
    r <- data.frame(
      Sample       = as.character(wide$Sample),
      Group        = as.character(wide$Group),
      Gene         = rep(g, n),
      Ct_target    = parse_ct(wide[[g]]),
      Ct_reference = parse_ct(wide$Ct_reference),
      stringsAsFactors = FALSE
    )
    if (has_sub) r$Subgroup <- as.character(wide$Subgroup)
    r$.row <- seq_len(n)
    r
  })
  out <- do.call(rbind, blocks)
  # Keep original row order (all genes of row 1, then row 2, ...), matching
  # the order the user pasted.
  out <- out[order(out$.row, match(out$Gene, gcols)), ]
  out$.row <- NULL
  rownames(out) <- NULL
  out
}

# Long -> wide. Rows sharing the same Sample/Group/Subgroup (technical
# replicates, or accidental duplicates) are kept as separate wide rows in
# their original order; each keeps its own Ct_reference.
to_wide_df <- function(long) {
  if (is.null(long) || nrow(long) == 0) return(empty_wide())
  genes <- unique(long$Gene[nchar(trimws(long$Gene)) > 0])
  if (length(genes) == 0) return(empty_wide())
  has_sub <- "Subgroup" %in% names(long)
  key_cols <- if (has_sub) c("Sample", "Group", "Subgroup") else c("Sample", "Group")
  key <- do.call(paste, c(long[key_cols], sep = "\t"))
  # Replicate index within (key, gene) so duplicates stay distinct rows.
  long$.rep <- ave(seq_len(nrow(long)), paste(key, long$Gene, sep = "\t"),
                   FUN = seq_along)
  first_idx <- !duplicated(paste(key, long$.rep, sep = "\t"))
  wide <- long[first_idx, c(key_cols, ".rep", "Ct_reference"), drop = FALSE]
  wide$.order <- seq_len(nrow(wide))
  for (g in genes) {
    gene_sub <- long[long$Gene == g, c(key_cols, ".rep", "Ct_target"), drop = FALSE]
    names(gene_sub)[ncol(gene_sub)] <- g
    wide <- merge(wide, gene_sub, by = c(key_cols, ".rep"), all.x = TRUE, sort = FALSE)
  }
  wide <- wide[order(wide$.order), , drop = FALSE]
  wide$.rep <- NULL
  wide$.order <- NULL
  rownames(wide) <- NULL
  wide
}

# Guess the field separator of pasted text: tab (Excel), else comma, else
# semicolon, else runs of 2+ spaces.
guess_separator <- function(header_line) {
  if (grepl("\t", header_line, fixed = TRUE)) return("\t")
  if (grepl(",", header_line, fixed = TRUE))  return(",")
  if (grepl(";", header_line, fixed = TRUE))  return(";")
  " {2,}"
}

split_fields <- function(line, sep) {
  f <- if (sep == " {2,}") strsplit(trimws(line), sep)[[1]]
       else strsplit(line, sep, fixed = TRUE)[[1]]
  trimws(f)
}

looks_numeric <- function(x) !is.na(parse_ct(x)) | tolower(trimws(x)) %in% UNDETERMINED_TOKENS[-1]

# Parse pasted spreadsheet text into a wide data frame.
#
# Expected layout: Sample | Group | [Subgroup] | Reference Ct | Gene1 | Gene2 ...
# The Subgroup column is recognised by its header (see SUBGROUP_ALIASES) in
# column 3 or 4. The reference column is recognised by header (REFERENCE_ALIASES)
# or, failing that, taken as the first column after the text columns.
#
# Returns list(ok = TRUE, wide = <data.frame>, warnings = <character>) or
# list(ok = FALSE, error = <message>).
parse_pasted_text <- function(txt) {
  fail <- function(msg) list(ok = FALSE, error = msg)
  if (is.null(txt) || !nzchar(trimws(txt))) return(fail("Nothing to paste."))
  txt   <- gsub("\r\n", "\n", txt, fixed = TRUE)
  txt   <- gsub("\r",   "\n", txt, fixed = TRUE)
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  lines <- lines[nchar(trimws(lines)) > 0]
  if (length(lines) < 2)
    return(fail("Paste needs a header row plus at least one data row."))

  sep     <- guess_separator(lines[1])
  headers <- split_fields(lines[1], sep)
  if (length(headers) < 4)
    return(fail(paste0("Header has only ", length(headers),
                       " column(s); expected at least 4: Sample, Group, Reference Ct, Gene1.")))

  # Header sanity: a header row should not be mostly numbers.
  if (sum(looks_numeric(headers[-(1:2)])) >= 2)
    return(fail("The first row looks like data, not a header. Include a header row (Sample, Group, Ct_ref, Gene1, ...)."))

  h_low <- tolower(trimws(headers))
  sub_idx <- which(h_low %in% SUBGROUP_ALIASES)
  sub_idx <- sub_idx[sub_idx %in% c(3, 4)]
  has_sub <- length(sub_idx) > 0
  sub_idx <- if (has_sub) sub_idx[1] else NA_integer_

  ref_idx <- which(h_low %in% REFERENCE_ALIASES)
  ref_idx <- ref_idx[ref_idx <= 5]
  if (length(ref_idx) == 0) ref_idx <- if (has_sub) 4 else 3 else ref_idx <- ref_idx[1]
  if (has_sub && ref_idx == sub_idx) ref_idx <- ref_idx + 1

  meta_idx  <- c(1, 2, if (has_sub) sub_idx, ref_idx)
  gene_idx  <- setdiff(seq_along(headers), meta_idx)
  gene_names <- headers[gene_idx]
  gene_names <- gene_names[nzchar(gene_names)]
  gene_idx   <- gene_idx[nzchar(headers[gene_idx])]
  if (length(gene_names) == 0)
    return(fail("No gene columns found after the reference column."))
  if (any(duplicated(gene_names)))
    return(fail(paste0("Duplicate gene column names: ",
                       paste(unique(gene_names[duplicated(gene_names)]), collapse = ", "))))

  warnings <- character(0)
  short_rows <- 0L
  rows <- lapply(lines[-1], function(line) {
    fields <- split_fields(line, sep)
    if (length(fields) < max(meta_idx)) { short_rows <<- short_rows + 1L; return(NULL) }
    fields <- c(fields, rep("", max(0, max(gene_idx) - length(fields))))
    r <- data.frame(Sample = fields[1], Group = fields[2], stringsAsFactors = FALSE)
    if (has_sub) r$Subgroup <- fields[sub_idx]
    r$Ct_reference <- parse_ct(fields[ref_idx])
    for (i in seq_along(gene_idx)) r[[gene_names[i]]] <- parse_ct(fields[gene_idx[i]])
    r
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(fail("No complete data rows found."))
  wide <- do.call(rbind, rows)
  rownames(wide) <- NULL

  if (short_rows > 0)
    warnings <- c(warnings, paste0(short_rows, " row(s) had too few columns and were skipped."))
  blank_meta <- !nzchar(wide$Sample) | !nzchar(wide$Group)
  if (any(blank_meta)) {
    warnings <- c(warnings, paste0(sum(blank_meta), " row(s) with a blank Sample or Group were skipped."))
    wide <- wide[!blank_meta, , drop = FALSE]
  }
  if (nrow(wide) == 0) return(fail("No usable data rows (Sample and Group must not be blank)."))
  n_na <- sum(is.na(wide$Ct_reference)) + sum(vapply(gene_names, function(g) sum(is.na(wide[[g]])), integer(1)))
  if (n_na > 0)
    warnings <- c(warnings, paste0(n_na, " Ct cell(s) were blank, 'Undetermined' or non-numeric and are treated as missing."))

  list(ok = TRUE, wide = wide, warnings = warnings, has_subgroup = has_sub)
}

# Human-readable list of missing Ct cells in a long data frame, for the
# validation panel. Returns character(0) when nothing is missing.
missing_ct_cells <- function(long, max_items = 6) {
  bad <- which(is.na(long$Ct_target) | is.na(long$Ct_reference))
  if (length(bad) == 0) return(character(0))
  lab <- paste0(long$Sample[bad], " / ", long$Gene[bad],
                ifelse(is.na(long$Ct_reference[bad]), " (reference)", ""))
  lab <- unique(lab)
  if (length(lab) > max_items) lab <- c(lab[seq_len(max_items)], paste0("... and ", length(lab) - max_items, " more"))
  lab
}
