options(stringsAsFactors = FALSE)

SOURCE_FILE <- Sys.getenv(
  "SOURCE_DATA_FILE",
  unset = file.path("data", "Additional_file_2_Source_data.xlsx")
)


# Source-data workbook sheet names are retained from the archived Additional file 2.
# Manuscript figure numbers were subsequently reordered when Methods was moved to
# the end of the main text; scripts therefore use semantic keys rather than
# hard-coding manuscript figure numbers.
SOURCE_SHEETS <- c(
  climate = "Fig1_climate",
  microclimate = "Fig2_microclimate",
  plant = "Plant_source",
  lmm = "LMM_results",
  ptdi = "Fig3_PTDI",
  module_trait = "Fig4_module_trait",
  functional = "Fig5_source",
  hub_expression = "Fig6_hub_expression"
)

source_sheet <- function(key) {
  if (!key %in% names(SOURCE_SHEETS)) {
    stop("Unknown source-data key: ", key, call. = FALSE)
  }
  unname(SOURCE_SHEETS[[key]])
}

assert_file <- function(path, label = path) {
  if (!file.exists(path)) {
    stop(sprintf("%s was not found at '%s'.", label, path), call. = FALSE)
  }
  invisible(path)
}

ensure_dirs <- function() {
  dir.create("results", showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path("results", "expression"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path("results", "intermediate"), showWarnings = FALSE, recursive = TRUE)
  dir.create("figures", showWarnings = FALSE, recursive = TRUE)
}

# Every simple source-data sheet contains two description rows and one blank row
# before the actual table header.
read_source_table <- function(sheet) {
  assert_file(SOURCE_FILE, "Additional file 2")
  readxl::read_excel(SOURCE_FILE, sheet = sheet, skip = 3)
}

read_count_matrix <- function(path) {
  assert_file(path, "Raw count matrix")

  if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    x <- readxl::read_excel(path)
  } else {
    x <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  }
  x <- as.data.frame(x, check.names = FALSE)
  if (ncol(x) < 2) stop("Count matrix must contain one gene-ID column and sample columns.")

  gene_id <- as.character(x[[1]])
  x[[1]] <- NULL
  mat <- as.matrix(x)
  storage.mode(mat) <- "integer"
  rownames(mat) <- gene_id
  mat
}

read_expression_matrix <- function(path) {
  assert_file(path, "Normalized expression matrix")

  if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    x <- readxl::read_excel(path)
  } else {
    x <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  }
  x <- as.data.frame(x, check.names = FALSE)
  if (ncol(x) < 2) stop("Expression matrix must contain one gene-ID column and sample columns.")

  gene_id <- as.character(x[[1]])
  x[[1]] <- NULL
  mat <- as.matrix(x)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- gene_id
  mat
}

prepare_plant_source <- function() {
  d <- read_source_table(source_sheet("plant"))
  d |>
    dplyr::mutate(
      Sample = as.character(Sample),
      Block = factor(Block),
      Warming = factor(Warming, levels = c("T0", "T1", "T2", "T3")),
      Watering = factor(Watering, levels = c("Ambient", "Added")),
      Species = factor(Species, levels = c("Aj", "Ce"))
    )
}

prepare_environment_source <- function() {
  d <- read_source_table(source_sheet("microclimate"))
  d |>
    dplyr::mutate(
      Sample = as.character(Sample),
      Block = factor(Block),
      Warming = factor(Warming, levels = c("T0", "T1", "T2", "T3")),
      Watering = factor(Watering, levels = c("Ambient", "Added"))
    )
}

species_metadata <- function(species, complete_traits = NULL) {
  d <- prepare_plant_source() |>
    dplyr::filter(Species == species)

  if (!is.null(complete_traits)) {
    d <- d |> tidyr::drop_na(dplyr::all_of(complete_traits))
  }

  d <- d |>
    dplyr::arrange(Watering, Warming, Block) |>
    dplyr::mutate(
      Group_Combo = factor(
        paste(Warming, Watering, sep = "_"),
        levels = c(
          "T0_Ambient", "T1_Ambient", "T2_Ambient", "T3_Ambient",
          "T0_Added", "T1_Added", "T2_Added", "T3_Added"
        )
      )
    )

  as.data.frame(d)
}

match_matrix_to_metadata <- function(mat, meta) {
  missing_samples <- setdiff(meta$Sample, colnames(mat))
  if (length(missing_samples) > 0) {
    stop(
      "Expression matrix is missing sample(s): ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }
  mat[, meta$Sample, drop = FALSE]
}

clean_lmer_anova <- function(tab, trait = NULL, analysis = NULL, species = NULL) {
  x <- as.data.frame(tab) |> tibble::rownames_to_column("Term")
  names(x) <- gsub("NumDF", "NumDF", names(x), fixed = TRUE)
  if ("F value" %in% names(x)) names(x)[names(x) == "F value"] <- "F_value"
  if ("Pr(>F)" %in% names(x)) names(x)[names(x) == "Pr(>F)"] <- "P_value"
  if ("DenDF" %in% names(x)) x$DenDF <- as.numeric(x$DenDF)

  if (!is.null(trait)) x <- dplyr::mutate(x, Trait = trait, .before = 1)
  if (!is.null(species)) x <- dplyr::mutate(x, Species = species, .before = 1)
  if (!is.null(analysis)) x <- dplyr::mutate(x, Analysis = analysis, .before = 1)
  x
}

extract_functional_source_tables <- function() {
  assert_file(SOURCE_FILE, "Additional file 2")
  raw <- readxl::read_excel(SOURCE_FILE, sheet = source_sheet("functional"), col_names = FALSE)
  raw <- as.data.frame(raw, check.names = FALSE)

  panel_b_row <- which(raw[[1]] == "Panel B. Functional gene subsets used for eigengene/PC1 analyses")
  go_header <- which(raw[[1]] == "Species" & raw[[2]] == "Module" & raw[[3]] == "GO_term")
  gene_header <- which(raw[[1]] == "Species" & raw[[2]] == "Module" & raw[[3]] == "Functional_subset")

  if (length(panel_b_row) != 1 || length(go_header) != 1 || length(gene_header) != 1) {
    stop("Could not locate both functional-source-data panels.")
  }

  go <- raw[(go_header + 1):(panel_b_row - 1), 1:6, drop = FALSE]
  go <- go[!is.na(go[[1]]) & go[[1]] != "", , drop = FALSE]
  colnames(go) <- c("Species", "Module", "GO_term", "GeneRatio", "Q_value", "Count")

  genes <- raw[(gene_header + 1):nrow(raw), 1:4, drop = FALSE]
  genes <- genes[!is.na(genes[[1]]) & genes[[1]] != "", , drop = FALSE]
  colnames(genes) <- c("Species", "Module", "Functional_subset", "Unigene_ID")

  list(go = go, genes = genes)
}

zscore_rows <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

mean_se_df <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x)
  se <- stats::sd(x) / sqrt(length(x))
  data.frame(y = m, ymin = m - se, ymax = m + se)
}
