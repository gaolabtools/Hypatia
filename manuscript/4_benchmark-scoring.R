# Score Hypatia and satuRn against simulation ground truth

suppressPackageStartupMessages({
  library(DEXSeq)
  library(precrec)
  library(dplyr)
})

# Precision-recall uses an absolute usage change of at least 0.10; null calibration uses zero change.
score_benchmark <- function(truth_file, rundiu_file, saturn_file, run_id) {
  effect_threshold <- 0.10
  safe_min <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else min(x)
  }

  simes_pvalue <- function(pvalue) {
    pvalue <- sort(pvalue[is.finite(pvalue)])
    n <- length(pvalue)
    if (n == 0L) {
      return(NA_real_)
    }
    min(1, min(n * pvalue / seq_len(n)))
  }

  safe_ratio <- function(numerator, denominator) {
    if (denominator == 0L) NA_real_ else numerator / denominator
  }

  message("Reading benchmark inputs.")
  truth_input <- readr::read_tsv(truth_file, show_col_types = FALSE)
  rundiu_input <- readRDS(rundiu_file)
  if (!"cramers.v" %in% names(rundiu_input$stats) &&
      "effect.size" %in% names(rundiu_input$stats)) {
    rundiu_input$stats <- dplyr::rename(
      rundiu_input$stats,
      cramers.v = effect.size
    )
  }
  saturn_input <- readRDS(saturn_file)
  if (!identical(saturn_input$run_id, run_id)) stop("satuRn result run_id does not match the benchmark.")

  filter_min_gene_pct <- 0.05
  filter_min_gene_cts <- 15
  filter_min_tx_cts <- 1
  filtered_genes <- rundiu_input$data %>%
    filter(cts.1 >= filter_min_tx_cts | cts.2 >= filter_min_tx_cts) %>%
    group_by(gene) %>%
    filter(n() > 1L) %>%
    summarize(
      gene.pct.1 = dplyr::first(gene.pct.1),
      gene.pct.2 = dplyr::first(gene.pct.2),
      gene.cts.1 = sum(cts.1),
      gene.cts.2 = sum(cts.2),
      .groups = "drop"
    ) %>%
    filter(
      gene.pct.1 >= filter_min_gene_pct,
      gene.pct.2 >= filter_min_gene_pct,
      gene.cts.1 >= filter_min_gene_cts,
      gene.cts.2 >= filter_min_gene_cts
    ) %>%
    pull(gene)

  required_truth_columns <- c("testID", "catID", "gt.delta")
  if (!all(required_truth_columns %in% names(truth_input))) {
    stop(
      "Ground truth is missing required columns: ",
      paste(setdiff(required_truth_columns, names(truth_input)), collapse = ", ")
    )
  }

  has_explicit_exact_null <- "gt.exact.null" %in% names(truth_input)
  if (has_explicit_exact_null) {
    truth_input$gt.exact.null <- as.logical(truth_input$gt.exact.null)
    if (anyNA(truth_input$gt.exact.null)) {
      stop("gt.exact.null must contain only TRUE/FALSE or 1/0 values.")
    }
  }

  truth_input <- truth_input %>%
    mutate(
      gene_id = paste0("test", testID),
      transcript_id = paste0("test", testID, ".cat", catID)
    )

  truth <- truth_input %>%
    summarize(
      max_abs_true_delta = max(abs(gt.delta)),
      meaningful_diu = any(abs(gt.delta) >= effect_threshold),
      exact_null = if (has_explicit_exact_null) {
        all(gt.exact.null)
      } else {
        all(gt.delta == 0)
      },
      exact_null_consistent = if (has_explicit_exact_null) {
        dplyr::n_distinct(gt.exact.null) == 1L
      } else {
        TRUE
      },
      .by = gene_id
    ) %>%
    mutate(
      statistical_diu = !exact_null,
      true_diu = meaningful_diu
    )

  if (anyDuplicated(truth$gene_id) || nrow(truth) == 0L) {
    stop("Unexpected ground-truth gene schema.")
  }
  if (!all(truth$exact_null_consistent)) {
    stop("gt.exact.null is not constant within one or more genes.")
  }
  if (any(truth$exact_null & truth$max_abs_true_delta != 0)) {
    stop("Exact-null genes must have zero true usage deltas.")
  }

  if (anyDuplicated(truth_input$transcript_id) || any(!is.finite(truth_input$gt.delta))) {
    stop("Truth requires unique transcript IDs and finite deltas.")
  }
  if (!setequal(truth_input$transcript_id, saturn_input$feature_filter$transcript_id) ||
      anyDuplicated(saturn_input$feature_filter$transcript_id) ||
      anyDuplicated(saturn_input$results$transcript_id) ||
      !all(rundiu_input$stats$gene %in% truth$gene_id)) {
    stop("Benchmark results and truth do not describe the same simulation.")
  }

  # Recover gene counts before satuRn filtering so untested genes remain in the evaluation.
  gene_counts <- saturn_input$feature_filter %>%
    summarize(total_count = sum(total_count), .by = gene_id)

  rundiu_gene <- rundiu_input$stats %>%
    transmute(
      gene_id = gene,
      rundiu_pvalue = pval,
      rundiu_fdr = padj,
      rundiu_permutation_pvalue = pval.perm,
      rundiu_permutation_fdr = padj.perm,
      rundiu_abs_delta = abs(max.prop.diff),
      rundiu_cramers_v = cramers.v
    )

  if (anyDuplicated(rundiu_gene$gene_id)) {
    stop("RunDIU results contain multiple rows per gene.")
  }
  if (!any(is.finite(rundiu_gene$rundiu_permutation_pvalue))) {
    stop(
      "RunDIU results do not contain finite permutation p-values. ",
      "Rerun RunDIU with permutation = TRUE."
    )
  }

  # Recalculate adjusted p-values within the filtered gene set.
  rundiu_filtered_fdr <- rundiu_gene %>%
    filter(gene_id %in% filtered_genes) %>%
    transmute(
      gene_id,
      rundiu_filtered_fdr = p.adjust(rundiu_pvalue, method = "BH"),
      rundiu_filtered_permutation_fdr = p.adjust(
        rundiu_permutation_pvalue,
        method = "BH"
      )
    )
  rundiu_gene <- rundiu_gene %>%
    left_join(rundiu_filtered_fdr, by = "gene_id")

  # Rank genes by their smallest empirical transcript p-value and use DEXSeq q-values for FDR calls.
  saturn_screen_input <- saturn_input$results %>%
    filter(is.finite(empirical_pval)) %>%
    dplyr::select(gene_id, empirical_pval)

  if (nrow(saturn_screen_input) == 0L) {
    stop("No finite satuRn empirical p-values are available for gene screening.")
  }
  if (any(
    saturn_screen_input$empirical_pval < 0 |
      saturn_screen_input$empirical_pval > 1
  )) {
    stop("satuRn empirical p-values must lie on [0, 1].")
  }

  calculate_saturn_gene_screen <- function(screen_input) {
    if (nrow(screen_input) == 0L) {
      return(data.frame(
        gene_id = character(),
        gene_screen_qvalue = numeric(),
        stringsAsFactors = FALSE
      ))
    }

    gene_factor <- factor(screen_input$gene_id)
    gene_split <- split(seq_len(nrow(screen_input)), gene_factor)
    screen_min_pvalue <- vapply(
      gene_split,
      function(i) min(screen_input$empirical_pval[i]),
      numeric(1)
    )
    screen_theta <- unique(sort(screen_min_pvalue))
    screen_qvalue_by_theta <- DEXSeq:::perGeneQValueExact(
      pGene = screen_min_pvalue,
      theta = screen_theta,
      geneSplit = gene_split
    )

    data.frame(
      gene_id = names(gene_split),
      gene_screen_qvalue = pmin(
        1,
        screen_qvalue_by_theta[match(screen_min_pvalue, screen_theta)]
      ),
      stringsAsFactors = FALSE
    )
  }

  saturn_gene_screen <- saturn_screen_input %>%
    calculate_saturn_gene_screen() %>%
    dplyr::rename(saturn_gene_screen_qvalue = gene_screen_qvalue)
  saturn_filtered_gene_screen <- saturn_screen_input %>%
    filter(gene_id %in% filtered_genes) %>%
    calculate_saturn_gene_screen() %>%
    dplyr::rename(
      saturn_filtered_gene_screen_qvalue = gene_screen_qvalue
    )

  # Keep raw and Simes-combined p-values for sensitivity analyses.
  saturn_gene <- saturn_input$results %>%
    summarize(
      saturn_min_pvalue = safe_min(pval),
      saturn_min_empirical_pvalue = safe_min(empirical_pval),
      saturn_min_empirical_fdr = safe_min(empirical_FDR),
      saturn_simes_empirical_pvalue = simes_pvalue(empirical_pval),
      saturn_tested_transcripts = sum(is.finite(empirical_pval)),
      .by = gene_id
    ) %>%
    mutate(
      saturn_simes_fdr = p.adjust(
        saturn_simes_empirical_pvalue,
        method = "BH"
      )
    ) %>%
    left_join(saturn_gene_screen, by = "gene_id")

  saturn_filtered_adjustments <- saturn_gene %>%
    filter(gene_id %in% filtered_genes) %>%
    transmute(
      gene_id,
      saturn_filtered_simes_fdr = p.adjust(
        saturn_simes_empirical_pvalue,
        method = "BH"
      )
    )
  saturn_filtered_transcript_fdr <- saturn_input$results %>%
    filter(gene_id %in% filtered_genes) %>%
    mutate(
      filtered_empirical_fdr = p.adjust(empirical_pval, method = "BH")
    ) %>%
    summarize(
      saturn_filtered_min_empirical_fdr = safe_min(
        filtered_empirical_fdr
      ),
      .by = gene_id
    )
  saturn_gene <- saturn_gene %>%
    left_join(saturn_filtered_gene_screen, by = "gene_id") %>%
    left_join(saturn_filtered_adjustments, by = "gene_id") %>%
    left_join(saturn_filtered_transcript_fdr, by = "gene_id")

  gene_scores <- truth %>%
    left_join(gene_counts, by = "gene_id") %>%
    left_join(rundiu_gene, by = "gene_id") %>%
    left_join(saturn_gene, by = "gene_id") %>%
    mutate(
      total_count = coalesce(total_count, 0),
      count_group = case_when(
        total_count < 100 ~ "<100",
        total_count < 200 ~ "100-200",
        TRUE ~ ">=200"
      ),
      rundiu_tested = is.finite(rundiu_pvalue),
      rundiu_permutation_tested = is.finite(rundiu_permutation_pvalue),
      saturn_tested = is.finite(saturn_min_empirical_pvalue),
      rundiu_pvalue_score = -log10(
        pmax(rundiu_pvalue, .Machine$double.xmin)
      ),
      rundiu_permutation_pvalue_score = -log10(
        pmax(rundiu_permutation_pvalue, .Machine$double.xmin)
      ),
      saturn_min_pvalue_score = -log10(
        pmax(saturn_min_pvalue, .Machine$double.xmin)
      ),
      saturn_min_empirical_pvalue_score = -log10(
        pmax(saturn_min_empirical_pvalue, .Machine$double.xmin)
      ),
      saturn_simes_empirical_pvalue_score = -log10(
        pmax(saturn_simes_empirical_pvalue, .Machine$double.xmin)
      )
    )

  # Give untested genes a score below every tested gene so they remain in the precision-recall analysis.
  ranking_methods <- data.frame(
    method = c("Hypatia p-value", "Hypatia empirical p-value",
               "satuRn minimum empirical p-value",
               "satuRn minimum raw p-value (sensitivity)",
               "satuRn Simes empirical p-value (sensitivity)"),
    pvalue = c("rundiu_pvalue", "rundiu_permutation_pvalue",
               "saturn_min_empirical_pvalue", "saturn_min_pvalue",
               "saturn_simes_empirical_pvalue"),
    score = c("rundiu_pvalue_score", "rundiu_permutation_pvalue_score",
              "saturn_min_empirical_pvalue_score", "saturn_min_pvalue_score",
              "saturn_simes_empirical_pvalue_score")
  )
  ranking_data <- bind_rows(lapply(seq_len(nrow(ranking_methods)), function(i) {
    definition <- ranking_methods[i, ]
    gene_scores %>% transmute(
      gene_id, true_diu, count_group, method = definition$method,
      tested = is.finite(.data[[definition$pvalue]]),
      score = coalesce(.data[[definition$score]], -1)
    )
  }))

  call_data <- bind_rows(
    gene_scores %>%
      transmute(
        gene_id, true_diu, count_group,
        method = "RunDIU BH FDR < 0.05",
        tested = rundiu_tested,
        called = coalesce(rundiu_fdr < 0.05, FALSE),
        called_filtered = coalesce(
          rundiu_filtered_fdr < 0.05,
          FALSE
        )
      ),
    gene_scores %>%
      transmute(
        gene_id, true_diu, count_group,
        method = "RunDIU BH FDR < 0.05 + Cramer's V >= 0.20",
        # Apply Cramer's V to the calls after BH adjustment, keeping all genes in the evaluation.
        tested = rundiu_tested,
        called = coalesce(
          rundiu_fdr < 0.05 & is.finite(rundiu_cramers_v) &
            rundiu_cramers_v >= 0.20,
          FALSE
        ),
        called_filtered = coalesce(
          rundiu_filtered_fdr < 0.05 & is.finite(rundiu_cramers_v) &
            rundiu_cramers_v >= 0.20,
          FALSE
        )
      ),
    gene_scores %>%
      transmute(
        gene_id, true_diu, count_group,
        method = "RunDIU permutation BH FDR < 0.05",
        tested = rundiu_permutation_tested,
        called = coalesce(rundiu_permutation_fdr < 0.05, FALSE),
        called_filtered = coalesce(
          rundiu_filtered_permutation_fdr < 0.05,
          FALSE
        )
      ),
    gene_scores %>%
      transmute(
        gene_id, true_diu, count_group,
        method = "satuRn gene-screen FDR < 0.05",
        tested = is.finite(saturn_gene_screen_qvalue),
        called = coalesce(saturn_gene_screen_qvalue < 0.05, FALSE),
        called_filtered = coalesce(
          saturn_filtered_gene_screen_qvalue < 0.05,
          FALSE
        )
      ),
    gene_scores %>%
      transmute(
        gene_id, true_diu, count_group,
        method = paste(
          "satuRn any transcript empirical FDR < 0.05",
          "(sensitivity)"
        ),
        tested = saturn_tested,
        called = coalesce(saturn_min_empirical_fdr < 0.05, FALSE),
        called_filtered = coalesce(
          saturn_filtered_min_empirical_fdr < 0.05,
          FALSE
        )
      ),
    gene_scores %>%
      transmute(
        gene_id, true_diu, count_group,
        method = "satuRn Simes BH FDR < 0.05 (sensitivity)",
        tested = saturn_tested,
        called = coalesce(saturn_simes_fdr < 0.05, FALSE),
        called_filtered = coalesce(
          saturn_filtered_simes_fdr < 0.05,
          FALSE
        )
      )
  )

  # Compare gene-level FDR calls using exact-null genes; keep Simes as a sensitivity analysis.
  calibration_call_data <- call_data %>%
    filter(method %in% c("RunDIU BH FDR < 0.05", "RunDIU permutation BH FDR < 0.05",
      "satuRn gene-screen FDR < 0.05", "satuRn Simes BH FDR < 0.05 (sensitivity)")) %>%
    left_join(dplyr::select(gene_scores, gene_id, meaningful_diu, statistical_diu, exact_null),
              by = "gene_id") %>%
    dplyr::select(gene_id, meaningful_diu, statistical_diu, exact_null, count_group,
                  method, tested, called, called_filtered)

  populations <- c(
    "All genes",
    "Filtered genes"
  )
  count_groups <- c("Overall", "<100", "100-200", ">=200")

  pr_curves <- list()
  auprc_summary <- list()
  result_index <- 1L

  for (population_name in populations) {
    for (count_group_name in count_groups) {
      for (method_name in unique(ranking_data$method)) {
        evaluation_data <- ranking_data %>%
          filter(method == method_name)

        if (population_name == "Filtered genes") {
          evaluation_data <- evaluation_data %>%
            filter(gene_id %in% filtered_genes)
        }
        if (count_group_name != "Overall") {
          evaluation_data <- evaluation_data %>%
            filter(count_group == count_group_name)
        }

        if (nrow(evaluation_data) == 0L ||
            length(unique(evaluation_data$true_diu)) < 2L) {
          next
        }

        evaluation <- precrec::evalmod(
          scores = evaluation_data$score,
          labels = as.integer(evaluation_data$true_diu)
        )
        pr <- evaluation$prcs[[1]]
        auprc <- attr(pr, "auc")

        pr_curves[[result_index]] <- data.frame(
          recall = pr$x,
          precision = pr$y,
          method = method_name,
          population = population_name,
          count_group = count_group_name,
          auprc = auprc,
          stringsAsFactors = FALSE
        )
        auprc_summary[[result_index]] <- data.frame(
          method = method_name,
          population = population_name,
          count_group = count_group_name,
          genes = nrow(evaluation_data),
          positives = sum(evaluation_data$true_diu),
          prevalence = mean(evaluation_data$true_diu),
          auprc = auprc,
          stringsAsFactors = FALSE
        )
        result_index <- result_index + 1L
      }
    }
  }

  pr_curves <- dplyr::bind_rows(c(list(data.frame(
    recall = numeric(), precision = numeric(), method = character(),
    population = character(), count_group = character(), auprc = numeric()
  )), pr_curves))
  auprc_summary <- dplyr::bind_rows(c(list(data.frame(
    method = character(), population = character(), count_group = character(),
    genes = integer(), positives = integer(), prevalence = numeric(), auprc = numeric()
  )), auprc_summary))

  binary_metrics <- list()
  result_index <- 1L

  for (population_name in populations) {
    for (count_group_name in count_groups) {
      for (method_name in unique(call_data$method)) {
        evaluation_data <- call_data %>% filter(method == method_name)

        if (population_name == "Filtered genes") {
          evaluation_data <- evaluation_data %>%
            filter(gene_id %in% filtered_genes) %>%
            mutate(called = called_filtered)
        }
        if (count_group_name != "Overall") {
          evaluation_data <- evaluation_data %>%
            filter(count_group == count_group_name)
        }
        if (nrow(evaluation_data) == 0L) {
          next
        }

        true_positive <- sum(evaluation_data$called & evaluation_data$true_diu)
        false_positive <- sum(evaluation_data$called & !evaluation_data$true_diu)
        true_negative <- sum(!evaluation_data$called & !evaluation_data$true_diu)
        false_negative <- sum(!evaluation_data$called & evaluation_data$true_diu)

        precision <- safe_ratio(true_positive, true_positive + false_positive)
        recall <- safe_ratio(true_positive, true_positive + false_negative)

        binary_metrics[[result_index]] <- data.frame(
          method = method_name,
          population = population_name,
          count_group = count_group_name,
          genes = nrow(evaluation_data),
          tested = sum(evaluation_data$tested),
          true_positive = true_positive,
          false_positive = false_positive,
          true_negative = true_negative,
          false_negative = false_negative,
          precision = precision,
          recall = recall,
          specificity = safe_ratio(
            true_negative,
            true_negative + false_positive
          ),
          f1 = if (
            is.na(precision) || is.na(recall) || precision + recall == 0
          ) {
            NA_real_
          } else {
            2 * precision * recall / (precision + recall)
          },
          stringsAsFactors = FALSE
        )
        result_index <- result_index + 1L
      }
    }
  }

  binary_metrics <- dplyr::bind_rows(binary_metrics)

  # Assess Hypatia and Simes at the gene level, and native satuRn tests at the transcript level.
  saturn_transcript_scored <- saturn_input$results %>%
    inner_join(
      truth_input %>%
        dplyr::select(transcript_id, gene_id, gt.delta),
      by = c("transcript_id", "gene_id")
    ) %>%
    left_join(
      gene_scores %>% dplyr::select(gene_id, exact_null, count_group),
      by = "gene_id"
    ) %>%
    mutate(
      transcript_exact_null = gt.delta == 0,
      transcript_statistical_diu = !transcript_exact_null,
      transcript_meaningful_diu = abs(gt.delta) >= effect_threshold
    )
  if (nrow(saturn_transcript_scored) != nrow(saturn_input$results)) {
    stop("Could not align every satuRn transcript result to simulation truth.")
  }

  null_pvalue_data <- bind_rows(
    gene_scores %>%
      transmute(
        gene_id, exact_null, count_group,
        method = "Hypatia p-value",
        hypothesis_level = "gene",
        pvalue = rundiu_pvalue
      ),
    gene_scores %>%
      transmute(
        gene_id, exact_null, count_group,
        method = "Hypatia empirical p-value",
        hypothesis_level = "gene",
        pvalue = rundiu_permutation_pvalue
      ),
    gene_scores %>%
      transmute(
        gene_id, exact_null, count_group,
        method = "satuRn Simes empirical p-value (sensitivity)",
        hypothesis_level = "gene",
        pvalue = saturn_simes_empirical_pvalue
      ),
    saturn_transcript_scored %>%
      transmute(
        gene_id, exact_null = transcript_exact_null, count_group,
        method = "satuRn raw p-value (transcript)",
        hypothesis_level = "transcript",
        pvalue = pval
      ),
    saturn_transcript_scored %>%
      transmute(
        gene_id, exact_null = transcript_exact_null, count_group,
        method = "satuRn empirical p-value (transcript)",
        hypothesis_level = "transcript",
        pvalue = empirical_pval
      )
  ) %>%
    filter(exact_null, is.finite(pvalue))

  alpha_levels <- c(0.01, 0.05, 0.10)
  if (nrow(null_pvalue_data) > 0L) {
    null_type1 <- list()
    result_index <- 1L
    for (count_group_name in count_groups) {
      for (method_name in unique(null_pvalue_data$method)) {
        evaluation_data <- null_pvalue_data %>% filter(method == method_name)
        if (count_group_name != "Overall") {
          evaluation_data <- evaluation_data %>%
            filter(count_group == count_group_name)
        }
        if (nrow(evaluation_data) == 0L) {
          next
        }
        for (alpha in alpha_levels) {
          null_type1[[result_index]] <- data.frame(
            method = method_name,
            hypothesis_level = unique(evaluation_data$hypothesis_level),
            count_group = count_group_name,
            alpha = alpha,
            null_hypotheses = nrow(evaluation_data),
            rejections = sum(evaluation_data$pvalue <= alpha),
            type_i_error = mean(evaluation_data$pvalue <= alpha),
            stringsAsFactors = FALSE
          )
          result_index <- result_index + 1L
        }
      }
    }
    null_type1 <- bind_rows(null_type1)
    null_qq_data <- null_pvalue_data %>%
      arrange(method, pvalue) %>%
      mutate(
        expected = (row_number() - 0.5) / n(),
        .by = c(method, hypothesis_level)
      )
  } else {
    null_type1 <- data.frame(
      method = character(), hypothesis_level = character(),
      count_group = character(), alpha = numeric(),
      null_hypotheses = integer(), rejections = integer(),
      type_i_error = numeric(), stringsAsFactors = FALSE
    )
    null_qq_data <- data.frame(
      gene_id = character(), exact_null = logical(), count_group = character(),
      method = character(), hypothesis_level = character(), pvalue = numeric(),
      expected = numeric(), stringsAsFactors = FALSE
    )
  }

  if (any(gene_scores$exact_null)) {
    statistical_metrics <- list()
    result_index <- 1L
    for (population_name in populations) {
      for (count_group_name in count_groups) {
        for (method_name in unique(calibration_call_data$method)) {
          evaluation_data <- calibration_call_data %>%
            filter(method == method_name, tested)
          if (population_name == "Filtered genes") {
            evaluation_data <- evaluation_data %>%
              filter(gene_id %in% filtered_genes) %>%
              mutate(called = called_filtered)
          }
          if (count_group_name != "Overall") {
            evaluation_data <- evaluation_data %>%
              filter(count_group == count_group_name)
          }
          if (nrow(evaluation_data) == 0L) {
            next
          }

          discoveries <- sum(evaluation_data$called)
          false_discoveries <- sum(
            evaluation_data$called & evaluation_data$exact_null
          )
          true_discoveries <- sum(
            evaluation_data$called & evaluation_data$statistical_diu
          )
          tested_null <- sum(evaluation_data$exact_null)
          tested_non_null <- sum(evaluation_data$statistical_diu)
          tested_meaningful <- sum(evaluation_data$meaningful_diu)
          tested_weak <- sum(
            evaluation_data$statistical_diu & !evaluation_data$meaningful_diu
          )

          statistical_metrics[[result_index]] <- data.frame(
            method = method_name,
            population = population_name,
            count_group = count_group_name,
            tested = nrow(evaluation_data),
            tested_null = tested_null,
            tested_non_null = tested_non_null,
            discoveries = discoveries,
            false_discoveries = false_discoveries,
            true_discoveries = true_discoveries,
            null_rejection_rate = safe_ratio(false_discoveries, tested_null),
            observed_fdr = safe_ratio(false_discoveries, discoveries),
            power_non_null = safe_ratio(true_discoveries, tested_non_null),
            power_meaningful = safe_ratio(
              sum(evaluation_data$called & evaluation_data$meaningful_diu),
              tested_meaningful
            ),
            power_weak = safe_ratio(
              sum(
                evaluation_data$called & evaluation_data$statistical_diu &
                  !evaluation_data$meaningful_diu
              ),
              tested_weak
            ),
            stringsAsFactors = FALSE
          )
          result_index <- result_index + 1L
        }
      }
    }
    statistical_metrics <- bind_rows(statistical_metrics)
  } else {
    statistical_metrics <- data.frame(
      method = character(), population = character(), count_group = character(),
      tested = integer(), tested_null = integer(), tested_non_null = integer(),
      discoveries = integer(), false_discoveries = integer(),
      true_discoveries = integer(), null_rejection_rate = numeric(),
      observed_fdr = numeric(), power_non_null = numeric(),
      power_meaningful = numeric(), power_weak = numeric(),
      stringsAsFactors = FALSE
    )
  }

  primary_gene_level_fdr_methods <- c(
    "RunDIU BH FDR < 0.05",
    "satuRn gene-screen FDR < 0.05"
  )
  gene_level_fdr_metrics <- statistical_metrics %>%
    filter(method %in% primary_gene_level_fdr_methods)
  gene_level_effect_metrics <- binary_metrics %>%
    filter(method %in% primary_gene_level_fdr_methods)

  if (any(saturn_transcript_scored$transcript_exact_null)) {
    saturn_transcript_metrics <- list()
    result_index <- 1L
    for (count_group_name in count_groups) {
      evaluation_data <- saturn_transcript_scored %>%
        filter(is.finite(empirical_FDR))
      if (count_group_name != "Overall") {
        evaluation_data <- evaluation_data %>%
          filter(count_group == count_group_name)
      }
      if (nrow(evaluation_data) == 0L) {
        next
      }

      called <- evaluation_data$empirical_FDR < 0.05
      discoveries <- sum(called)
      false_discoveries <- sum(
        called & evaluation_data$transcript_exact_null
      )
      true_discoveries <- sum(
        called & evaluation_data$transcript_statistical_diu
      )
      tested_null <- sum(evaluation_data$transcript_exact_null)
      tested_non_null <- sum(evaluation_data$transcript_statistical_diu)
      tested_meaningful <- sum(evaluation_data$transcript_meaningful_diu)

      saturn_transcript_metrics[[result_index]] <- data.frame(
        method = "satuRn empirical FDR < 0.05",
        hypothesis_level = "transcript",
        count_group = count_group_name,
        tested = nrow(evaluation_data),
        tested_null = tested_null,
        tested_non_null = tested_non_null,
        discoveries = discoveries,
        false_discoveries = false_discoveries,
        true_discoveries = true_discoveries,
        null_rejection_rate = safe_ratio(false_discoveries, tested_null),
        observed_fdr = safe_ratio(false_discoveries, discoveries),
        power_non_null = safe_ratio(true_discoveries, tested_non_null),
        power_meaningful = safe_ratio(
          sum(called & evaluation_data$transcript_meaningful_diu),
          tested_meaningful
        ),
        stringsAsFactors = FALSE
      )
      result_index <- result_index + 1L
    }
    saturn_transcript_metrics <- bind_rows(saturn_transcript_metrics)
  } else {
    saturn_transcript_metrics <- data.frame(
      method = character(), hypothesis_level = character(),
      count_group = character(), tested = integer(), tested_null = integer(),
      tested_non_null = integer(), discoveries = integer(),
      false_discoveries = integer(), true_discoveries = integer(),
      null_rejection_rate = numeric(), observed_fdr = numeric(),
      power_non_null = numeric(), power_meaningful = numeric(),
      stringsAsFactors = FALSE
    )
  }

  population_summary <- bind_rows(
    gene_scores %>%
      summarize(
        population = "All genes",
        genes = n(),
        meaningful_diu = sum(meaningful_diu),
        statistical_diu = sum(statistical_diu),
        exact_null = sum(exact_null)
      ),
    gene_scores %>%
      filter(gene_id %in% filtered_genes) %>%
      summarize(
        population = "Filtered genes",
        genes = n(),
        meaningful_diu = sum(meaningful_diu),
        statistical_diu = sum(statistical_diu),
        exact_null = sum(exact_null)
      )
  )

  scoring_result <- list(
    run_id = run_id,
    truth_definition = paste0(
      "any transcript with absolute true delta >= ", effect_threshold
    ),
    truth_definitions = list(
      meaningful_diu = paste0(
        "any transcript with absolute true delta >= ", effect_threshold
      ),
      exact_null = "all transcript true deltas equal zero",
      statistical_diu = "not an exact-null gene"
    ),
    has_explicit_exact_null = has_explicit_exact_null,
    input_files = list(
      truth = truth_file,
      RunDIU = rundiu_file,
      satuRn = saturn_file
    ),
    calibration_notes = c(
      "PRC uses the practical meaningful-DIU effect definition.",
      paste(
        "RunDIU BH and satuRn DEXSeq gene-screen fixed-FDR metrics are",
        "gene-level."
      ),
      paste(
        "Adjusted p-values and gene-screen q-values are recomputed within the",
        "filtered-gene testing family."
      ),
      paste(
        "satuRn native null calibration is transcript-level; raw and Simes",
        "results are stored only as sensitivity analyses."
      ),
      paste(
        "DEXSeq gene-screen q-values are excluded from null calibration because",
        "adjusted q-values are not expected to be uniform under the null."
      )
    ),
    saturn_gene_screening = list(
      method = "DEXSeq:::perGeneQValueExact",
      dexseq_version = as.character(utils::packageVersion("DEXSeq")),
      source_pvalue = "satuRn empirical_pval"
    ),
    population_summary = population_summary,
    gene_scores = gene_scores,
    ranking_data = ranking_data,
    call_data = call_data,
    auprc = auprc_summary,
    pr_curves = pr_curves,
    effect_metrics = binary_metrics,
    binary_metrics = binary_metrics,
    null_pvalues = null_pvalue_data,
    null_qq = null_qq_data,
    null_type1 = null_type1,
    calibration_call_data = calibration_call_data,
    gene_level_fdr_metrics = gene_level_fdr_metrics,
    gene_level_effect_metrics = gene_level_effect_metrics,
    statistical_metrics = statistical_metrics,
    saturn_transcript_metrics = saturn_transcript_metrics,
    session_info = utils::sessionInfo()
  )

  scoring_result
}

run_id <- "run_3"
data_dir <- "manuscript/figure1/data"
artifact <- function(suffix) file.path(data_dir, paste0(run_id, suffix))
truth_file <- artifact("_results.tsv")
hypatia_file <- artifact("_hypatia_res.rds")
saturn_file <- artifact("_saturn_res.rds")
input_files <- c(truth_file, hypatia_file, saturn_file)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) stop("Missing benchmark input(s): ", paste(missing_files, collapse = ", "))

scores <- score_benchmark(truth_file, hypatia_file, saturn_file, run_id)
saveRDS(scores, artifact("_benchmark_scores.rds"), compress = "gzip")
for (name in c("population_summary", "auprc", "binary_metrics", "null_type1",
               "statistical_metrics", "saturn_transcript_metrics")) {
  readr::write_tsv(scores[[name]], artifact(paste0("_", name, ".tsv")))
}
saveRDS(utils::sessionInfo(), artifact("_benchmark_scoring_session.rds"))
print(scores$auprc[scores$auprc$count_group == "Overall", ])
