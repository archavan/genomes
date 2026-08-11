# Arun Chavan and Claude
# Started: 2026-02-22

# Utility functions for constructing informative HOG display labels
# Source this file to use the functions in other scripts.
#
# Main entry point: format_hog_label()
#
# Label format: m5b1_Nr1i(2,3).Thr(a,b).Vdr
#   - m5b1 = 5 mouse genes, 1 bflo gene
#   - clusters separated by .
#   - each cluster: Stem(suffixes) or bare Stem
#   - consecutive numeric suffixes collapsed to ranges (1-3)
#   - Gm\d+ and \d+Rik names collapsed to Gm(xN) / Rik(xN)

library(stringr)

# Longest common prefix of two strings
lcp2 <- function(a, b) {
  a_chars <- utf8ToInt(a)
  b_chars <- utf8ToInt(b)
  n <- min(length(a_chars), length(b_chars))
  if (n == 0) {
    return("")
  }
  match <- a_chars[1:n] == b_chars[1:n]
  prefix_len <- if (all(match)) n else which(!match)[1] - 1
  substr(a, 1, prefix_len)
}

# Longest common prefix of a character vector
lcp <- function(x) {
  if (length(x) == 0) {
    return("")
  }
  if (length(x) == 1) {
    return(x)
  }
  Reduce(lcp2, x)
}

# Cluster a vector of symbols by longest common prefix (LCP).
# Returns a list of lists, each with $stem and $suffixes.
# Symbols sharing a prefix of >= min_prefix chars are grouped.
# Ungrouped symbols become singletons with suffix = "".
cluster_symbols <- function(syms, min_prefix = 3L) {
  syms <- sort(unique(syms))
  if (length(syms) == 0) {
    return(list())
  }
  if (length(syms) == 1) {
    return(list(list(stem = syms, suffixes = "")))
  }

  # Pre-filter uninformative names
  is_gm <- str_detect(syms, "^Gm\\d+$")
  is_rik <- str_detect(syms, "Rik\\d*$")
  informative <- syms[!is_gm & !is_rik]
  gm_syms <- syms[is_gm]
  rik_syms <- syms[is_rik]

  clusters <- list()

  # Cluster informative symbols by LCP
  remaining <- sort(informative)
  while (length(remaining) > 0) {
    # Take the first remaining symbol, find all others sharing LCP >= min_prefix
    anchor <- remaining[1]
    prefixes <- vapply(remaining[-1], function(s) lcp2(anchor, s), character(1))
    shared <- c(TRUE, nchar(prefixes) >= min_prefix)

    if (sum(shared) == 1) {
      # No partner found: singleton
      clusters <- c(clusters, list(list(stem = anchor, suffixes = "")))
      remaining <- remaining[-1]
    } else {
      # Found a group — refine the stem to the LCP of all group members
      group <- remaining[shared]
      stem <- lcp(group)

      # Now expand: check if any remaining non-group symbols also share this stem
      not_in_group <- remaining[!shared]
      also_match <- nchar(stem) >= min_prefix &
        vapply(
          not_in_group,
          function(s) {
            nchar(lcp2(stem, s)) >= nchar(stem)
          },
          logical(1)
        )
      if (any(also_match)) {
        group <- c(group, not_in_group[also_match])
        group <- sort(group)
        stem <- lcp(group)
        not_in_group <- not_in_group[!also_match]
      }

      suffixes <- str_remove(group, fixed(stem))
      # If all suffixes are empty (identical symbols?), keep as-is
      clusters <- c(clusters, list(list(stem = stem, suffixes = suffixes)))
      remaining <- not_in_group
    }
  }

  # Add Gm and Rik summary clusters
  if (length(gm_syms) > 0) {
    clusters <- c(
      clusters,
      list(list(
        stem = "Gm",
        suffixes = paste0("x", length(gm_syms))
      ))
    )
  }
  if (length(rik_syms) > 0) {
    clusters <- c(
      clusters,
      list(list(
        stem = "Rik",
        suffixes = paste0("x", length(rik_syms))
      ))
    )
  }

  clusters
}

# Natural sort key: pad numeric portions so lexicographic sort = numeric sort.
# "1a" -> "00001a", "10" -> "00010", "3a" -> "00003a"
natural_sort_key <- function(x) {
  str_replace_all(x, "\\d+", function(m) sprintf("%05d", as.integer(m)))
}

# Collapse a sorted vector of integers into range notation: 1,2,3,5,7,8,9 -> "1-3,5,7-9"
collapse_runs <- function(x) {
  if (length(x) == 0) {
    return("")
  }
  if (length(x) == 1) {
    return(as.character(x))
  }

  breaks <- c(1, which(diff(x) != 1) + 1)
  ends <- c(breaks[-1] - 1, length(x))

  parts <- vapply(
    seq_along(breaks),
    function(i) {
      start <- x[breaks[i]]
      end <- x[ends[i]]
      if (start == end) {
        as.character(start)
      } else if (end - start == 1) {
        paste0(start, ",", end)
      } else {
        paste0(start, "-", end)
      }
    },
    character(1)
  )

  paste(parts, collapse = ",")
}

# Format a vector of suffixes into a compact string.
# Purely numeric suffixes: collapse consecutive runs into ranges (1-3).
# Otherwise: comma-separate.
format_suffixes <- function(suffixes) {
  if (length(suffixes) == 1 && suffixes == "") {
    return(NULL)
  }

  # Check for xN summary (Gm/Rik)
  if (length(suffixes) == 1 && str_detect(suffixes, "^x\\d")) {
    return(suffixes)
  }

  all_numeric <- all(str_detect(suffixes, "^\\d+$"))
  if (all_numeric) {
    nums <- sort(as.integer(suffixes))
    collapse_runs(nums)
  } else {
    # Natural sort: extract leading numeric part for sorting, then alpha remainder
    suffixes <- suffixes[order(natural_sort_key(suffixes))]
    paste(suffixes, collapse = ",")
  }
}

# Format a single cluster into a string.
# Bare stem if single symbol with no suffix.
# Stem(formatted_suffixes) otherwise.
format_cluster <- function(cluster) {
  fmt <- format_suffixes(cluster$suffixes)
  if (is.null(fmt)) {
    cluster$stem
  } else {
    paste0(cluster$stem, "(", fmt, ")")
  }
}

# Main entry point: build a HOG display label.
#
# Arguments:
#   syms     - character vector of gene symbols in the HOG (for the focal species)
#   n_sp1    - number of unique genes for species 1
#   n_sp2    - number of unique genes for species 2
#   sp1_abb  - single-letter or short abbreviation for species 1 (default "m")
#   sp2_abb  - single-letter or short abbreviation for species 2 (default "b")
#   min_prefix - minimum LCP length for clustering (default 3)
#
# Returns: a single string, e.g. "m5b1_Nr1i(2,3).Thr(a,b).Vdr"
format_hog_label <- function(
  syms,
  n_sp1,
  n_sp2,
  sp1_abb = "m",
  sp2_abb = "b",
  min_prefix = 3L
) {
  syms <- unique(na.omit(syms))

  if (length(syms) == 0) {
    return(paste0(sp1_abb, n_sp1, sp2_abb, n_sp2))
  }

  clusters <- cluster_symbols(syms, min_prefix = min_prefix)
  formatted <- vapply(clusters, format_cluster, character(1))
  formatted <- sort(formatted)
  symbol_part <- paste(formatted, collapse = ".")

  paste0(sp1_abb, n_sp1, sp2_abb, n_sp2, "_", symbol_part)
}
