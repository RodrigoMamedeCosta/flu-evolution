## ----setup, include=FALSE------------

library(runner)
library(tidyverse)
library(igraph)


## ------------------------------------
############### ROLLING SPHERE #############

#This updated version does performs the randomizations at the raw data (reads) levels, instead of the sphere means. The results are the same, showing their robustness.

### 1) READ IN DISTANCE MATRICES ###
HAmat <- read.csv("Input Data/dm-HAtrimer.csv")
M1mat <- read.csv("Input Data/dm-M1.csv")
M2mat <- read.csv("Input Data/dm-M2.csv")
NPmat <- read.csv("Input Data/dm-NP.csv")
NRmat <- read.csv("Input Data/dm-NR.csv")
NS1mat <- read.csv("Input Data/dm-NS1dimer.csv")
POLmat <- read.csv("Input Data/dm-POL-6QNW.csv") #6RR7 deprecated, missing some C-terminus aminoacids.


### 2) READ IN ANC/DER COUNTS ###
csv_file_path <- "Input Data/residue.ancdercounts.R10.NOINDELS.csv" 
allcounts <- read.csv(csv_file_path)

      #RECOMMENDED: SET MINIMUM COVERAGE TO 1
      #This will convert values below 1 into 1. Each converted site just looks like ancestral so it biases but to a more conservative number (derived frequency = 0) rather than omitting a sample.
      # Function to update values in .Anc columns to remove ancestral zeros
      update_min_value <- function(allcounts) {
        # Identify columns ending with .Anc and .Der
        anc_columns <- grep("\\.Anc$", colnames(allcounts), value = TRUE)
        der_columns <- gsub("\\.Anc$", ".Der", anc_columns)
        
        # Update values in .Anc columns
        for (i in seq_along(anc_columns)) {
          anc_col <- anc_columns[i]
          der_col <- der_columns[i]
          total_sum <- allcounts[[anc_col]] + allcounts[[der_col]]
          allcounts[[anc_col]] <- ifelse(total_sum < 1, 1, allcounts[[anc_col]])
        }
        
        return(allcounts)
      }

# Apply the function to the DataFrame
allcounts <- update_min_value(allcounts)

### 3) DEFINE PROTEINS
protnms <- c("HA", "M1", "M2", "NP", "NR", "NS1", "POL")

### 4) DEFINE PARAMETERS
cdist_values <- c(6,12)    # Sphere radii

qtile  <- 0.999 #Quantile

covmin <- 20 # Minimum MEDIAN coverage across segment (sample excluded from analysis of that protein)

reps <- 20 # Number of randomizations

### 5) FIX COLUMN NAMES ###
allcounts <- allcounts %>%
  mutate(HostStrain = case_when(
    HostStrain == "BALB/C" ~ "BALBC",
    HostStrain == "C57BL/6" ~ "C57BL6",
    TRUE ~ HostStrain
  )) %>% 
  mutate(Sex = case_when(
    Sex == "F" ~ "FEMALE",
    Sex == "M" ~ "MALE",
    TRUE ~ Sex
  ))

### 6) DEFINE GROUPS ###
hosts       <- unique(allcounts$Host)        # e.g. "BALB", "BL6", etc.
hoststrains <- unique(allcounts$HostStrain)  # e.g. "BALBC", "C57BL6"
sexes       <- unique(allcounts$Sex)         # "FEMALE", "MALE"
treatments  <- unique(allcounts$Treatment)   # 1X, 2X, 4X

  categories <- list(
  Host       = hosts,
  HostStrain = hoststrains,
  Sex        = sexes,
  Treatment  = treatments)

### EXAMPLE: run only "AllHosts" minus a custom group = "treatments"
groups <- c("AllHosts", hosts, hoststrains, sexes)

#groups <- c("AllHosts", hosts, hoststrains, sexes,treatments)

### 7) DETECT DATA TYPE ###
data_type <- ifelse(grepl("ALLDATA", csv_file_path, ignore.case = TRUE), "All_Data", "No_Indels")

### 8)RESET DFS ###
  # Initialize/reset data frames
  results_df <- data.frame(
    group    = character(),
    quantile = numeric(),
    cdist    = numeric(),
    site     = character(),
    stringsAsFactors = FALSE
  )
  
  removed_summary <- data.frame(
    cdist    = numeric(),
    group    = character(),
    protnm   = character(),
    SampleID = character(),
    stringsAsFactors = FALSE
  )

   all_results <- list()  

### 8) START LOOPS  ###################################################

### 8) START LOOPS  ###################################################

## Helper: group subsetting
subset_indices_for_group <- function(df, group, hosts, hoststrains, sexes, treatments) {
  if (group == "AllHosts") {
    seq_len(nrow(df))
  } else if (group %in% df$Host) {
    which(df$Host == group)
  } else if (group %in% hoststrains) {
    which(df$HostStrain == group)
  } else if (group %in% sexes) {
    which(df$Sex == group)
  } else if (group %in% treatments) {
    which(df$Treatment == group)
  } else {
    integer(0)
  }
}

## Helpers: structural id and residue id from distance matrix names
## Examples:
##   "E.M1.002.SER" -> struct id "E.M1.002", residue id "M1.002"
##   "A.NP.021.ASN" -> "A.NP.021", "NP.021"
##   "A.PA.001.MET" -> "A.PA.001", "PA.001"
extract_struct_id <- function(x) {
  x2 <- gsub("\\.+", ".", x)
  parts <- strsplit(x2, "\\.")[[1]]
  parts <- parts[nzchar(parts)]
  if (length(parts) >= 3) {
    paste(parts[1], parts[2], parts[3], sep = ".")
  } else {
    x2
  }
}

extract_residue_id <- function(x) {
  x2 <- gsub("\\.+", ".", x)
  parts <- strsplit(x2, "\\.")[[1]]
  parts <- parts[nzchar(parts)]
  if (length(parts) >= 3 && nchar(parts[1]) == 1) {
    paste(parts[2], parts[3], sep = ".")
  } else if (length(parts) >= 2) {
    paste(parts[1], parts[2], sep = ".")
  } else {
    x2
  }
}

## Helper: normalize count column names to "SEG.POS.Suffix"
##   "HA.025.M.Anc" -> "HA.025.Anc"
##   "NP.021.A.Anc" -> "NP.021.Anc"
normalize_count_colname <- function(nm) {
  if (!grepl("\\.(Anc|Der)$", nm)) return(nm)
  suffix <- sub(".*\\.(Anc|Der)$", "\\1", nm)
  core   <- sub("\\.(Anc|Der)$", "", nm)
  core2  <- gsub("\\.+", ".", core)
  parts  <- strsplit(core2, "\\.")[[1]]
  parts  <- parts[nzchar(parts)]
  parts  <- parts[!parts %in% c("Anc", "Der")]
  if (length(parts) >= 3 && nchar(parts[1]) == 1) {
    parts <- parts[-1]
  }
  if (length(parts) >= 2) {
    prefix <- paste(parts[1], parts[2], sep = ".")
  } else {
    prefix <- core2
  }
  paste0(prefix, ".", suffix)
}

normalize_count_columns <- function(df) {
  cn <- colnames(df)
  new_cn <- vapply(cn, normalize_count_colname, character(1))
  colnames(df) <- new_cn
  df
}

## Helper: randomize residue assignment within each sample
## justspace: columns ordered [site1.Anc, site1.Der, site2.Anc, site2.Der, ...]
randomize_dataset <- function(justspace, sites) {
  out <- justspace
  for (i in seq_len(nrow(out))) {
    shuf <- sample(sites)
    newcols <- as.vector(rbind(paste0(shuf, ".Anc"),
                               paste0(shuf, ".Der")))
    out[i, ] <- justspace[i, newcols, drop = FALSE]
  }
  out
}

## Helper: spatial neighborhood on full structural ids
compute_spacedist_full <- function(spacemat, structIDs, cdist) {
  n <- length(structIDs)
  spacedist <- vector("list", n)
  names(spacedist) <- structIDs
  for (j in seq_len(n)) {
    v <- spacemat[, j + 1]          # distances from residue j to all rows
    neigh_idx <- which(v < cdist)
    spacedist[[j]] <- structIDs[neigh_idx]
  }
  spacedist
}

## Helper: compute stats at residue level using full structural neighborhoods
compute_spaced_stats <- function(spaceA, spaceD,
                                 spacedist_full,
                                 struct_to_residue) {
  # spaceA / spaceD:
  #   rows = samples
  #   cols = residue ids like "M1.002"
  # spacedist_full:
  #   named list, names = structural ids like "E.M1.002"
  # struct_to_residue:
  #   named vector, names = structural ids, values = residue ids

  sites <- colnames(spaceA)  # residue-level ids
  m <- nrow(spaceA)
  lfrac <- if (m > 0) 1 / sqrt(m) else 0

  stat1 <- stat2 <- stat3 <- numeric(length(sites))
  names(stat1) <- names(stat2) <- names(stat3) <- sites

  # map residue id -> all structural ids that collapse to it
  residue_to_struct <- split(names(struct_to_residue),
                             struct_to_residue)

  for (resID in sites) {
    structIDs_site <- residue_to_struct[[resID]]
    if (is.null(structIDs_site) || length(structIDs_site) == 0) next

    # gather all structural neighbors across chains
    neigh_full <- unique(unlist(spacedist_full[structIDs_site]))
    if (length(neigh_full) == 0) next

    # convert structural neighbors to residue ids
    neigh_res <- unique(struct_to_residue[neigh_full])
    neigh_res <- intersect(neigh_res, sites)
    if (length(neigh_res) == 0) next

    sumA <- rowSums(spaceA[, neigh_res, drop = FALSE])
    sumD <- rowSums(spaceD[, neigh_res, drop = FALSE])

    p <- sumD / (sumA + sumD)
    p[is.nan(p)] <- 0
    ord <- order(p)
    s <- m:1

    stat1[resID] <- mean(p)
    stat2[resID] <- sum(s * p[ord]) / sum(s)

    # safe logistic regression
    valid_for_glm <- (
      sum(sumA) > 0 &&
      sum(sumD) > 0 &&
      any(p > 0 & p < 1) &&
      length(unique(p)) > 1
    )

    if (valid_for_glm) {
      tmp <- data.frame(a = sumA, d = sumD, p = p, s = s)

      lg <- suppressWarnings(
        glm(cbind(d, a) ~ s, data = tmp, family = binomial)
      )

      pred <- suppressWarnings(
        predict(lg, newdata = data.frame(s = lfrac * m))
      )

      prob <- exp(pred) / (1 + exp(pred))
      if (!is.finite(prob)) prob <- 0
      stat3[resID] <- prob
    } else {
      stat3[resID] <- 0
    }
  }

  data.frame(site = sites,
             stat1 = stat1,
             stat2 = stat2,
             stat3 = stat3,
             stringsAsFactors = FALSE)
}


cspacedat_group <- list()
crspacedat_group <- list()

## Main loop: per protein
## Main loop: per protein
for (protnm in protnms) {
  cat("\n========== Processing protein:", protnm, "==========\n")

  ## 1. Structural matrix and mapping
  spacemat <- get(paste0(protnm, "mat"))

  structIDs_full  <- vapply(spacemat$Residue, extract_struct_id,  character(1))
  residueIDs_full <- vapply(spacemat$Residue, extract_residue_id, character(1))
  struct_to_residue <- setNames(residueIDs_full, structIDs_full)

  ## 2. Normalize counts once
  counts_norm <- normalize_count_columns(allcounts)

  ## 3. Extract counts only and restrict to residues present in structure
  counts_only_all <- counts_norm[, -(1:9), drop = FALSE]
  if (ncol(counts_only_all) == 0) {
    cat("No count columns found for", protnm, "\n")
    next
  }

  site_prefix_counts <- sub("\\.(Anc|Der)$", "", colnames(counts_only_all))
  keep_idx_cols <- which(site_prefix_counts %in% unique(residueIDs_full))
  if (length(keep_idx_cols) == 0) {
    cat("No matching residue ids between counts and structure for", protnm, "\n")
    next
  }

  counts_only_global <- counts_only_all[, keep_idx_cols, drop = FALSE]

  anc_cols <- grep("\\.Anc$", colnames(counts_only_global), value = TRUE)
  der_cols <- gsub("\\.Anc$", ".Der", anc_cols)

  spaceA_global <- counts_only_global[, anc_cols, drop = FALSE]
  spaceD_global <- counts_only_global[, der_cols, drop = FALSE]

  sites <- sub("\\.(Anc|Der)$", "", anc_cols)   # residue ids like "M1.002"
  colnames(spaceA_global) <- sites
  colnames(spaceD_global) <- sites

  ## 4. Build justspace_global with interleaved Anc/Der for randomization
  new_order <- as.vector(rbind(paste0(sites, ".Anc"),
                               paste0(sites, ".Der")))
  justspace_global <- counts_only_global[, new_order, drop = FALSE]

  ## 5. For each radius: compute structural spacedist, then per group
  for (cdist in cdist_values) {
    cat("  Sphere radius:", cdist, "Angstrom\n")

    spacedist_full <- compute_spacedist_full(spacemat, structIDs_full, cdist)

    for (group in groups) {
      idx_all <- subset_indices_for_group(counts_norm, group,
                                          hosts, hoststrains, sexes, treatments)
      if (length(idx_all) == 0) next

      ## 5a. Coverage filtering for this group, using real data only
      dat_real_pre <- justspace_global[idx_all, , drop = FALSE]

      tmpA_pre <- dat_real_pre[, grep("\\.Anc$", colnames(dat_real_pre)), drop = FALSE]
      tmpD_pre <- dat_real_pre[, grep("\\.Der$", colnames(dat_real_pre)), drop = FALSE]
      colnames(tmpA_pre) <- sites
      colnames(tmpD_pre) <- sites

      sample_cov <- apply(tmpA_pre + tmpD_pre, 1, median, na.rm = TRUE)
      good <- which(sample_cov >= covmin)
      bad  <- which(sample_cov < covmin)

      if (length(bad) > 0) {
        removed_sample_ids <- counts_norm$SampleID[idx_all[bad]]
        removed_df <- data.frame(
          cdist    = cdist,
          group    = group,
          protnm   = protnm,
          SampleID = removed_sample_ids,
          stringsAsFactors = FALSE
        )
        removed_summary <- rbind(removed_summary, removed_df)
        cat("Removed", length(bad), "samples for", protnm,
            "in group", group, "due to median coverage <", covmin, "\n")
      }

      if (length(good) == 0) {
        cat("No samples kept for", protnm, "in group", group, "\n")
        next
      }

      idx_keep <- idx_all[good]
      cat("Kept", length(idx_keep), "samples for", protnm, "in group", group, "\n")

      ## 5b. Build per group real dataset and randomizations after filtering
      dat_real <- justspace_global[idx_keep, , drop = FALSE]

      randomized_list_group <- vector("list", reps + 1)
      randomized_list_group[[1]] <- dat_real
      for (r in 1:reps) {
        randomized_list_group[[r + 1]] <- randomize_dataset(dat_real, sites)
      }

      key <- paste(group, cdist, sep = "_")

      ## 5c. Real data stats
      dat_real_use <- randomized_list_group[[1]]
      tmpA <- dat_real_use[, grep("\\.Anc$", colnames(dat_real_use)), drop = FALSE]
      tmpD <- dat_real_use[, grep("\\.Der$", colnames(dat_real_use)), drop = FALSE]
      colnames(tmpA) <- sites
      colnames(tmpD) <- sites

      spacedat_real <- compute_spaced_stats(tmpA, tmpD,
                                            spacedist_full,
                                            struct_to_residue)
      spacedat_real$protnm <- protnm

      cspacedat_group[[key]] <- rbind(cspacedat_group[[key]], spacedat_real)

      ## 5d. Randomization stats
      allrspacedat <- data.frame()

      for (r in 1:reps) {
        dat_rand <- randomized_list_group[[r + 1]]
        tmpA_r <- dat_rand[, grep("\\.Anc$", colnames(dat_rand)), drop = FALSE]
        tmpD_r <- dat_rand[, grep("\\.Der$", colnames(dat_rand)), drop = FALSE]
        colnames(tmpA_r) <- sites
        colnames(tmpD_r) <- sites

        rspaced <- compute_spaced_stats(tmpA_r, tmpD_r,
                                        spacedist_full,
                                        struct_to_residue)
        rspaced$protnm <- protnm
        rspaced <- subset(rspaced, stat1 + stat2 + stat3 > 0)
        allrspacedat <- rbind(allrspacedat, rspaced)
      }

      crspacedat_group[[key]] <- rbind(crspacedat_group[[key]], allrspacedat)

      cat("    Finished", protnm, "in", group, "at", cdist, "Angstrom\n")
    }  # end group loop
  }    # end cdist loop
}      # end protein loop


## 8. Final analysis and output per group and radius

full_output_path <- "Outputs/Raw_Outputs/Spatial_outputs/"
if (!dir.exists(full_output_path)) dir.create(full_output_path, recursive = TRUE)

for (group in groups) {
  for (cdist in cdist_values) {
    key <- paste(group, cdist, sep = "_")
    cspacedat <- cspacedat_group[[key]]
    crspacedat <- crspacedat_group[[key]]
    
    if (is.null(cspacedat) || is.null(crspacedat) ||
        nrow(cspacedat) == 0 || nrow(crspacedat) == 0) {
      next
    }
    
    ## Write raw outputs
    output_file_path_spacedat <- paste0(full_output_path,
                                        group, "_", cdist, "A_", qtile, "_spacedat.csv")
    write.csv(cspacedat, output_file_path_spacedat, row.names = FALSE)
    
    output_file_path_randomspacedat <- paste0(full_output_path,
                                              group, "_", cdist, "A_", qtile, "_randomspacedat.csv")
    write.csv(crspacedat, output_file_path_randomspacedat, row.names = FALSE)
    
    ## Hotspot calls
    Protspace1 <- cspacedat$site[
      cspacedat$stat1 > quantile(crspacedat$stat1, qtile, na.rm = TRUE)
    ]
    Protspace2 <- cspacedat$site[
      cspacedat$stat2 > quantile(crspacedat$stat2, qtile, na.rm = TRUE)
    ]
    Protspace3 <- cspacedat$site[
      cspacedat$stat3 > quantile(crspacedat$stat3, qtile, na.rm = TRUE)
    ]
    
    Protspace1 <- unique(Protspace1)
    Protspace2 <- unique(Protspace2)
    Protspace3 <- unique(Protspace3)
    
    ## HA signal peptide correction
    Protspace_list <- list(Protspace1 = Protspace1,
                           Protspace2 = Protspace2,
                           Protspace3 = Protspace3)
    
    Protspace_adjusted_list <- list()
    
    for (i in seq_along(Protspace_list)) {
      v <- Protspace_list[[i]]
      Protspace_df <- data.frame(
        Prefix = sapply(v, function(x) strsplit(x, "\\.")[[1]][1]),
        Number = sapply(v, function(x) strsplit(x, "\\.")[[1]][2]),
        stringsAsFactors = FALSE
      )
      Protspace_df$Number <- as.numeric(Protspace_df$Number)
      Protspace_df$AdjustedNumber <- ifelse(Protspace_df$Prefix == "HA",
                                            Protspace_df$Number - 16,
                                            Protspace_df$Number)
      Protspace_adjusted <- paste(
        Protspace_df$Prefix,
        sprintf("%03d",
                ifelse(is.na(Protspace_df$AdjustedNumber),
                       Protspace_df$Number,
                       Protspace_df$AdjustedNumber)),
        sep = "."
      )
      Protspace_adjusted_list[[names(Protspace_list)[i]]] <- Protspace_adjusted
    }
    
    Protspace1 <- Protspace_adjusted_list$Protspace1
    Protspace2 <- Protspace_adjusted_list$Protspace2
    Protspace3 <- Protspace_adjusted_list$Protspace3
    
    ## Choose stat2 as the main result
    results <- Protspace2
    if (!exists("results") || is.null(results)) {
      results <- character(0)
    }
    
    all_results[[group]] <- results
    
    if (length(results) > 0) {
      results_temp_df <- data.frame(
        group    = rep(group, length(results)),
        quantile = rep(qtile, length(results)),
        cdist    = rep(cdist, length(results)),
        site     = results,
        stringsAsFactors = FALSE
      )
      results_df <- rbind(results_df, results_temp_df)
    }
    
    rm(results)
    
    cat("\nSaved data for group:", group,
        "for quantile =", qtile,
        "with sphere radius =", cdist,
        "minimum median coverage =", covmin, "\n\n")
  }
}


##################################### END OF DATA ANALYSIS AND RANDOMIZATIONS ###############################

   
### 10) AFTER ALL GROUPS/POLYMERS are processed, do final intersections 

if ("AllHosts" %in% names(all_results)) {
  # Make a directory to hold intersections
  if (!dir.exists("Outputs/Intersect_Outputs")) dir.create("Outputs/Intersect_Outputs")
  
  for (g in names(all_results)) {
    if (g == "AllHosts") next  # skip self

    intersect_with_all <- intersect(all_results[[g]], all_results[["AllHosts"]])
    
    
    #name the csv
    output_file <- sprintf("Outputs/Intersect_Outputs/%s_vs_AllHosts.csv", g)
    
    # Write the CSV
    write.csv(intersect_with_all, output_file, row.names = FALSE, quote = TRUE)
    
    cat("Intersected", g, "with AllHosts. Found", length(intersect_with_all), "common sites.\n")
  }
}

### 11) OPTIONAL: If doing “shared by multiple groups” logic,
### intersections can be combined using combn(...):

# Example: find sites shared by *all* groups (besides "AllHosts"):
other_groups <- setdiff(names(all_results), "AllHosts")
if (length(other_groups) > 1) {
  shared_all <- Reduce(intersect, all_results[other_groups])
  cat("\nSites shared by all groups:\n", paste(shared_all, collapse=", "), "\n")
}

### 12) WRITE OUT THE results df
 
# Define the output path for the CSV file
  csv_output_path <- paste0("Outputs/table_output.csv")
  
  # Write the results_df data frame to a CSV file
  write.csv(results_df, csv_output_path, row.names = FALSE)
  
### 13) WRITE OUT REMOVED SUMMARY ###
  write.csv(removed_summary, paste0("Outputs/removed_samples.csv"), row.names = FALSE)
  
  cat("\nAnalysis complete\n")



###############################################################################
### 14) WITHIN-CATEGORY COMPARISONS (UNIQUE VS. SHARED AMONG HOSTS, STRAINS, ETC.)
###############################################################################

# Create an output directory for the within-group reports:
dir.create("Outputs/WithinGroup_Comparisons", showWarnings = FALSE, recursive = TRUE)


# 2) For each category (e.g. "Host"), gather the subgroups (like c("Dog","Cat","Human"))
#    from all_results. Then do the shared/unique logic, and write a report.

for (cat_name in names(categories)) {
  # The subgroups for this category:
  subgroups <- categories[[cat_name]]  
  # e.g. c("MALE","FEMALE") if cat_name=="Sex"
  
  # Filter out any subgroups that might not appear in `all_results` 
  # (in case some are missing coverage or were skipped).
  subgroups <- intersect(subgroups, names(all_results))
  
  # If there's only 0 or 1 subgroup found, no need for "within-group" comparisons:
  if (length(subgroups) < 2) {
    cat(sprintf("\n[SKIP] Category '%s' has <2 subgroups in all_results. No within-group comparison.\n", cat_name))
    next
  }
  
  # Extract hits for each subgroup in this category:
  sub_hits <- all_results[subgroups]
  # Now sub_hits is a named list, e.g. sub_hits[["BALBC"]], sub_hits[["C57BL6"]]
  
  # ----------------------------------------------------------------------------
  # Step A: Gather shared sites for all combos
  # ----------------------------------------------------------------------------
  shared_results <- list()
  
  # We'll do combos from largest to smaller. For example, 
  # if subgroups = c("1X","2X","4X"), we do combos of length 3, then combos of length 2.
  for (n in seq(from = length(subgroups), to = 2)) {
    # All combos of size n
    combos <- combn(subgroups, m = n, simplify = FALSE)
    
    for (combo in combos) {
      # Intersection for these subgroups
      shared_vals <- Reduce(intersect, sub_hits[combo])
      
      if (length(shared_vals) > 0) {
        # Remove hits already captured by a previously recorded (larger) combo
        for (existing_combo in names(shared_results)) {
          shared_vals <- setdiff(shared_vals, shared_results[[existing_combo]])
        }
        
        if (length(shared_vals) > 0) {
          combo_name <- paste(combo, collapse = ", ")
          shared_results[[combo_name]] <- shared_vals
        }
      }
    }
  }
  
  # ----------------------------------------------------------------------------
  # Step B: Identify unique sites by excluding all that appear in any shared combos
  # ----------------------------------------------------------------------------
  unique_results <- list()
  already_shared <- unique(unlist(shared_results))
  
  for (sg in subgroups) {
    # Subgroup's hits minus any that appear in shared combos
    unique_vals <- setdiff(sub_hits[[sg]], already_shared)
    unique_results[[sg]] <- unique_vals
  }
  
  # ----------------------------------------------------------------------------
  # Step C: Write the report to a text file
  # ----------------------------------------------------------------------------
  report_path <- file.path("Outputs","WithinGroup_Comparisons", 
                           paste0(cat_name, "_WithinGroup_Report.txt"))
  
  sink(report_path)
  cat("=== WITHIN-GROUP REPORT FOR CATEGORY:", cat_name, "===\n\n")
  cat("Subgroups:", paste(subgroups, collapse=", "), "\n\n")
  
  cat("UNIQUE SITES:\n\n")
  for (sg in subgroups) {
    cat(sprintf("Unique to %s:\n", sg))
    if (length(unique_results[[sg]]) == 0) {
      cat("  None\n\n")
    } else {
      cat("  ", paste(unique_results[[sg]], collapse=", "), "\n\n")
    }
  }
  
  cat("SHARED SITES:\n\n")
  for (combo_name in names(shared_results)) {
    cat(sprintf("Shared by %s:\n", combo_name))
    cat("  ", paste(shared_results[[combo_name]], collapse=", "), "\n\n")
  }
  
  sink()  # close the report file
  
  # Print message to console
  cat(sprintf("\n[INFO] Wrote within-group report for '%s' to:\n  %s\n\n", 
              cat_name, report_path))
}

cat("\nAll within-category comparisons completed.\n\n")


cat("\nRolling Sphere Complete.\n\n")


## ------------------------------------
############ HOTSPOT NEXUS (MEDOID) — SIMPLE, FREQ & RECURRENCE ################
# This section:
# 1) Reads in significant centers detected by the rolling sphere.
# 2) Uses site-level distance matrices to cluster nearby centers.
# 3) Builds the union of all residues within the radius of any center in a cluster.
# 4) Assigns per-site weights (frequency or recurrence).
# 5) Picks a medoid (geometric center, optionally weighted by mutation frequency).
# 6) Identifies the driver site and hotspot category (single / multiple / dispersed).
# 7) Defines contributing sites as union members with non-zero weight for that group.
# 8) Produces summary tables and a clusters matrix with per-host hotspot types and frequencies.
#################################################################################

## ====================== High level parameters and toggles =====================

single_driver_threshold <- 0.75   # driver share > 0.70 => "single", otherwise "dispersed"
use_jaccard             <- FALSE  # optionally connect centers with high neighborhood overlap
jaccard_cut             <- 0.50   # Jaccard overlap threshold if use_jaccard = TRUE
out_dir_base            <- "Outputs/Summary Tables"

# Weighting mode for hotspot summaries:
#   "freq"       = group mean frequency per site (from freqs_grp)
#   "recurrence" = per-sample mean frequency per site (from freqs_per_sample) after exclusions
weight_mode <- "recurrence"

# Medoid definition:
#   TRUE  = medoid minimizes weighted sum of distances (weights from chosen weight_mode)
#   FALSE = medoid minimizes unweighted sum of distances (purely geometric center)
use_weighted_medoid <- TRUE

# Sample IDs to exclude from per-sample recurrence calculations when duplicates exist
EXCLUDE_SAMPLE_IDS <- c("9770", "97700", "10318")

## ========================== Lightweight convenience ===========================

fmt3           <- function(n) sprintf("%03d", as.integer(n))
num_from_cdist <- function(x) as.numeric(gsub("[^0-9.]+", "", as.character(x)))
fmt_lab        <- function(r) ifelse(abs(r - round(r)) < 1e-9, as.character(as.integer(r)), as.character(r))
site_protein   <- function(s) sub("\\..*$", "", s)

# Assign hotspot category based on driver share thresholds
hotspot_type_from_share <- function(share,
                                    single_thr  = single_driver_threshold) {
  dplyr::case_when(
    is.na(share)        ~ NA_character_,
    share > single_thr  ~ "single",
    TRUE                ~ "dispersed"
  )
}

## ============================== Load inputs ===================================

# Significant centers from rolling-sphere (table_output.csv from upstream code)
sig_path <- "Outputs/table_output.csv"
if (!file.exists(sig_path)) stop("Missing significant centers at: ", sig_path)
sig_df <- readr::read_csv(sig_path, show_col_types = FALSE)
stopifnot(all(c("site", "group", "cdist") %in% names(sig_df)))

# Normalize cdist and add protein
sig_df <- dplyr::mutate(
  sig_df,
  protein   = site_protein(site),
  cdist_num = num_from_cdist(cdist),
  cdist     = cdist_num
)

groups       <- sort(unique(sig_df$group))
cdist_values <- sort(unique(sig_df$cdist))

# Group mean frequencies per site (used in "freq" mode)
freqs_grp_path <- "Input Data/Mean nonsyn variant frequency by strain.csv"
if (!file.exists(freqs_grp_path)) stop("Missing group freqs at: ", freqs_grp_path)
freqs_grp <- read.csv(freqs_grp_path, check.names = FALSE) %>%
  dplyr::mutate(site = sprintf("%s.%03d", Protein, as.integer(AA.pos)))

# Per-sample frequencies per site (used in "recurrence" mode)
recur_path <- "Input Data/freqs.allsamples.csv"
if (!file.exists(recur_path)) stop("Missing per-sample freqs at: ", recur_path)
freqs_per_sample <- readr::read_csv(recur_path, show_col_types = FALSE) %>%
  dplyr::mutate(SampleID = as.character(SampleID)) %>%
  dplyr::filter(!(SampleID %in% EXCLUDE_SAMPLE_IDS))

# Mice to exclude for recurrence mode (group × protein family × radius)
removed_samples <- tryCatch(
  readr::read_csv("Outputs/removed_samples.csv", show_col_types = FALSE),
  error = function(e) tibble::tibble(
    cdist   = numeric(),
    group   = character(),
    protnm  = character(),
    SampleID = character()
  )
)

## ====================== Load and standardize distance maps ====================

load_dm_if_missing <- function(obj_name, file_hint) {
  if (!exists(obj_name, inherits = TRUE)) {
    path <- file.path("Input Data", file_hint)
    if (file.exists(path)) {
      assign(obj_name, read.csv(path, check.names = FALSE), envir = .GlobalEnv)
    }
  }
}

load_dm_if_missing("HAmat",  "dm-HAtrimer.csv")
load_dm_if_missing("M1mat",  "dm-M1.csv")
load_dm_if_missing("M2mat",  "dm-M2.csv")
load_dm_if_missing("NPmat",  "dm-NP.csv")
load_dm_if_missing("NRmat",  "dm-NR.csv")
load_dm_if_missing("NS1mat", "dm-NS1dimer.csv")
load_dm_if_missing("POLmat", "dm-POL-6QNW.csv")

# Coerce distance matrices to numeric with aligned row/column labels
coerce_dm <- function(X) {
  if (is.data.frame(X)) {
    rn <- as.character(X[[1]])
    M  <- suppressWarnings(as.matrix(X[, -1, drop = FALSE]))
    storage.mode(M) <- "double"
    rownames(M) <- rn
    if (!is.null(colnames(X)) && length(colnames(X)) == ncol(M) + 1) {
      colnames(M) <- colnames(X)[-1]
    }
    if (is.null(colnames(M))) colnames(M) <- rownames(M)
    M
  } else if (is.matrix(X)) {
    M <- X
    storage.mode(M) <- "double"
    if (is.null(rownames(M)) && !is.null(colnames(M))) rownames(M) <- colnames(M)
    if (is.null(colnames(M))) colnames(M) <- rownames(M)
    M
  } else {
    NULL
  }
}

# Convert atom-level labels to site labels (e.g. "A.NS1.183.GLU" → "NS1.183")
raw_to_site <- function(v) {
  v <- sub("^\\w*\\.", "", v)
  sub("\\.[^.]*$", "", v)
}

# Collapse atom-level matrix to site-level by minimum atom–atom distance
collapse_to_site_dm <- function(M) {
  if (is.null(M)) return(NULL)
  r_site <- raw_to_site(rownames(M))
  c_site <- raw_to_site(colnames(M))
  sites  <- sort(unique(c(r_site, c_site)))
  if (!length(sites)) return(NULL)

  DM <- matrix(NA_real_, length(sites), length(sites),
               dimnames = list(sites, sites))
  r_grp <- split(seq_len(nrow(M)), r_site)
  c_grp <- split(seq_len(ncol(M)), c_site)

  for (si in sites) {
    ri <- r_grp[[si]]
    if (is.null(ri)) next
    for (sj in sites) {
      cj <- c_grp[[sj]]
      if (is.null(cj)) next
      DM[si, sj] <- suppressWarnings(min(M[ri, cj, drop = FALSE], na.rm = TRUE))
    }
  }
  pmin(DM, t(DM), na.rm = TRUE)
}

mat_names <- ls(pattern = "mat$")
if (!length(mat_names)) stop("No distance matrices in memory (objects ending with 'mat').")
raw_objs <- mget(mat_names, inherits = TRUE)

distance_mats <- list()
for (nm in names(raw_objs)) {
  M0 <- coerce_dm(raw_objs[[nm]])
  if (is.null(M0)) next
  DM <- collapse_to_site_dm(M0)
  if (is.null(DM)) next
  key <- sub("mat$", "", nm)
  key <- sub("\\.$", "", key)
  distance_mats[[key]] <- DM[rownames(DM), rownames(DM), drop = FALSE]
}

keep_keys <- c("HA", "M1", "M2", "NP", "NR", "NS1", "POL")
distance_mats <- distance_mats[names(distance_mats) %in% keep_keys]
if (!length(distance_mats)) stop("No usable protein distance matrices after collapse.")

# Adjust HA numbering to post-cleavage coordinates (subtract 16) so it matches site labels
if ("HA" %in% names(distance_mats)) {
  DM <- distance_mats[["HA"]]
  rn <- rownames(DM)
  is_ha <- grepl("^HA\\.\\d{3}$", rn)
  pos   <- as.integer(sub("^HA\\.(\\d{3})$", "\\1", rn))
  pos_new <- ifelse(is_ha, pos - 16L, NA_integer_)
  rn_new <- rn
  rn_new[is_ha] <- sprintf("HA.%03d", pos_new[is_ha])
  keep <- !(is_ha & pos_new < 1L)
  DM   <- DM[keep, keep, drop = FALSE]
  rownames(DM) <- rn_new[keep]
  colnames(DM) <- rn_new[keep]
  distance_mats[["HA"]] <- DM
}

## ======================= Weighting schemes for sites ==========================

# Group mean frequency per site ("freq" mode)
weight_freq <- function(protein, sites, grp) {
  gcol <- if (grp %in% names(freqs_grp)) grp else {
    nn  <- gsub("[^A-Za-z0-9]+", "", tolower(names(freqs_grp)))
    idx <- which(nn == gsub("[^A-Za-z0-9]+", "", tolower(grp)))[1]
    if (is.na(idx)) return(tibble::tibble(site = sites, w = 0))
    names(freqs_grp)[idx]
  }

  freqs_grp %>%
    dplyr::filter(Protein == protein, site %in% sites) %>%
    dplyr::transmute(site, w = .data[[gcol]]) %>%
    dplyr::right_join(tibble::tibble(site = sites), by = "site") %>%
    dplyr::mutate(w = dplyr::if_else(is.na(w), 0, w))
}

# Map aggregate groups to underlying single-host groups
expand_groups <- function(grp) {
  switch(
    grp,
    "BALBF"    = c("BALBF"),
    "BALBM"    = c("BALBM"),
    "BL6F"     = c("BL6F"),
    "BL6M"     = c("BL6M"),
    "BALBC"    = c("BALBF", "BALBM"),
    "C57BL6"   = c("BL6F", "BL6M"),
    "FEMALE"   = c("BALBF", "BL6F"),
    "MALE"     = c("BALBM", "BL6M"),
    "AllHosts" = c("BALBF", "BALBM", "BL6F", "BL6M"),
    grp
  )
}

# Per-sample mean frequency per site ("recurrence" mode), excluding low-coverage mice
weight_recurrence <- function(protein, sites, grp, r_num) {
  if (!length(sites)) return(tibble::tibble(site = character(), w = numeric()))
  prot_key <- if (protein %in% c("PB1", "PB2", "PA")) "POL" else protein

  df <- freqs_per_sample %>%
    dplyr::filter(
      site %in% sites,
      if (protein %in% c("PB1", "PB2", "PA")) Protein %in% c("PB1", "PB2", "PA") else Protein == protein,
      Group %in% expand_groups(grp)
    )

  if (nrow(removed_samples)) {
    rmv_ids <- removed_samples %>%
      dplyr::filter(group == grp, protnm == prot_key, as.numeric(cdist) == as.numeric(r_num)) %>%
      dplyr::pull(SampleID) %>%
      as.character()
    if (length(rmv_ids)) {
      df <- df %>% dplyr::filter(!(as.character(SampleID) %in% rmv_ids))
    }
  }
  if (!nrow(df)) return(tibble::tibble(site = sites, w = 0))

  df %>%
    dplyr::group_by(site) %>%
    dplyr::summarise(w = mean(freq, na.rm = TRUE), .groups = "drop") %>%
    dplyr::right_join(tibble::tibble(site = sites), by = "site") %>%
    dplyr::mutate(w = dplyr::if_else(is.na(w), 0, w)) %>%
    dplyr::select(site, w)
}

## ==================== Neighborhood maps and medoid helper ====================

# For each site, list sites within radius r (including the site itself)
build_within_map <- function(DM, r) {
  sites <- rownames(DM)
  A     <- DM <= r & !is.na(DM)
  diag(A) <- TRUE
  setNames(lapply(seq_along(sites), function(i) sites[A[i, ]]), sites)
}

# Select medoid from candidate sites given a distance matrix and an optional weight table
weighted_medoid <- function(DM, candidates, w_tbl, use_weighted = TRUE) {
  U <- intersect(candidates, rownames(DM))
  if (!length(U)) return(NA_character_)

  D <- DM[U, U, drop = FALSE]

  if (!use_weighted) {
    s <- rowSums(D, na.rm = TRUE)
    return(U[which.min(s)])
  }

  w <- w_tbl$w[match(U, w_tbl$site)]
  w[is.na(w)] <- 0
  if (sum(w) <= 0) {
    s <- rowSums(D, na.rm = TRUE)
    return(U[which.min(s)])
  }

  s <- as.numeric(D %*% w)
  U[which.min(s)]
}

## ======================= Core pipeline for one weighting mode =================

compute_annotations_for_mode <- function(weight_mode = c("freq", "recurrence")) {
  weight_mode <- match.arg(weight_mode)

  get_weights <- switch(
    weight_mode,
    "freq" = function(prot, sites, grp, r) weight_freq(prot, sites, grp),
    "recurrence" = function(prot, sites, grp, r) weight_recurrence(prot, sites, grp, r)
  )

  within_cache <- new.env(parent = emptyenv())
  result_rows  <- list()

  for (grp in groups) {
    for (r in cdist_values) {
      sub_df <- dplyr::filter(sig_df, group == grp, cdist == r)
      if (!nrow(sub_df)) next

      for (prot in unique(sub_df$protein)) {
        DM_key <- if (prot %in% c("PB1", "PB2", "PA")) "POL" else prot
        if (!(DM_key %in% names(distance_mats))) next

        DM    <- distance_mats[[DM_key]]
        sub_p <- dplyr::filter(sub_df, protein == prot)
        centers <- intersect(unique(sub_p$site), rownames(DM))
        if (!length(centers)) next

        cache_key <- paste0(DM_key, "|", r)
        nh_map <- if (exists(cache_key, envir = within_cache)) {
          get(cache_key, envir = within_cache)
        } else {
          res <- build_within_map(DM, r)
          assign(cache_key, res, envir = within_cache)
          res
        }

        # Build graph of centers: edges if within radius r; optional Jaccard overlap
        edges <- tibble::tibble(a = character(), b = character())
        if (length(centers) >= 2) {
          A <- DM[centers, centers, drop = FALSE] <= r
          diag(A) <- FALSE
          e <- which(A, arr.ind = TRUE)
          if (nrow(e)) {
            edges <- tibble::tibble(
              a = centers[e[, 1]],
              b = centers[e[, 2]]
            ) %>%
              dplyr::distinct()
          }

          if (use_jaccard) {
            for (i in seq_along(centers)) {
              for (j in (i + 1):length(centers)) {
                si  <- centers[i]
                sj  <- centers[j]
                nsi <- nh_map[[si]]
                nsj <- nh_map[[sj]]
                if (length(nsi) && length(nsj)) {
                  jac <- length(intersect(nsi, nsj)) / length(union(nsi, nsj))
                  if (jac >= jaccard_cut) {
                    edges <- dplyr::bind_rows(
                      edges,
                      tibble::tibble(a = si, b = sj)
                    )
                  }
                }
              }
            }
            edges <- dplyr::distinct(edges)
          }
        }

        g <- igraph::graph_from_data_frame(
          edges,
          directed = FALSE,
          vertices = tibble::tibble(name = centers)
        )
        comps    <- igraph::components(g)
        clusters <- split(names(comps$membership), comps$membership)

        for (k in seq_along(clusters)) {
          members <- clusters[[k]]

          # Union of all residues within radius r of any center in the cluster
          U_union <- unique(unlist(nh_map[members], use.names = FALSE))
          if (!length(U_union)) next

          # Candidate sites for medoid (entire union, independent of weights)
          U_cands <- U_union

          # Weights on union, used for contributing sites, driver and weighted medoid
          w_tbl <- get_weights(prot, U_union, grp, r)

          if (nrow(w_tbl)) {
            top <- w_tbl$site[which.max(w_tbl$w)]
            if (length(top) && !is.na(top) && !(top %in% U_cands)) {
              U_cands <- c(U_cands, top)
            }
          }

          med <- weighted_medoid(
            DM        = DM,
            candidates = U_cands,
            w_tbl      = w_tbl,
            use_weighted = use_weighted_medoid
          )

          # Contributing sites are union members with non-zero weight for this group
          contributing_set <- w_tbl$site[w_tbl$w > 0]
          contributing_set <- intersect(contributing_set, U_union)
          contributing_set <- sort(unique(contributing_set))

          # Driver site and share among contributing sites
          tot <- sum(w_tbl$w[w_tbl$site %in% contributing_set], na.rm = TRUE)
          if (tot > 0) {
            # Restrict to contributing set for driver
            w_contrib  <- w_tbl[w_tbl$site %in% contributing_set, , drop = FALSE]
            top_idx    <- which.max(w_contrib$w)
            driver_site  <- w_contrib$site[top_idx]
            driver_share <- w_contrib$w[top_idx] / tot
          } else {
            driver_site  <- NA_character_
            driver_share <- NA_real_
          }

          hotspot_type <- hotspot_type_from_share(driver_share)

          result_rows[[length(result_rows) + 1]] <- tibble::tibble(
            group           = grp,
            cdist           = r,
            protein         = prot,
            site            = members,
            cluster_id      = sprintf("CLU%s_%s_%s_%03d", fmt_lab(r), prot, grp, k),
            cluster_medoid  = med,
            hotspot_type    = hotspot_type,
            driver_site     = driver_site,
            driver_share    = driver_share,
            contributing_sites = if (length(contributing_set)) {
              paste(contributing_set, collapse = ";")
            } else {
              NA_character_
            }
          )
        }
      }
    }
  }

  cluster_rows <- dplyr::bind_rows(result_rows)

  extra_cols <- setdiff(names(sig_df), c("site", "group", "cdist", "cdist_num", "protein"))
  sig_slim   <- dplyr::select(
    sig_df,
    site, group, cdist, protein,
    dplyr::any_of(extra_cols)
  )

  annotated <- sig_slim %>%
    dplyr::left_join(cluster_rows, by = c("site", "group", "cdist", "protein"))

  minimal <- annotated %>%
    dplyr::select(
      dplyr::any_of(c(
        "group", "quantile", "cdist", "site", "protein",
        "cluster_id", "cluster_medoid", "hotspot_type",
        "driver_site", "driver_share", "contributing_sites"
      ))
    ) %>%
    dplyr::arrange(group, cdist, protein, site)

  list(
    annotated = annotated,
    minimal   = minimal
  )
}

## ============================= Run chosen mode ================================

mode_suffix <- if (weight_mode == "freq") "freq" else "recurrence"
res         <- compute_annotations_for_mode(weight_mode)

ann_tbl     <- res$annotated
minimal_tbl <- res$minimal

dir.create(out_dir_base, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  minimal_tbl,
  file.path(out_dir_base, paste0("significant_centers_minimal-", mode_suffix, ".csv")),
  na = "-"
)

cat("\nSaved minimal significant-center table for mode '", weight_mode, "'.\n", sep = "")

## ========== Summaries by groups + medoid cross-group labels ==================

# Groups to summarise at the cluster level (both single-host and aggregate)
individual_groups_full <- c(
  "BALBF", "BALBM", "BL6F", "BL6M",
  "BALBC", "C57BL6", "FEMALE", "MALE", "AllHosts"
)

# Four primary host groups used for classification labels
primary_lineages <- c("BALBF", "BALBM", "BL6F", "BL6M")

# Compact cluster-level summary per group × radius
make_individual_summary <- function(annot_tbl) {
  annot_tbl %>%
    dplyr::filter(group %in% individual_groups_full) %>%
    dplyr::group_by(group, cdist, protein, cluster_id) %>%
    dplyr::summarise(
      cluster_medoid   = dplyr::first(na.omit(cluster_medoid)),
      hotspot_type     = dplyr::first(na.omit(hotspot_type)),
      driver_site      = dplyr::first(na.omit(driver_site)),
      driver_share     = dplyr::first(na.omit(driver_share)),
      n_centers        = dplyr::n_distinct(site),
      contributing_sites = {
        cs <- unique(na.omit(contributing_sites))
        if (!length(cs)) "-" else paste(sort(unique(unlist(strsplit(cs, ";")))), collapse = ";")
      },
      .groups = "drop"
    ) %>%
    dplyr::arrange(group, cdist, protein, cluster_id)
}

indiv_summary <- make_individual_summary(ann_tbl)

# Classification helper based on presence in the four primary lineages
classify_groups <- function(gs) {
  gs <- sort(unique(na.omit(gs)))
  k  <- length(gs)
  if (k == 0) return("-")
  if (k >= 3)                               return("mouse-specific")
  if (identical(gs, c("BALBF", "BL6F")))    return("sex-specific (female)")
  if (identical(gs, c("BALBM", "BL6M")))    return("sex-specific (male)")
  if (identical(gs, c("BALBF", "BALBM")))   return("genotype-specific (BALBC)")
  if (identical(gs, c("BL6F", "BL6M")))     return("genotype-specific (C57BL6)")
  if (k == 1)                               return("host-specific")
  "other"
}

# Presence map of medoids across the four primary lineages (any radius)
make_presence_map <- function(annot_tbl) {
  annot_tbl %>%
    dplyr::filter(group %in% primary_lineages) %>%
    dplyr::filter(!is.na(cluster_medoid), cluster_medoid != "-") %>%
    dplyr::distinct(protein, cluster_medoid, group) %>%
    dplyr::group_by(protein, cluster_medoid) %>%
    dplyr::summarise(
      groups_present = list(sort(unique(group))),
      classification = classify_groups(groups_present[[1]]),
      .groups        = "drop"
    )
}

presence_map <- make_presence_map(ann_tbl)

# Attach cross-group classification to individual group summaries
indiv_summary <- indiv_summary %>%
  dplyr::left_join(presence_map, by = c("protein", "cluster_medoid")) %>%
  dplyr::relocate(classification, .after = hotspot_type) %>%
  dplyr::mutate(
    groups_present = ifelse(
      lengths(groups_present) == 0,
      "-",
      vapply(groups_present, paste, character(1L), collapse = ",")
    )
  )

readr::write_csv(
  indiv_summary,
  file.path(out_dir_base, paste0("annotations-individual-groups-", mode_suffix, ".csv")),
  na = "-"
)

cat("\nSaved cluster annotations by group with medoid classification:\n  - ",
    paste0("annotations-individual-groups-", mode_suffix, ".csv"), "\n", sep = "")

## ===== Cluster × host × radii matrix + per-host driver mean frequencies ======

# Choose canonical driver per (protein, medoid) based on most common driver across hosts,
# with a tie-break by overall mean frequency.
choose_driver_per_medoid <- function(ann_tbl, freqs_ps) {
  cand <- ann_tbl %>%
    dplyr::filter(group %in% primary_lineages,
                  !is.na(cluster_medoid), cluster_medoid != "-") %>%
    dplyr::transmute(
      protein,
      cluster_medoid,
      driver_site = dplyr::na_if(driver_site, "-")
    )

  mode_tbl <- cand %>%
    dplyr::filter(!is.na(driver_site)) %>%
    dplyr::count(protein, cluster_medoid, driver_site, name = "n") %>%
    dplyr::group_by(protein, cluster_medoid) %>%
    dplyr::filter(n == max(n)) %>%
    dplyr::ungroup()

  if (!nrow(mode_tbl)) {
    return(
      ann_tbl %>%
        dplyr::distinct(protein, cluster_medoid) %>%
        dplyr::mutate(driver_for_medoid = cluster_medoid)
    )
  }

  mean_by_driver <- freqs_ps %>%
    dplyr::filter(Group %in% primary_lineages) %>%
    dplyr::group_by(site) %>%
    dplyr::summarise(mean_all = mean(freq, na.rm = TRUE), .groups = "drop")

  driver_choice <- mode_tbl %>%
    dplyr::left_join(mean_by_driver, by = c("driver_site" = "site")) %>%
    dplyr::group_by(protein, cluster_medoid) %>%
    dplyr::slice_max(dplyr::coalesce(mean_all, -Inf), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      protein,
      cluster_medoid,
      driver_for_medoid = dplyr::coalesce(driver_site, cluster_medoid)
    )

  all_medoids <- ann_tbl %>% dplyr::distinct(protein, cluster_medoid)

  all_medoids %>%
    dplyr::left_join(driver_choice, by = c("protein", "cluster_medoid")) %>%
    dplyr::mutate(driver_for_medoid = dplyr::coalesce(driver_for_medoid, cluster_medoid))
}

# Mean frequency per host at a given site, with low-coverage exclusions
per_host_mean_at_site <- function(site, freqs_ps, removed_samples, protein_label) {
  indiv_groups <- primary_lineages
  prot_of_site <- sub("\\..*$", "", site)

  df <- freqs_ps %>%
    dplyr::filter(
      site == !!site,
      if (prot_of_site %in% c("PB1", "PB2", "PA"))
        Protein %in% c("PB1", "PB2", "PA")
      else
        Protein == !!prot_of_site,
      Group %in% indiv_groups
    ) %>%
    dplyr::mutate(SampleID = as.character(SampleID))

  if (!nrow(df)) {
    return(tibble::tibble(
      `freq-BALBF` = 0,
      `freq-BALBM` = 0,
      `freq-BL6F`  = 0,
      `freq-BL6M`  = 0
    ))
  }

  if (nrow(removed_samples)) {
    rmv <- removed_samples %>%
      dplyr::filter(protnm == protein_label) %>%
      dplyr::pull(SampleID) %>%
      as.character() %>%
      unique()
    if (length(rmv)) {
      df <- df %>% dplyr::filter(!(SampleID %in% rmv))
    }
  }

  agg <- df %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(mean_freq = mean(freq, na.rm = TRUE), .groups = "drop")

  tibble::tibble(
    `freq-BALBF` = dplyr::coalesce(agg$mean_freq[agg$Group == "BALBF"], 0),
    `freq-BALBM` = dplyr::coalesce(agg$mean_freq[agg$Group == "BALBM"], 0),
    `freq-BL6F`  = dplyr::coalesce(agg$mean_freq[agg$Group == "BL6F"],  0),
    `freq-BL6M`  = dplyr::coalesce(agg$mean_freq[agg$Group == "BL6M"],  0)
  )
}

# Construct cluster matrix with per-host radii, driver frequencies, and hotspot types
make_cluster_matrix_with_freqs <- function(annot_tbl, presence_tbl, freqs_ps, removed_samples) {

  driver_map <- choose_driver_per_medoid(annot_tbl, freqs_ps)
  indiv_groups <- primary_lineages   # BALBF BALBM BL6F BL6M

  core <- annot_tbl %>%
    dplyr::filter(group %in% indiv_groups,
                  !is.na(cluster_medoid),
                  cluster_medoid != "-") %>%
    dplyr::mutate(cdist = as.numeric(cdist)) %>%
    dplyr::distinct(protein, cluster_medoid, group, cdist)

  grp_r <- core %>%
    dplyr::group_by(protein, cluster_medoid, group) %>%
    dplyr::summarise(
      radii = paste(sort(unique(cdist)), collapse = ","),
      .groups = "drop"
    )

  wide <- tidyr::pivot_wider(
    grp_r,
    id_cols     = c(protein, cluster_medoid),
    names_from  = group,
    values_from = radii
  )

  # Guarantee presence columns for all four hosts
  for (g in indiv_groups) {
    if (!g %in% names(wide)) wide[[g]] <- NA_character_
  }

  wide <- wide %>%
    dplyr::mutate(
      BALBF = dplyr::coalesce(.data$BALBF, "-"),
      BALBM = dplyr::coalesce(.data$BALBM, "-"),
      BL6F  = dplyr::coalesce(.data$BL6F,  "-"),
      BL6M  = dplyr::coalesce(.data$BL6M,  "-")
    ) %>%
    dplyr::left_join(driver_map, by = c("protein", "cluster_medoid"))

  # Driver frequencies
  freq_cols <- purrr::pmap_dfr(
    list(wide$driver_for_medoid, wide$protein),
    function(driver_site, prot) {
      prot_key <- if (prot %in% c("PB1","PB2","PA")) "POL" else prot
      per_host_mean_at_site(driver_site, freqs_ps, removed_samples, prot_key)
    }
  )
  wide2 <- dplyr::bind_cols(wide, freq_cols)

# Hotspot type table (may miss columns for hosts with no hotspots)
htypes <- annot_tbl %>%
  dplyr::filter(
    group %in% indiv_groups,
    !is.na(cluster_medoid),
    cluster_medoid != "-"
  ) %>%
  dplyr::group_by(protein, cluster_medoid, group) %>%
  dplyr::summarise(
    hotspot_type = dplyr::first(na.omit(hotspot_type)),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    id_cols     = c(protein, cluster_medoid),
    names_from  = group,
    values_from = hotspot_type,
    names_prefix = "type-"
  )

# Ensure all type-* columns exist even when the host had no hotspots
for (g in indiv_groups) {
  cname <- paste0("type-", g)
  if (!cname %in% names(htypes)) {
    htypes[[cname]] <- NA_character_
  }
}

wide3 <- wide2 %>%
  dplyr::left_join(htypes, by = c("protein", "cluster_medoid")) %>%
  dplyr::left_join(presence_tbl, by = c("protein", "cluster_medoid")) %>%
  dplyr::select(
    protein, cluster_medoid,
    BALBF, BALBM, BL6F, BL6M,
    `freq-BALBF`, `freq-BALBM`, `freq-BL6F`, `freq-BL6M`,
    `type-BALBF`, `type-BALBM`, `type-BL6F`, `type-BL6M`,
    classification
  ) %>%
  dplyr::arrange(protein, cluster_medoid)

return(wide3) 

}

clusters_matrix <- make_cluster_matrix_with_freqs(
  ann_tbl,
  presence_map,
  freqs_per_sample,
  removed_samples
)

readr::write_csv(
  clusters_matrix,
  file.path(out_dir_base, paste0("clusters-matrix-", mode_suffix, ".csv")),
  na = "-"
)

cat("\nSaved cluster matrix with per-host mean driver-site frequencies and hotspot types:\n  - ",
    paste0("clusters-matrix-", mode_suffix, ".csv"), "\n", sep = "")

## ================== XLSX with bolded per-host freqs and types =================
suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(readr)
  library(tidyr)
})

# --------- File paths ----------
in_csv   <- file.path(out_dir_base, paste0("clusters-matrix-", mode_suffix, ".csv"))
out_xlsx <- file.path(out_dir_base, paste0("clusters-matrix-", mode_suffix, "-significant.xlsx"))

# --------- Load matrix ----------
mat <- readr::read_csv(in_csv, show_col_types = FALSE)

# Required minimal columns for XLSX export
required <- c("protein", "cluster_medoid",
              "freq-BALBF", "freq-BALBM", "freq-BL6F", "freq-BL6M",
              "classification")

# Ensure minimal columns exist; missing → fill with defaults
for (col in required) {
  if (!col %in% names(mat)) {
    if (grepl("^freq-", col)) {
      mat[[col]] <- 0
    } else {
      mat[[col]] <- "-"
    }
  }
}

# --------- Build the combined hotspot-type-per-host column ----------
host_order <- c("BALBF", "BALBM", "BL6F", "BL6M")

# ensure type columns exist
for (g in host_order) {
  col <- paste0("type-", g)
  if (!col %in% names(mat)) {
    mat[[col]] <- NA_character_
  }
}

mat$hotspot_type_per_host <- apply(mat, 1, function(row){
  parts <- c()
  for (g in host_order) {
    tcol <- paste0("type-", g)
    val  <- row[[tcol]]
    if (!is.na(val) && val != "-") {
      parts <- c(parts, paste0(g, ": ", val))
    }
  }
  if (length(parts) == 0) return("-")
  paste(parts, collapse = "; ")
})

# --------- Remove unused columns for XLSX output ----------
out <- mat %>%
  select(
    protein, cluster_medoid,
    `freq-BALBF`, `freq-BALBM`, `freq-BL6F`, `freq-BL6M`,
    classification,
    hotspot_type_per_host
  )

# --------- Determine bolding based on radii presence ----------
# Presence columns may not exist; treat non-existing or "-" as FALSE
presence_cols <- c("BALBF", "BALBM", "BL6F", "BL6M")
for (pc in presence_cols) {
  if (!pc %in% names(mat)) mat[[pc]] <- "-"
}

has_numeric_value <- function(x) {
  x <- ifelse(is.na(x), "-", trimws(as.character(x)))
  grepl("\\d", x)
}

b_BALBF <- has_numeric_value(mat$BALBF)
b_BALBM <- has_numeric_value(mat$BALBM)
b_BL6F  <- has_numeric_value(mat$BL6F)
b_BL6M  <- has_numeric_value(mat$BL6M)

# --------- Write XLSX ----------
wb <- createWorkbook()
addWorksheet(wb, "matrix")

# Write data
writeData(wb, "matrix", out, headerStyle = createStyle(textDecoration = "bold"))

# numeric formatting
freq_cols <- which(startsWith(names(out), "freq-"))
num_style <- createStyle(numFmt = "0.000")
addStyle(
  wb, "matrix", num_style,
  rows = 2:(nrow(out) + 1), cols = freq_cols,
  gridExpand = TRUE, stack = TRUE
)

# bolding
bold_style <- createStyle(textDecoration = "bold")
r_idx <- 2:(nrow(out) + 1)

c_BALBF <- which(names(out) == "freq-BALBF")
c_BALBM <- which(names(out) == "freq-BALBM")
c_BL6F  <- which(names(out) == "freq-BL6F")
c_BL6M  <- which(names(out) == "freq-BL6M")

if (any(b_BALBF)) addStyle(wb, "matrix", bold_style, rows = r_idx[b_BALBF], cols = c_BALBF, stack = TRUE)
if (any(b_BALBM)) addStyle(wb, "matrix", bold_style, rows = r_idx[b_BALBM], cols = c_BALBM, stack = TRUE)
if (any(b_BL6F))  addStyle(wb, "matrix", bold_style, rows = r_idx[b_BL6F],  cols = c_BL6F,  stack = TRUE)
if (any(b_BL6M))  addStyle(wb, "matrix", bold_style, rows = r_idx[b_BL6M],  cols = c_BL6M,  stack = TRUE)

setColWidths(wb, "matrix", cols = 1:ncol(out), widths = "auto")
saveWorkbook(wb, out_xlsx, overwrite = TRUE)

cat("Wrote XLSX: ", out_xlsx, "\n")



