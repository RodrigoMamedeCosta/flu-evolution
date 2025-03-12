## ----setup, include=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(
	echo = FALSE,
	message = FALSE,
	warning = FALSE
)
library(runner)
library(tidyverse)



## ---------------------------------------------------------------------------------------
### 1) READ IN DISTANCE MATRICES ###
HAmat <- read.csv("Input Data/dm-HAtrimer.csv")
M1mat <- read.csv("Input Data/dm-M1.csv")
M2mat <- read.csv("Input Data/dm-M2.csv")
NPmat <- read.csv("Input Data/dm-NP.csv")
NRmat <- read.csv("Input Data/dm-NR.csv")
NS1mat <- read.csv("Input Data/dm-NS1dimer.csv")

pol.6rr7 <- read.csv("Input Data/dm-POL-6RR7.csv")
pol.6qnw <- read.csv("Input Data/dm-POL-6QNW.csv")

#Define polymerase to analyze.
#POLmat_list <- list("6RR7" = pol.6rr7) - deprecated, missing some C-terminus aminoacids.
POLmat_list <- list("6QNW" = pol.6qnw)

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

### EXAMPLE: run only "AllHosts" plus a custom group = "treatments"
groups <- c("AllHosts", hosts, hoststrains, sexes, treatments)

#groups <- c("AllHosts", hosts, hoststrains, sexes)

### 7) DETECT DATA TYPE ###
data_type <- ifelse(grepl("ALLDATA", csv_file_path, ignore.case = TRUE), "All_Data", "No_Indels")

### 9) MAIN ANALYSIS LOOPS ###
for (pol.type in names(POLmat_list)) {
  
  cat("\n=============================\n")
  cat("Processing Polymerase Type:", pol.type, "\n")
  cat("=============================\n\n")
  
  
  
  # Initialize/reset data frames for each polymerase type
  results_df <- data.frame(
    group    = character(),
    quantile = numeric(),
    cdist    = numeric(),
    site     = character(),
    pol.dm   = character(),
    stringsAsFactors = FALSE
  )
  
  removed_summary <- data.frame(
    pol_type = character(),
    cdist    = numeric(),
    group    = character(),
    protnm   = character(),
    SampleID = character(),
    stringsAsFactors = FALSE
  )
  
  all_results <- list()  # Reset for each pol.type
  
  POLmat <- POLmat_list[[pol.type]]
  
  for (cdist in cdist_values) {
    for (group in groups) {
      
      # (A) SUBSET THE DATA
      if (group == "AllHosts") {
        counts <- allcounts
      } else if (group %in% allcounts$Host) {
        counts <- allcounts %>% filter(Host == group)
      } else if (group %in% hoststrains) {
        counts <- allcounts %>% filter(HostStrain == group)
      } else if (group %in% sexes) {
        counts <- allcounts %>% filter(Sex == group)
      } else if (group %in% treatments) {
        counts <- allcounts %>% filter(Treatment == group)
      } else {
        cat("Warning: Group '", group, "' not found in columns. Skipping.\n")
        next
      }
      
      # Define a function to find close residues
      findclosespace <- function(v) spacemat$Residue[v < cdist]
      
      #Create spacedat and crspacedat that will have the concatenated spatial data.
      cspacedat <- NULL
      crspacedat <- NULL
      
      #Find hotspots and do randomization
      # Loop over every protein
      for (protnm in protnms) {
        spacemat <- get(paste0(protnm,"mat"))
      
      cat("\n","Analyzing",protnm, "in", group, "\n\n")
        
      # Extract unique chains from spacemat and alphabetically order them
      unique_chains <- sort(unique(sub("\\..*$", "", spacemat$Residue)))
      
      # Check if spacemat string format has 4 sections (3 period separators)
      if (length(unlist(strsplit(spacemat$Residue[1], "\\."))) == 4) {
        # Initialize an empty dataframe for the results
        counts_chain <- data.frame(matrix(NA, nrow = nrow(counts), ncol = 0))
      
        # Include the first 9 columns directly
        counts_chain[, 1:9] <- counts[, 1:9]
      
        # Iterate through each column in counts starting from the 10th column
        for (col in names(counts)[10:length(counts)]) {
          # Extract the common part of the column name (e.g., NR.748)
          common_part <- sub("\\.[^.]*\\.[^.]*$", "", col)  # Updated to extract only protein and amino acid number
      
          # Check if the common part matches any residue in spacemat
          matching_residues <- grep(common_part, spacemat$Residue)
          
          if (length(matching_residues) > 0) {
            # Create new columns for each chain
            for (chain in unique_chains) {
              new_col_name <- paste0(chain, ".", col)
              counts_chain[[new_col_name]] <- counts[[col]]
            }
          }
        }
      cat(" Data formatting successful\n")
      } else {
         # If spacemat column names have only 3 substrings (no chain), save original counts file as counts_chain and print a message
         counts_chain <- counts
        cat("WARNING: NOT HOMOPOLYMER OR FORMAT INCORRECT\n")
      }
      
      #Sort names alphabetically to prevent problems downstream
      # Get names from the 10th column onwards
      relevant_colnames <- names(counts_chain)[10:ncol(counts_chain)]
      
      # Sort the column names alphabetically
      sorted_colnames <- sort(relevant_colnames)
      
      # Concatenate the first 9 column names with the sorted column names
      all_colnames <- c(names(counts_chain)[1:9], sorted_colnames)
      
      # Subset counts_chain based on the concatenated column names
      counts_chain <- counts_chain[, all_colnames]
      
      
      # Define the function to convert format
      convert_format <- function(string) {
          sub("\\.([^.]+)\\.(Anc|Der)$", ".\\2", string)
      }
      
      # Apply the function to all column names of counts_chain
      colnames(counts_chain) <- lapply(colnames(counts_chain), convert_format)
      
      #Show format
      #cat("\nCurrent counts_chain format:", names(counts_chain)[10:13], "\n")
      
      # Shorten the names FOR BOTH ROWS AND COLUMNS TO REMOVE AMINOACID NAME FROM SPACEMAT
          #shorten row names - this operation removes anything after and including the last period.
          spacemat$Residue <- sub("\\.[^.]*$", "", spacemat$Residue)
      
          # Get the current column names
          current_colnames <- colnames(spacemat)
          # Remove everything after the last period in each column name
          new_colnames <- sub("\\.[^.]*$", "", current_colnames)
          # Assign the modified column names back to the data frame
          colnames(spacemat) <- new_colnames
      
      # Find close residues for the third column of spacemat, to check if it's working
      findclosespace(spacemat[,40])
      #NOTE - THIS IS NOW ALSO FINDING AMINOACIDS IN OTHER CHAINS, WHICH IS GOOD, FOR EXAMPLE WHEN WE LOOK AT COLUMN 130
      
      # Convert spacemat columns to a list and find close residues for each
      spacelist <- as.list(spacemat[,2:ncol(spacemat)])
      spacedist <- lapply(spacelist, findclosespace)
      
      # Plot the distribution of the number of close residues for each site
      #plot(unlist(lapply(spacedist, length)), type="b")
      
      #Remove sample data (first 9 columns) from counts_only
      counts_only <- counts_chain[, -c(1:9)]
      
      # Create a vector of column names with everything before the third period in the counts dataframe
      # i.e. remove the Anc/Der
      colnames_noancder <- sapply(strsplit(names(counts_only), "\\."), function(x) paste(x[1:3], collapse = "."))
      
      
      # Subset the counts data frame based on the matching prefixes
      tmp <- counts_only[, colnames_noancder %in% spacemat$Residue]
      
      
      # Separate ancestral and derived
      spaceAnccounts <- tmp[grep("\\.Anc$", names(tmp))]
      spaceDercounts <- tmp[grep("\\.Der$", names(tmp))]
      
      # Updating the sub patterns for column renaming
      names(spaceAnccounts) <- sub("^([^\\.]+\\.[^\\.]+\\.[^\\.]+)\\..*$", "\\1", names(spaceAnccounts))
      names(spaceDercounts) <- sub("^([^\\.]+\\.[^\\.]+\\.[^\\.]+)\\..*$", "\\1", names(spaceDercounts))
      
      #Reorder columns
      spaceAnccounts <- spaceAnccounts[, names(spaceAnccounts) %in% names(spacedist)]
      spaceDercounts <- spaceDercounts[, names(spaceDercounts) %in% names(spacedist)]
      
      # Local copies, so low-coverage removal doesn't affect other proteins
      counts_chain_local   <- counts_chain
      spaceAnccounts_local <- spaceAnccounts
      spaceDercounts_local <- spaceDercounts
      
      # Calculate the **median** coverage per row (sample) for this protein
      # Using na.rm = TRUE just in case there are any NA cells
      sample_coverage <- apply(spaceAnccounts_local + spaceDercounts_local, 1, median, na.rm = TRUE)
      
      # Identify which samples are below the threshold
      lowcov_idx <- which(sample_coverage < covmin)
      
      if (length(lowcov_idx) > 0) {
        # Assuming counts_chain_local has a column named SampleID in the first column:
        removed_sample_ids <- counts_chain_local[[1]][lowcov_idx]
        cat("Removed samples", paste(removed_sample_ids, collapse = ", "),
            "due to median coverage below", covmin, "in", protnm, "\n")
        
        #  Build a dataframe to store this info
        removed_df <- data.frame(
          pol_type = pol.type,
          cdist    = cdist,
          group    = group,
          protnm   = protnm,
          SampleID = removed_sample_ids,
          stringsAsFactors = FALSE
        )
      
        # Append to the cumulative removed_summary
        removed_summary <- rbind(removed_summary, removed_df)
        
      } else {
        cat("No samples removed in", protnm, "for median coverage <", covmin, "\n")
      }
      
      # Keep only rows (samples) that pass coverage
      keep_idx <- (sample_coverage >= covmin)
      counts_chain_local   <- counts_chain_local[keep_idx, ]
      spaceAnccounts_local <- spaceAnccounts_local[keep_idx, ]
      spaceDercounts_local <- spaceDercounts_local[keep_idx, ]
      
      # Update the local sample count
      m_local <- nrow(counts_chain_local)
      
      # ALSO subset counts_only with the same keep_idx!
      counts_only <- counts_only[keep_idx, ]
      
      # Overwrite the originals so the existing code below remains unchanged
      counts_chain   <- counts_chain_local
      spaceAnccounts <- spaceAnccounts_local
      spaceDercounts <- spaceDercounts_local
      m              <- m_local
      
      # Apply space function along rows using spacedist
      sspaceAnccounts <- 0 * spaceAnccounts
      sspaceDercounts <- 0 * spaceDercounts
      
      # Identify common columns in spaceAnccounts and spaceDercounts
      common_cols <- intersect(names(spaceAnccounts), names(spaceDercounts))
      
      # Filter spacedist to include only indices that reference common columns
      filtered_spacedist <- lapply(spacedist, function(cols) cols[cols %in% common_cols])
      
      # Loop through each residue and apply the space function along rows
      for (inm in 1:length(filtered_spacedist)) {
          # Only proceed if filtered_spacedist at index inm is not empty
          if (length(filtered_spacedist[[inm]]) > 0) {
              tmp <- as.data.frame(spaceAnccounts[, filtered_spacedist[[inm]]])
              sspaceAnccounts[, inm] <- apply(tmp, 1, sum)
              tmp <- as.data.frame(spaceDercounts[, filtered_spacedist[[inm]]])
              sspaceDercounts[, inm] <- apply(tmp, 1, sum)
          }
      }
      
      # Calculate statistics for each residue
      spacedat <- data.frame(site=names(spaceAnccounts), stat1=0, stat2=0, stat3=0)
      
      # Calculate the number of rows in the dataset and the square root of the reciprocal.
      #In this case, we don't calculate the number of rows, we just state it, but if it's different from the number of samples in our counts dataframe it will not work.
      
      #Define function.
      lfrac <- 1/sqrt(m)
      
      # Loop through each residue's data and calculate summary statistics
      for (inm in 1:ncol(spaceAnccounts)) {
        nm <- names(spaceAnccounts)[inm]
        tmp <- data.frame(a=unlist(sspaceAnccounts[,inm]),
                          d=unlist(sspaceDercounts[,inm]))
        tmp$p <- tmp$d/(tmp$a + tmp$d)
        tmp$p[is.nan(tmp$p)] <- 0
        tmp <- tmp[order(tmp$p),]
        tmp <- cbind(s=m:1,tmp)
        
        # Calculate statistics and store them in the spacedat data frame
        spacedat$stat1[inm] <- mean(tmp$p)
        spacedat$stat2[inm] <- sum(tmp$s*tmp$p)/sum(tmp$s)
        newdat <- data.frame(s=lfrac*m)
        stat3 <- 0
        
        # Perform logistic regression and predict probabilities
        if (tmp$p[m-1] > 0) {
          lg <- glm(cbind(d,a) ~ s,tmp,family=binomial)
          stat3 <- predict(lg,newdata=newdat)
          spacedat$stat3[inm] <- exp(stat3)/(1 + exp(stat3))
        }
      }
      
      #Extract data for the protein
      Protdat <- spacedat
      
      # Set up for the randomization, names fixed for variable number of characters.
      spacenms <- Protdat$site
      
      #Ignore crap, take names only
      relevant_colnames <- names(counts_only)
      
      # Sort the column names alphabetically again. Probably not required but just in case
      sorted_colnames <- sort(relevant_colnames)
      
      # Subset counts_chain based on the sorted column names
      justspace <- counts_only[, sorted_colnames]
      
      # Update column names in justspace
      new_names <- paste0(rep(spacenms, each = 2), c(".Anc", ".Der"))
      names(justspace) <- new_names
      
      # Identify columns with NA names (introduced by the jury rigged way to populate columns at the beginning)
      na_column_indices <- is.na(names(justspace))
      
      # Exclude columns with NA names
      justspace <- justspace[, !na_column_indices]
      
      # justspace_cleaned now contains the DataFrame without the columns that had NA names
      
      
      #If needed, confirm integrity of dataframe by accessing csv
      #path_out = '.\\Outputs\\'
      #write.csv(justspace,paste(path_out,"justspacealpha.csv"),row.names=FALSE)
      

      # Initialize an empty data frame to store results from all repetitions
      allrspacedat <- NULL
      
      cat("\n","Beginning randomization for",protnm, "in", group, "\n\n")
      
      
      for (irep in 1:reps) {
        randspace <- justspace
        
        for (i in 1:nrow(justspace)) {
          # Shuffle the residue names
          rnms <- sample(spacenms)
      
      # Append .Anc and .Der to the shuffled names
          rnms_with_suffixes <- paste0(rep(rnms,each=2), rep(c(".Anc", ".Der"),length(rnms)))
      
          # Replace the row with shuffled values, ensuring the columns exist
          if(all(rnms_with_suffixes %in% names(justspace))) {
            randspace[i,] <- justspace[i, rnms_with_suffixes]
          }
        }
        # Load the "randspaceout.R" script   #########################################SPACEOUT 
        # Extract Ancestral and Derived counts for selected residues
        tmp <- randspace[, sapply(strsplit(names(randspace), "\\."), function(x) paste(x[1:3], collapse = ".")) %in% spacenms]
      
        spaceAnccounts <- tmp[, grep("Anc$", names(tmp))]
        spaceDercounts <- tmp[, grep("Der$", names(tmp))]
      
        # Set the names of Ancestral and Derived counts to match residue names
        names(spaceAnccounts) <- names(spaceDercounts) <- spacenms
      
        sspaceAnccounts <- matrix(0, nrow = nrow(spaceAnccounts), ncol = length(filtered_spacedist))
        sspaceDercounts <- matrix(0, nrow = nrow(spaceDercounts), ncol = length(filtered_spacedist))
      
        for (inm in 1:length(filtered_spacedist)) {
          if (length(filtered_spacedist[[inm]]) > 0) {
            if (!is.null(spaceAnccounts[, filtered_spacedist[[inm]], drop = FALSE]) && ncol(spaceAnccounts[, filtered_spacedist[[inm]], drop = FALSE]) > 0) {
              tmp <- spaceAnccounts[, filtered_spacedist[[inm]], drop = FALSE]
              sspaceAnccounts[, inm] <- apply(tmp, 1, sum)
            } else {
              cat("No sites found for Ancestral at index", inm, "- skipping\n")
            }
      
            if (!is.null(spaceDercounts[, filtered_spacedist[[inm]], drop = FALSE]) && ncol(spaceDercounts[, filtered_spacedist[[inm]], drop = FALSE]) > 0) {
              tmp <- spaceDercounts[, filtered_spacedist[[inm]], drop = FALSE]
              sspaceDercounts[, inm] <- apply(tmp, 1, sum)
            } else {
              cat("No sites found for Derived at index", inm, "- skipping\n")
            }
          }
        }
      
      # Calculate statistics for the randomized data
      rspacedat <- data.frame(site = names(spaceAnccounts),
                              stat1 = 0, stat2 = 0, stat3 = 0)
      
      for (inm in 1:ncol(spaceAnccounts)) {
        nm <- names(spaceAnccounts)[inm]
        tmp <- data.frame(a = unlist(sspaceAnccounts[, inm]),
                          d = unlist(sspaceDercounts[, inm]))
        
        # Calculate the proportion of Derived mutations
        tmp$p <- tmp$d / (tmp$a + tmp$d)
        tmp$p[is.nan(tmp$p)] <- 0
        
        # Sort data by mutation proportion and add a column for sequence position
        tmp <- tmp[order(tmp$p),]
        tmp <- cbind(s = m:1, tmp)
        
        # Calculate the three statistics for each residue
        rspacedat$stat1[inm] <- mean(tmp$p)  # Mean mutation proportion
        rspacedat$stat2[inm] <- sum(tmp$s * tmp$p) / sum(tmp$s)  # Weighted sum of positions
        newdat <- data.frame(s = lfrac * m)
        stat3 <- 0
        
        # If mutation proportion is not zero, fit a logistic regression model
        if (tmp$p[m - 1] > 0) {
          lg <- glm(cbind(d, a) ~ s, tmp, family = binomial)
          stat3 <- predict(lg, newdata = newdat)
          rspacedat$stat3[inm] <- exp(stat3) / (1 + exp(stat3))
        }
      }
       # End the "randspaceout.R" script
        
        # Print the current iteration
        cat("Iteration", irep, "\n")
      
        # Filter out rows with non-zero statistics
        rspacedat <- subset(rspacedat, stat1 + stat2 + stat3 > 0)
        
        # Append the results from the current iteration to the overall results
        allrspacedat <- rbind(allrspacedat, rspacedat)
      }
      
      # Print the current iteration
      cat("\n","Randomization complete for ",protnm,"\n\n")
      
      cspacedat <- rbind(cspacedat,spacedat)
      crspacedat <- rbind(crspacedat,allrspacedat)
      
      }
      
      # Define the full directory path including the data type and the 'space_outputs' subdirectory
      full_output_path <- paste0("Outputs/Raw_Outputs/",pol.type,"/","/Spatial_outputs/")
      
      # Create the directory if it doesn't exist, recursively creating all necessary subdirectories
      if (!dir.exists(full_output_path)) {
          dir.create(full_output_path, recursive = TRUE)
      }
      
      # Construct the complete output file path for spacedat
      output_file_path_spacedat <- paste0(full_output_path, group, "_", cdist, "A_", qtile, "_spacedat.csv")
      
      # Write the spacedat data to the CSV file
      write.csv(cspacedat, output_file_path_spacedat, row.names = FALSE)
      
      # Construct the complete output file path for randomspacedat
      output_file_path_randomspacedat <- paste0(full_output_path, group, "_", cdist, "A_", qtile, "_randomspacedat.csv")
      
      # Write the randomspacedat data to the CSV file
      write.csv(crspacedat, output_file_path_randomspacedat, row.names = FALSE)
      
      
      #End message
      cat(paste("\n\n",
                reps, "randomizations successfully performed for",group,"\n",
                "Sphere radius = ", cdist, "Å\n",
                "Proceed with analysis"))
      
      
      ###### Final analysis and save data #####
      #qtile <- 0.999 # Set the quantile threshold if different from the parameters section
      
      #Test real data against randomized data
      Protspace1 <- cspacedat$site[cspacedat$stat1 > quantile(crspacedat$stat1, qtile)]
      Protspace2 <- cspacedat$site[cspacedat$stat2 > quantile(crspacedat$stat2, qtile)]
      Protspace3 <- cspacedat$site[cspacedat$stat3 > quantile(crspacedat$stat3, qtile)]
      
      # Print the sites in each Protspace subset
      #print(Protspace1)
      #print(Protspace2)
      #print(Protspace3)
      
      ####Clean up data####
      
      #Remove chain identifiers by removing the first substring and separator for each element
      Protspace1 <- sub("^[^.]*\\.", "", Protspace1)
      Protspace2 <- sub("^[^.]*\\.", "", Protspace2)
      Protspace3 <- sub("^[^.]*\\.", "", Protspace3)
      
      #Remove duplicate values
      Protspace1 <- unique(Protspace1)
      Protspace2 <- unique(Protspace2)
      Protspace3 <- unique(Protspace3)
      
      #This part now will correct the count of the HA aminoacids, to account for the cleaving of the signal peptide(16 aminoacids). Other proteins are correct.
      # Create a list containing the three Protspace vectors
      Protspace_list <- list(Protspace1 = Protspace1, Protspace2 = Protspace2, Protspace3 = Protspace3)
      
      # Initialize an empty list to store adjusted vectors
      Protspace_adjusted_list <- list()
      
      # Loop through each Protspace vector in the list
      for (i in seq_along(Protspace_list)) {
        # Convert the current vector into a data frame
        Protspace_df <- data.frame(
          Prefix = sapply(Protspace_list[[i]], function(x) strsplit(x, "\\.")[[1]][1]),
          Number = sapply(Protspace_list[[i]], function(x) strsplit(x, "\\.")[[1]][2]),
          stringsAsFactors = FALSE
        )
        
        # Ensure that Number is treated as numeric, handling NAs properly
        Protspace_df$Number <- as.numeric(Protspace_df$Number)
        
        # Apply the adjustment conditionally
        Protspace_df$AdjustedNumber <- ifelse(Protspace_df$Prefix == "HA", Protspace_df$Number - 16, Protspace_df$Number)
        
        # Reassemble the adjusted identifiers, ensuring no NAs are introduced
        # Directly reference dataframe columns to avoid scope issues
        Protspace_adjusted <- paste(Protspace_df$Prefix, sprintf("%03d", ifelse(is.na(Protspace_df$AdjustedNumber), Protspace_df$Number, Protspace_df$AdjustedNumber)), sep = ".")
        
        # Store the adjusted vector in the new list
        Protspace_adjusted_list[[names(Protspace_list)[i]]] <- Protspace_adjusted
      }
      # Update the original Protspace vectors
      Protspace1 <- Protspace_adjusted_list$Protspace1
      Protspace2 <- Protspace_adjusted_list$Protspace2
      Protspace3 <- Protspace_adjusted_list$Protspace3
      
      
      ########SELECT STATISTIC FOR REPORTS HERE######
      
      ###Make Report based on preferred statistic####
      #Define the results we're interested in. Stat 2 was the weighted sum of positions.
      results <- Protspace2
      
      # Finally, store it in the list:

      if (!exists("results")) {
        # If your script does not define "results" for some reason:
        results <- character(0)
      }

      # (C) STORE RESULTS IN NAMED LIST
      # Make sure you do this AFTER results is created in your script:
      var_name <- paste0("results_", group)  # optional
      all_results[[group]] <- results

      # (D) ALSO APPEND INTO results_df
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

      # Clean up for next iteration
      rm(results)
      
      cat("\nSaved data for group:", group,
          "for quantile =", qtile,
          " with sphere radius =", cdist,
          " minimum median coverage =", covmin,
          "and polymerase type =", pol.type, "\n\n")
    }
    # End of cdist loop
  }


### 10) AFTER ALL GROUPS/POLYMERS are processed, do final intersections 

if ("AllHosts" %in% names(all_results)) {
  # Make a directory to hold intersections
  if (!dir.exists("Outputs/Intersect_Outputs")) dir.create("Outputs/Intersect_Outputs")
  
  for (g in names(all_results)) {
    if (g == "AllHosts") next  # skip self

    intersect_with_all <- intersect(all_results[[g]], all_results[["AllHosts"]])
    
    
    #name the csv
    output_file <- sprintf("Outputs/Intersect_Outputs/%s_vs_AllHosts_%s.csv", g, pol.type)
    
    # Write the CSV
    write.csv(intersect_with_all, output_file, row.names = FALSE, quote = TRUE)
    
    cat("Intersected", g, "with AllHosts. Found", length(intersect_with_all), "common sites.\n")
  }
}

### 11) OPTIONAL: If you want to do “shared by multiple groups” logic,
### you can do that by combining intersections or using combn(...) 
### similarly to how you were doing it. For example:

# Example: find sites shared by *all* groups (besides "AllHosts"):
other_groups <- setdiff(names(all_results), "AllHosts")
if (length(other_groups) > 1) {
  shared_all <- Reduce(intersect, all_results[other_groups])
  cat("\nSites shared by all groups:\n", paste(shared_all, collapse=", "), "\n")
}

### 12) WRITE OUT THE results df
#Add polymerase column
results_df <- results_df %>% mutate(pol.dm = pol.type)
 
# Define the output path for the CSV file
  csv_output_path <- paste0("Outputs/table_output_", pol.type, ".csv")
  
  # Write the results_df data frame to a CSV file
  write.csv(results_df, csv_output_path, row.names = FALSE)
  
### 13) WRITE OUT REMOVED SUMMARY ###
  write.csv(removed_summary, paste0("Outputs/removed_samples_", pol.type, ".csv"), row.names = FALSE)
  
  cat("\nAnalysis complete for polymerase type:", pol.type, "\n")



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
                           paste0(cat_name, "_WithinGroup_Report_", pol.type,".txt"))
  
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

  # End of pol.type loop
}

cat("\nRolling Sphere Complete.\n\n")


## ---------------------------------------------------------------------------------------
# Group names
hosts        <- c("BALBF", "BALBM", "BL6F", "BL6M")
hoststrains  <- c("BALBC", "C57BL6")
sexes        <- c("FEMALE", "MALE")

# Define all groups including "AllHosts"
groups <- c("AllHosts", hosts, hoststrains, sexes)
cdists <- c(6, 12)  # The two distance values

# Generate all expected column names (_6 and _12 for each group)
expected_cols <- c("site", paste0(rep(groups, each = length(cdists)), "_", rep(cdists, times = length(groups))), "radius")

# Read and combine CSV files
outputs_dir <- "Outputs/"
table_output_files <- list.files(path = outputs_dir, pattern = "^table_output_.*\\.csv$", full.names = TRUE)

if (length(table_output_files) == 0) stop("No 'table_output_pol.type.csv' files found in the Outputs directory.")

combined_results <- bind_rows(lapply(table_output_files, read_csv))

# Keep only desired cdists and polymerase type 6QNW
target_pol <- "6QNW"
filtered_results <- combined_results %>%
  filter(cdist %in% cdists & pol.dm == target_pol) %>%
  distinct(site, group, cdist, .keep_all = TRUE)

# Create presence matrix
presence_df <- filtered_results %>%
  mutate(presence = "SIGNIF") %>%
  select(site, group, cdist, presence) %>%
  pivot_wider(
    names_from = c(group, cdist),
    values_from = presence,
    values_fill = list(presence = "-")
  )

# Add 'radius' column
radius_df <- filtered_results %>%
  group_by(site) %>%
  summarize(radius = paste(sort(unique(cdist)), collapse = ", "), .groups = "drop")

presence_df <- presence_df %>%
  left_join(radius_df, by = "site")

# Ensure all expected columns exist
missing_cols <- setdiff(expected_cols, names(presence_df))
for (col in missing_cols) presence_df[[col]] <- "-"

# Reorder columns
presence_df <- presence_df %>%
  select(all_of(expected_cols))

# Save the output
write.csv(presence_df, "Outputs/combined.rolling.output.table.csv", row.names = FALSE)

cat("\nCombined rolling sphere output table saved.\n")




## ---------------------------------------------------------------------------------------
# Use the already-created presence dataframe
df <- presence_df

# Function to find amino acids within a given distance
find_within_distance <- function(distance_matrices, site, distance) {
  protein <- sub("\\..*", "", site)
  position <- sub(".*\\.(\\d+)", "\\1", site)
  search_id <- paste0(protein, ".", position)

  results <- character(0)

  get_neighbors <- function(matrix, search_id, distance) {
    matching_row <- which(grepl(search_id, matrix$Residue))[1]
    if (!is.na(matching_row)) {
      distances <- matrix[matching_row, -1]
      within_distance <- distances <= distance
      sites_within_distance <- colnames(distances)[within_distance]
      return(unique(sites_within_distance[!is.na(sites_within_distance)]))
    }
    return(character(0))
  }

  if (protein %in% c("PB1", "PB2", "PA")) {
    if ("POL" %in% names(distance_matrices)) {
      results <- c(results, get_neighbors(distance_matrices[["POL"]], search_id, distance))
    }
  } else {
    for (matrix_name in names(distance_matrices)) {
      result <- get_neighbors(distance_matrices[[matrix_name]], search_id, distance)
      if (length(result) > 0) {
        results <- c(results, result)
        break
      }
    }
  }

  return(paste(unique(results), collapse = ", "))
}

# Distance matrices list (only 6QNW polymerase)
distance_matrices <- list(
  HA = HAmat,
  NS1 = NS1mat,
  M1 = M1mat,
  M2 = M2mat,
  NR = NRmat,
  NP = NPmat,
  POL = pol.6qnw  # Ensure only 6QNW is used
)

# Define cdist distances
cdist <- c(6, 12)

# Add distance-based columns
df <- df %>%
  rowwise() %>%
  mutate(
    Sites_within_6 = find_within_distance(distance_matrices, site, cdist[1]),
    Sites_within_12 = find_within_distance(distance_matrices, site, cdist[2])
  ) %>%
  ungroup()


# Read the frequency file
freqs <- read.csv("Input Data/Mean nonsyn variant frequency by strain.csv")

# Function to extract frequency at the exact df$site
get_freq_at_site <- function(site, group, freqs) {
  if (is.na(site) || site == "") return(NA)

  # Extract Protein and Position from df$site (e.g., HA.167 → Protein = HA, AA.pos = 167)
  parsed_site <- strsplit(site, "\\.")[[1]]
  if (length(parsed_site) < 2) return(NA)  # Handle parsing errors
  
  protein <- parsed_site[1]
  position <- as.numeric(parsed_site[2])

  # Retrieve the frequency for the given group at the exact site
  freq_value <- freqs %>%
    filter(Protein == protein, AA.pos == position) %>%
    pull(group)

  if (length(freq_value) == 0) return(NA) else return(round(freq_value, 4))
}

# Add frequency retrieval columns for each group
for (group in groups) {
  col_name <- paste0("Freq_", group)
  df[[col_name]] <- sapply(df$site, get_freq_at_site, group = group, freqs = freqs)
}

# Save the results to a new CSV file
write.csv(df, "Outputs/combined.rolling.output.table.with.driven.sites.by.mean.freq.csv", row.names = FALSE)

cat("\nCombined table with hotspot Nexus sites saved.\n")





## ---------------------------------------------------------------------------------------
# Load necessary libraries
library(openxlsx)

# Define output file
output_file <- "Outputs/processed_data.xlsx"

# Create a new workbook
wb <- createWorkbook()
addWorksheet(wb, "Sheet1")

# Write the dataframe
writeData(wb, "Sheet1", df)

# Save the workbook
saveWorkbook(wb, output_file, overwrite = TRUE)

cat("\nExcel file saved. Apply conditional formatting manually in Excel.\n")



## ---------------------------------------------------------------------------------------
# Function to retrieve the frequency at the exact site
get_freq_at_site <- function(sites_str, group, freqs) {
  if (is.na(sites_str) || sites_str == "") return(NA)
  
  # Extract protein and position from site strings
  site_list <- strsplit(sites_str, ", ")[[1]]
  parsed_sites <- do.call(rbind, strsplit(site_list, "\\."))

  if (ncol(parsed_sites) < 3) return(NA)  # Handle possible parsing errors
  
  proteins <- parsed_sites[,2]
  positions <- as.numeric(parsed_sites[,3])

  # Retrieve frequency for the given group at the first site in the list
  freq_value <- freqs %>%
    filter(Protein == proteins[1], AA.pos == positions[1]) %>%
    pull(group)

  if (length(freq_value) == 0) return(NA) else return(round(freq_value, 4))
}

# Add frequency retrieval columns for each group and cdist
for (group in groups) {
  for (cdist in cdists) {
    col_name <- paste0("Freq_", group, "_", cdist)
    sites_col <- paste0("Sites_within_", cdist)
    df[[col_name]] <- sapply(df[[sites_col]], get_freq_at_site, group = group, freqs = freqs)
  }
}


# Save the results to a new CSV file
write.csv(df, "Outputs/combined.rolling.output.table.with.driven.sites.by.mean.freq.csv", row.names = FALSE)

cat("\nCombined table with hotspot Nexus sites saved.\n")

