# DrugBank Dependencies
library(dbparser)
library(org.Hs.eg.db)
library(magrittr)
library(stringr)

# ClinicalTrials.gov Dependencies
library(ctrdata)
library(nodbi)
library(httr)

# PURPOSE ----
# 1. Retrieve drug name information from clinicaltrials.gov
# 2. Map clinicaltrials.gov names to drugbank identifiers. This code mainly attempts to clean up intervention drug names before mapping. 

# PARAMETERS ----
## Parameters for clinical trials search 
term <- "rhabdomyosarcoma" 
short_term <- term

## Parameters for DrugBank
# NOTE: Must be done once per session; this will take a while
load_data <- TRUE
# NOTE: Download Drugbank XML from: https://go.drugbank.com/releases/latest
drugbank_file <- "drugbank_all_full_database_5.1.8.xml"

# STEP 1: GET CLINICAL TRIAL DATA ----
term <- gsub(" ", "+", term)

# Download trial information
sqlite_file <- paste0(term, ".sql")

if(!file.exists(sqlite_file)) {
  db <- nodbi::src_sqlite(
    dbname = sqlite_file,
    collection = "trials"
  )
  
  url <- paste0("https://clinicaltrials.gov/ct2/results?cond=&term=", term, "&cntry=&state=&city=&dist=")
  ctrLoadQueryIntoDb(
    queryterm = url,
    con = db
  )  
} else {
  db <- nodbi::src_sqlite(
    dbname = sqlite_file,
    collection = "trials"
  )
}

# Convert to data.frame 
df <- nodbi::docdb_get(db, "trials")
tmp_results <- dbGetFieldsIntoDf(colnames(df), con = db, stopifnodata = FALSE)

idx <- lapply(tmp_results$condition, function(x) {
  any(grepl(short_term, tolower(x)))
}) %>% unlist
results <- tmp_results[idx,]

# Get trial drugs
drugs <- data.frame(nct_id=character(0), phase=character(0), overall_status=character(0), drug_name=character(0), stringsAsFactors=FALSE)

for(i in 1:nrow(results)) {
  #i <- 1
  cat("I: ", i, "\n")
  
  cur_entries <- results$intervention[[i]]
  
  if(!is.na(cur_entries) && is.null(nrow(cur_entries)) && (cur_entries$intervention_type == "Drug")) {
    tmp_drugs <- cur_entries$intervention_name
    
    if(length(tmp_drugs) > 0) {
      #drugs <- c(drugs, tmp_drugs)
      tmp <- data.frame(nct_id=results$`_id`[i], phase=results$phase[i], overall_status=results$overall_status[i], drug_name=tolower(tmp_drugs), stringsAsFactors=FALSE)
      drugs <- rbind(drugs, tmp)  
    }
  }
  
  if(!is.null(nrow(cur_entries))) {
    tmp_drugs <- cur_entries$intervention_name[cur_entries$intervention_type == "Drug"]
    
    if(length(tmp_drugs) > 0) {
      #drugs <- c(drugs, tmp_drugs)
      tmp <- data.frame(nct_id=results$`_id`[i], phase=results$phase[i], overall_status=results$overall_status[i], drug_name=tolower(tmp_drugs), stringsAsFactors=FALSE)
      drugs <- rbind(drugs, tmp)  
    }
  }  
}

# Summarize
drug_names <- drugs$drug_name %>% unique %>% sort(., decreasing=FALSE)

# STEP 2: EXTRACT DRUGBANK DATA ----
# NOTE: This can take a while
if(load_data) {
  tmp <- read_drugbank_xml_db(drugbank_file)
  drug_targets_actions <- targets_actions()
  drug_targets <- targets()
  #drug_target_links <- targets_links()
  drug_target_ids <- targets_polypep_ex_ident()
  
  #drug_ids <- drug_external_links()
  drug_ids <- drug_ex_identity()
  drug_synonyms <- drug_syn()
  drug_gen_info <- drug_general_information()
  #drug_brand <- drug_intern_brand()
  #drug_elem <- drug_element()  
}

# CONSTRUCT DRUG NAME-ID TABLE ----
tmp_drug_ids <- data.frame(db_id=character(0), drug_name=character(0), stringsAsFactors=FALSE)
# Wikipedia is the only resource in the table with human readable names 
t1 <- data.frame(db_id=drug_ids$parent_key[drug_ids$resource == "Wikipedia"], drug_name=drug_ids$identifier[drug_ids$resource == "Wikipedia"], stringsAsFactors=FALSE)
# Remove underscores from Wikipedia entries
t1$drug_name <- gsub("_", " ", t1$drug_name)
# Limit to "english"
t2 <- data.frame(db_id=drug_synonyms$`drugbank-id`[grepl("english", drug_synonyms$language)], drug_name=drug_synonyms$synonym[grepl("english", drug_synonyms$language)], stringsAsFactors=FALSE)
tmp_drug_ids <- rbind(t1, t2)

## General processing of drug names 
tmp_drug_ids$drug_name <- trimws(tmp_drug_ids$drug_name)

tmp_drug_ids <- unique(tmp_drug_ids)

# CLEAN DRUG NAMES ----
# From drugs and term come from: get_ctrdata.R
tmp_drugs <- drugs
tmp_drug_names <- tmp_drugs$drug_name

# Text strings used in clinicaltrials.gov intervention names and that cause problems with name matching
remove_patterns <- c("oral", "tablet", "tablets", "hydrochloride", "mesylate", 
                     "\\d+\\s?mg\\/kg", "\\d+\\s?mg", "\\d+\\s?mg/ml", "\\d+\\s?mg/day", "\\d+mg/kg/day", "sulfate", "s-malate", "hcl", 
                     "liposomal", "injection", "injectable", "mesilate", "phosphate", "sublingual spray", 
                     "malate", "plus", "capsule", "inhalant", "product", "prophylactic", 
                     "daily", "syringe", "inj", "spray", "unfractionated", "nebulized", "nebulised", "intranasal", "supplement",
                     "inhaled", "theraputic", "subcutaneous", "solution", "nasal", "aerosolized", "suspension",
                     "milligram", "intravenous", "nebuliser", "®", "\\d+.\\d+%",
                     "high", "low", "dose", "™", "dosage", "therapy", "vitro", "pegylated", "isolated", "chemotherapy", 
                     "colloidal", "dispersion", "emulsion", "neoadjuvant", "multiple", "agents",
                     "weekly", "standard", "formulation", "nanoparticle", "auc/ml\\*min", "interval", "combination of",
                     "combination", "once", "every", "\\d weeks", "in combination", "adjuvant", "surgery",
                     "administration of", "with", "\\d cycles of", "/m2", "/m^2", "/m²", "expansion", "escalation",
                     "single agent", "during surgery", "pill", "conventional", "intraperitoneal", "cream")

for(pattern in remove_patterns) {
  tmp_drug_names <- gsub(pattern, "", tmp_drug_names)
}

# Remove trailing punctuation
tmp_drug_names <- gsub("[\\.,]$", "", tmp_drug_names)
tmp_drug_names <- trimws(tmp_drug_names)

tmp_drugs$drug_name <- tmp_drug_names

# Markers indicating drug combinations
t1 <- c("+", "/", ",", "and", "&amp;", ";", "or")
t2 <- c("\\+", "/", ",", "and", "&amp;", ";", "or")
t3 <- c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE)
drug_combo_symbols <- data.frame(symbol=t1, symbol_escaped=t2, word_boundary_needed=t3, stringsAsFactors=FALSE)
      
d1 <- data.frame(
  nct_id=character(0),
  phase=character(0),
  overall_status=character(0),
  target_symbol=character(0), 
  target_id=character(0), 
  original_trial_drug_entry=character(0),
  combination_detected=logical(0),
  drug_name=character(0), 
  drug_id=character(0),
  stringsAsFactors=FALSE)

# Use regex for drug_name searching instead of exact matching
regex_flag <- TRUE
extract_combo_flag <- TRUE
extract_synonym_flag <- TRUE
min_char <- 3

for(i in 1:nrow(tmp_drugs)) {
  #i <- 1 
  cat("I: ", i, "\n")
  
  nct_id <- tmp_drugs$nct_id[i]
  drug_name <- tmp_drugs$drug_name[i]
  original_trial_drug_entry <- drug_name
  phase <- tmp_drugs$phase[i]
  overall_status = gsub(",", "", tmp_drugs$overall_status[i])
  combination_detected <- FALSE
  
  # PROCESS DRUG NAMES ----
  
  # "Weird" characters like this "）" can otherwise exist
  t1_drug_names <- iconv(drug_name, "UTF-8", "ASCII", sub=" ")

  if(extract_synonym_flag) {
    t1 <- stringr::str_match(t1_drug_names, "(.*)\\(.*\\)")[1,2]
    t2 <- stringr::str_match(t1_drug_names, "\\((.*)\\)")[1,2]

    if(!is.na(t1) && !is.na(t2) && nchar(t1) > min_char && nchar(t2) > min_char) {
      t1_drug_names <- c(t1, t2)      
    }
  }
  
  # NOTE: Do not run this block if the synonym block found something
  if(extract_combo_flag && length(t1_drug_names) == 1) {
    combo_symbol_idx <- sapply(1:nrow(drug_combo_symbols), function(i) {
      if(drug_combo_symbols$word_boundary_needed[i]) {
        grepl(paste0("\\b", drug_combo_symbols$symbol_escaped[i], "\\b"), t1_drug_names)        
      } else {
        grepl(drug_combo_symbols$symbol_escaped[i], t1_drug_names)
      }
    }, USE.NAMES=FALSE) %>% which
    
    if(length(combo_symbol_idx) > 0) {
      t1_drug_names <- stringr::str_split(t1_drug_names, drug_combo_symbols$symbol_escaped[combo_symbol_idx[1]])[[1]]
      combination_detected <- TRUE
    }
  }
  
  tmp_drug_names_processed <- sapply(t1_drug_names, function(drug_name) {
    # Remove regex meta-characters from string; square brackets [] are denoted with hex values 5B and 5D
    drug_name <- gsub("[\\(\\)\\*\\+\\?\\|\\x5B\\x5D]", " ", drug_name, perl=TRUE)
    drug_name <- trimws(drug_name)    
  }, USE.NAMES=FALSE)
  
  for(drug_name in tmp_drug_names_processed) {
    idx <- which(tolower(tmp_drug_ids$drug_name) == drug_name)    
    
    #if(length(idx) == 0 && regex_flag) {
    #  idx <- which(grepl(paste0("\\b", drug_name, "\\b"), tolower(tmp_drug_ids$drug_name)))
    #}
    
    db_id <- tmp_drug_ids$db_id[idx]
    
    if(length(db_id) > 0 && nchar(drug_name) > min_char) {
      if(length(db_id) > 1) {
        cat(paste0("Multiple DrugBank entries founds: ORG_NAME: ", original_trial_drug_entry, " NAME: ", drug_name, " DB-ID: ", db_id, "\n"))
      }
      
      db_id <- db_id[1]
      
      be_ids <- drug_targets[drug_targets$parent_key == db_id,]
      drug_gene_symbols <- drug_target_ids[drug_target_ids$resource == "GenAtlas",]
      t1 <- drug_gene_symbols[drug_gene_symbols$parent_key %in% be_ids$id, ]
      t1$drug_name <- drug_name
      t1$db_id <- db_id
      t2 <- t1[, c("identifier", "parent_key", "drug_name", "db_id")]
      colnames(t2) <- c("target_symbol", "target_id", "drug_name", "drug_id")  
      t2$nct_id <- nct_id
      t2$phase <- phase
      t2$overall_status <- overall_status
      t2$original_trial_drug_entry <- original_trial_drug_entry
      t2$combination_detected <- combination_detected
      
      d1 <- rbind(d1, t2)  
    } else {
      cat(paste0("No DrugBank entries founds: ORG_NAME: ", original_trial_drug_entry, " NAME: ", drug_name, "\n"))
    }
  }
}

# Find unique entries
d1 <- unique(d1)

# EXPORT RESULTS ----
output_file <- paste0(term, "_trial_drug_target.csv")
d1 <- d1[order(d1$target_symbol),]
write.csv(d1, output_file, quote=TRUE, row.names=FALSE)

# SUMMARIZE ----
table(d1$target_symbol) %>% sort(., decreasing=FALSE)
