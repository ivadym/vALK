# Libraries ---------------------------------------------------------------

options(repos=structure(c(CRAN="https://cran.rediris.es/")))

list.of.packages <- c("tcltk", "scales")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

suppressPackageStartupMessages({
  if (length(new.packages)) install.packages(new.packages, quiet = TRUE)
  library(tcltk) # TCL/TK package
  library(scales) # Scales package
})

# Constants ---------------------------------------------------------------

genes <- c("ALK") # Gene options

dest_cols_tsv <- c("Sample_ID", "Sample", "Extraction method", "Sample type",
               "Treatment progression", "NGS date","FUNC1.gene", "FUNC1.coding",
               "FUNC1.protein", "FUNC1.transcript", "INFO.1.DP", "INFO.1.FDP", 
               "INFO.1.FRO", "INFO.A.AO", "INFO.A.FAO", "INFO.1.FD", "rowtype",
               "call", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
               "INFO.1.FSRF", "INFO.1.FSRR", "INFO.A.FSAF", "INFO.A.FSAR",
               "INFO.A.LOD", "INFO.A.MAF", "FUNC1.codon", "FUNC1.exon", 
               "FUNC1.function", "FUNC1.location", "FUNC1.oncomineGeneClass",
               "FUNC1.oncomineVariantClass", "FUNC1.CLNSIG1", "MAPD",
               "Median read coverage", "Median molecular coverage",
               "Limits of detection (%)", "dPCR date", "MAF", "Repeated dPCR date",
               "Repeated MAF")

dest_cols_names <- c("Sample_ID", "Sample", "Extraction method", "Sample type",
               "Treatment progression", "NGS date","Gene", "Coding transcript change",
               "Protein change", "Transcript", "Total coverage", "Molecular coverage", 
               "WT molecular coverage", "Mutation reads", "Mutation molecular coverage",
               "CNV ratio", "Rowtype", "Call", "Chromosome", "Position", "ID", "Reference allele",
               "Alternate allele", "Quality", "Filter", "Reference allele forward strand",
               "Reference allele reverse strand", "Alternate allele forward strand",
               "Alternate allele reverse strand", "Detection limit (%)", "Minor allele frequency",
               "Codon", "Exon", "Functional classification", "Location", "Oncomine gene class",
               "Oncomine variant class", "Variant clinical significance", 
               "MAPD (Median Absolute Pairwise Difference)", "Median read coverage",
               "Median molecular coverage", "Limits of detection (%)", "dPCR date", "MAF",
               "Repeated dPCR date", "Repeated MAF")

lockBinding("genes", globalenv())
lockBinding("dest_cols_tsv", globalenv())
lockBinding("dest_cols_names", globalenv())

# Arguments ---------------------------------------------------------------

gene <- ""
var_path <- ""
dest_path <- ""

# Functions ---------------------------------------------------------------

# Extracts the possible fusions
# @param variants Genetic variants data frame to be analyzed
# @param gene     Gene to filter by
#
get.fusions <- function (variants, gene) {
  if (as.character(gene) == "ALK") {
    check_fus <- variants[as.character(variants$rowtype) == "ProcControl" & !is.na(variants$rowtype),]
    check_fus <- check_fus[as.character(check_fus$FILTER) == "PASS" & !is.na(check_fus$FILTER),]
    if (nrow(check_fus) > 1) {
      fus <- variants[as.character(variants$FUNC1.gene) == gene & !is.na(variants$FUNC1.gene),]
      fus <- fus[as.character(fus$rowtype) == "Fusion" & !is.na(fus$rowtype),]
      fus <- fus[as.numeric(as.character(fus$INFO...READ_COUNT)) > 25 &
                   !is.na(fus$INFO...READ_COUNT),]
      fus <- fus[as.numeric(as.character(fus$INFO...MOL_COUNT)) > 2 &
                   !is.na(fus$INFO...MOL_COUNT),]
    }
  }
}

MAPD = NULL # Median Absolute Pairwise Difference

# Extracts the copy number variations (duplication or deletion)
# @param variants Genetic variants data frame to be analyzed
# @param gene     Gene to filter by
#
get.CNVs <- function (variants, gene) {
  if (as.character(gene) == "ALK") {
    MAPD <<- as.numeric(as.character(read.csv(var_path, header = FALSE, sep = "=")[[2]][[4]]))
    if ((MAPD < 0.4) & !is.null(MAPD)) {
      cnvs <- variants[as.character(variants$FUNC1.gene) == gene & !is.na(variants$FUNC1.gene),]
      cnvs <- cnvs[as.character(cnvs$rowtype) == "CNV" & !is.na(cnvs$rowtype),]
      cnvs <- cnvs[(as.character(cnvs$FILTER) == "GAIN" | as.character(cnvs$FILTER) == "LOSS") & 
                     !is.na(cnvs$FILTER),]
    }
  }
}

# Extracts the single-nucleotide polymorphisms (SNPs)
# @param variants Genetic variants data frame to be analyzed
# @param gene     Gene to filter by
#
get.SNPs <- function (variants, gene) {
  if (as.character(gene) == "ALK") {
    
    # Applies the conditions c1
    # @param vars Variants to filter
    #
    get.c1 <- function(vars) {
      vars <- vars[as.character(vars$FILTER) == "PASS" & !is.na(vars$FILTER),]
      vars <- vars[as.character(vars$FUNC1.oncomineVariantClass) == "Hotspot" &
                     !is.na(vars$FUNC1.oncomineVariantClass),]
      vars <- vars[as.numeric(as.character(vars$INFO.A.FAO)) >= 1 &
                     !is.na(vars$INFO.A.FAO),]
      vars <- vars[as.numeric(as.character(vars$INFO.A.AO)) >= 15 &
                     !is.na(vars$INFO.A.AO),]
    }
    
    # Applies the conditions c2
    # @param vars Variants to filter
    #
    get.c2 <- function(vars) {
      vars <- vars[as.character(vars$FILTER) == "PASS" & !is.na(vars$FILTER),]
      vars <- vars[as.character(vars$FUNC1.oncomineVariantClass) == "Hotspot" &
                     !is.na(vars$FUNC1.oncomineVariantClass),]
      vars <- vars[as.numeric(as.character(vars$INFO.A.FAO)) >= 2 &
                     !is.na(vars$INFO.A.FAO),]
    }
    
    # Applies the conditions c3
    # @param vars Variants to filter
    #
    get.c3 <- function(vars) {
      vars <- vars[as.numeric(as.character(vars$INFO.A.FAO)) >= 2 &
                     !is.na(vars$INFO.A.FAO),]
      vars <- vars[as.numeric(as.character(vars$INFO.A.AO)) >= 25 &
                     !is.na(vars$INFO.A.AO),]
    }
    
    snps <- variants[as.character(variants$FUNC1.gene) == gene & !is.na(variants$FUNC1.gene),]
    snps <- snps[as.character(snps$rowtype) == "snp" & !is.na(snps$rowtype),]
    snps <- if ("FUNC1.CLNSIG1" %in% names(snps)) snps[((as.character(snps$FUNC1.CLNSIG1) != "Benign" &
                as.character(snps$FUNC1.CLNSIG1) != "Benign/Likely_benign") & !is.na(snps$FUNC1.CLNSIG1))
                                                       ,] else snps
   snps_c1 <- get.c1(snps)
   snps_c2 <- get.c2(snps)
   snps_c3 <- get.c3(snps)
   snps <- merge(snps_c1, snps_c2, all = TRUE)
   snps <- merge(snps, snps_c3, all = TRUE)
  }
}

# Extracts the insertion or deletion of bases 
# @param variants Genetic variants data frame to be analyzed
# @param gene     Gene to filter by
#
get.InDels <- function (variants, gene) {
  if (as.character(gene) == "ALK") {
    indels <- variants[variants$FUNC1.gene == gene & !is.na(variants$FUNC1.gene),]
    indels <- indels[(as.character(indels$rowtype) == "ins" | 
                        as.character(indels$rowtype) == "del") & !is.na(indels$rowtype),]
    indels <- indels[as.numeric(as.character(indels$INFO.A.FAO)) >= 25 &
                       !is.na(indels$INFO.A.FAO),]
    indels <- indels[as.numeric(as.character(indels$INFO.A.AO)) >= 200 &
                       !is.na(indels$INFO.A.AO),]
    indels <- indels[as.numeric(as.character(indels$INFO.A.AF)) > 0.03 &
                       !is.na(indels$INFO.A.AF),]
  }
}

# Extracts the multiple nucleotide polymorphisms (MNPs)
# @param variants Genetic variants data frame to be analyzed
# @param gene     Gene to filter by
#
get.MNPs <- function (variants, gene) {
  if (as.character(gene) == "ALK") {
    mnps <- variants[as.character(variants$FUNC1.gene) == gene & !is.na(variants$FUNC1.gene),]
    mnps <- mnps[as.character(mnps$rowtype) == "mnp" & !is.na(mnps$rowtype),]
    mnps <- mnps[as.numeric(as.character(mnps$INFO.A.FAO)) >= 50 &
                   !is.na(mnps$INFO.A.FAO),]
    mnps <- mnps[as.numeric(as.character(mnps$INFO.A.AO)) >= 200 &
                   !is.na(mnps$INFO.A.AO),]
    mnps <- mnps[as.numeric(as.character(mnps$INFO.A.AF)) > 0.03 &
                   !is.na(mnps$INFO.A.AF),]
    mnps <- mnps[as.numeric(as.character(mnps$POS)) != 29443609 &
                   as.numeric(as.character(mnps$POS)) != 29443610 &
                   as.numeric(as.character(mnps$POS)) != 29443611 &
                   !is.na(mnps$POS),] # 6 consecutive guanines
  }
}

# Saves in the specified .csv file the possible mutations
# @param file_path  Path of the destination .csv file
# @param mutations  Mutations to be saved
#
save.mutations <- function (file_path, mutations) {
  if (file.exists(file_path)) {
    write.table(mutations, file = file_path, append = TRUE, quote = TRUE, sep = ";", 
                eol = "\n", na = "", dec = ",", row.names = FALSE, col.names = FALSE)
  } else {
    write.table(mutations, file = file_path, append = FALSE, quote = TRUE, sep = ";", 
                eol = "\n", na = "", dec = ",", row.names = FALSE, col.names = TRUE)  
  }
}

# Displays the users interface and collects the input variables
#
get.variables<- function () {
  
  # Confirms the variable selection
  #
  accept <- function () {
    tkdestroy(tt)
  }
  
  # Cancels the variable input
  #
  cancel <- function () {
    tkdestroy(tt)
    # quit() # Execution from the terminal
  }
  
  # Opens a dialog to select a variants .tsv file
  #
  select.variants <- function() {
    Filters <- matrix(c("Tab-separated values", ".tsv"), 1, 2, byrow = TRUE)
    var <<- tclVar(tk_choose.files(caption = "Select the variants file", multi = FALSE, 
                                   filters = Filters))
    tkconfigure(v_f_e, textvariable = var)
  }
  
  # Opens a dialog to select a destination .csv file
  #
  select.destination <- function() {
    Filters <- matrix(c("Comma-separated values", ".csv"), 1, 2, byrow = TRUE)
    dest <<- tclVar(tk_choose.files(caption = "Select the output file", multi = FALSE,
                                    filters = Filters))
    tkconfigure(d_e, textvariable = dest)
  }
  
  # Updates the fusions filter selection and visualization flags
  #
  check.fus <- function () {
    value_fus <<- !value_fus
    value_var_prev <- value_var
    
    if (!value_var_prev & value_fus) {
      value_var <<- TRUE
    } else if (value_var_prev & !value_fus & !value_CNVs & !value_SNPs & 
               !value_InDels & !value_MNPs) {
      value_var <<- FALSE
    }
    
    if (!value_var_prev & value_var) {
      tkgrid(v_p_frm, column = 2, row = 1, sticky = 'nw', padx = 10)
      tkgrid(v_f_frm, column = 3, row = 1, sticky = 'nw')
      tkgrid(d_frm, column = 3, row = 2, sticky = 'nw')
      tkpack(ab_frm, side = "top", pady = 10)
    } else if (value_var_prev & !value_var) {
      tkgrid.forget(v_p_frm)
      tkgrid.forget(v_f_frm)
      tkgrid.forget(d_frm)
      tkpack.forget(ab_frm)
    }
  }
  
  # Updates the CNVs filter selection and visualization flags
  #
  check.CNVs <- function () {
    value_CNVs <<- !value_CNVs
    value_var_prev <- value_var
    
    if (!value_var_prev & value_CNVs) {
      value_var <<- TRUE
    } else if (value_var_prev & !value_fus & !value_CNVs & !value_SNPs & 
               !value_InDels & !value_MNPs) {
      value_var <<- FALSE
    }
    
    if (!value_var_prev & value_var) {
      tkgrid(v_p_frm, column = 2, row = 1, sticky = 'nw', padx = 10)
      tkgrid(v_f_frm, column = 3, row = 1, sticky = 'nw')
      tkgrid(d_frm, column = 3, row = 2, sticky = 'nw')
      tkpack(ab_frm, side = "top", pady = 10)
    } else if (value_var_prev & !value_var) {
      tkgrid.forget(v_p_frm)
      tkgrid.forget(v_f_frm)
      tkgrid.forget(d_frm)
      tkpack.forget(ab_frm)
    }
  }
  
  # Updates the SNPs filter selection and visualization flags
  #
  check.SNPs <- function () {
    value_SNPs <<- !value_SNPs
    value_var_prev <- value_var
    
    if (!value_var_prev & value_SNPs) {
      value_var <<- TRUE
    } else if (value_var_prev & !value_fus & !value_CNVs & !value_SNPs & 
               !value_InDels & !value_MNPs) {
      value_var <<- FALSE
    }
    
    if (!value_var_prev & value_var) {
      tkgrid(v_p_frm, column = 2, row = 1, sticky = 'nw', padx = 10)
      tkgrid(v_f_frm, column = 3, row = 1, sticky = 'nw')
      tkgrid(d_frm, column = 3, row = 2, sticky = 'nw')
      tkpack(ab_frm, side = "top", pady = 10)
    } else if (value_var_prev & !value_var) {
      tkgrid.forget(v_p_frm)
      tkgrid.forget(v_f_frm)
      tkgrid.forget(d_frm)
      tkpack.forget(ab_frm)
    }
  }
  
  # Updates the InDels filter selection and visualization flags
  #
  check.InDels <- function () {
    value_InDels <<- !value_InDels
    value_var_prev <- value_var
    
    if (!value_var_prev & value_InDels) {
      value_var <<- TRUE
    } else if (value_var_prev & !value_fus & !value_CNVs & !value_SNPs & 
               !value_InDels & !value_MNPs) {
      value_var <<- FALSE
    }
    
    if (!value_var_prev & value_var) {
      tkgrid(v_p_frm, column = 2, row = 1, sticky = 'nw', padx = 10)
      tkgrid(v_f_frm, column = 3, row = 1, sticky = 'nw')
      tkgrid(d_frm, column = 3, row = 2, sticky = 'nw')
      tkpack(ab_frm, side = "top", pady = 10)
    } else if (value_var_prev & !value_var) {
      tkgrid.forget(v_p_frm)
      tkgrid.forget(v_f_frm)
      tkgrid.forget(d_frm)
      tkpack.forget(ab_frm)
    }
  }
  
  # Updates the MNPs filter selection and visualization flags
  #
  check.MNPs <- function () {
    value_MNPs <<- !value_MNPs
    value_var_prev <- value_var
    
    if (!value_var_prev & value_MNPs) {
      value_var <<- TRUE
    } else if (value_var_prev & !value_fus & !value_CNVs & !value_SNPs & 
               !value_InDels & !value_MNPs) {
      value_var <<- FALSE
    }
    
    if (!value_var_prev & value_var) {
      tkgrid(v_p_frm, column = 2, row = 1, sticky = 'nw', padx = 10)
      tkgrid(v_f_frm, column = 3, row = 1, sticky = 'nw')
      tkgrid(d_frm, column = 3, row = 2, sticky = 'nw')
      tkpack(ab_frm, side = "top", pady = 10)
    } else if (value_var_prev & !value_var) {
      tkgrid.forget(v_p_frm)
      tkgrid.forget(v_f_frm)
      tkgrid.forget(d_frm)
      tkpack.forget(ab_frm)
    }
  }
  
  # Framework definition
  tt <- tktoplevel()
  tkwm.title(tt, "NGS filtering")
  
  # Input frame
  t_frm <- tkframe(tt, padx = 15, pady = 15)
  
  # Variants frame
  var <- tclVar("Variants file (.tsv)")
  gene <- tclVar("ALK")
  
  v_frm <- ttklabelframe(t_frm, text = "Variants", padding = 10)
  
  # Checkboxes
  c_frm <- ttklabelframe(v_frm, text = "Filter", padding = 10)
  value_var <- TRUE
  
  # Fusions checkbox
  label_fus <- tclVar("Fusions")
  value_fus <- TRUE

  tkpack(ttkcheckbutton(c_frm, variable = tclVar(value_fus), textvariable = label_fus, 
                        command = check.fus), side = "top", anchor = "nw")
  
  # CNVs checkbox
  label_CNVs <- tclVar("CNVs")
  value_CNVs <- TRUE
  tkpack(ttkcheckbutton(c_frm, variable = tclVar(value_CNVs), textvariable = label_CNVs, 
                        command = check.CNVs), side = "top", anchor = "nw")
  
  # SNPs checkbox
  label_SNPs <- tclVar("SNPs")
  value_SNPs <- TRUE
  tkpack(ttkcheckbutton(c_frm, variable = tclVar(value_SNPs), textvariable = label_SNPs, 
                        command = check.SNPs), side = "top", anchor = "nw")
  
  # InDels checkbox
  label_InDels <- tclVar("InDels")
  value_InDels <- TRUE
  tkpack(ttkcheckbutton(c_frm, variable = tclVar(value_InDels), textvariable = label_InDels, 
                        command = check.InDels), side = "top", anchor = "nw")
  
  # MNPs checkbox
  label_MNPs <- tclVar("MNPs")
  value_MNPs <- TRUE
  tkpack(ttkcheckbutton(c_frm, variable = tclVar(value_MNPs), textvariable = label_MNPs, 
                        command = check.MNPs), side = "top", anchor = "nw")
  
  tkgrid(c_frm, column = 1, row = 1,  rowspan = 2, sticky = 'nw')
  
  # Parameteres
  v_p_frm <- ttklabelframe(v_frm, text = "Parameters", padding = 10)
  
  # Gene selection
  tkpack(tklabel(v_p_frm, text = "Gene", padx = 5), side = "left")
  tkpack(ttkcombobox(v_p_frm, width = 15, values = genes, textvariable = gene, 
                     state = "readonly"), side = "left")
  
  tkgrid(v_p_frm, column = 2, row = 1, rowspan = 2, sticky = 'nw', padx = 10)
  
  # Variants file
  v_f_frm <- ttklabelframe(v_frm, text = "Variants file", padding = 10)
  v_f_e <- tkentry(v_f_frm, width = 50, textvariable = var)
  tkpack(v_f_e, side = "left")
  tkpack(tkbutton(v_f_frm, text = "Select file", command = select.variants), side = "left")
  tkgrid(v_f_frm, column = 3, row = 1, sticky = 'nw')
  
  # Destination frame
  dest <- tclVar("Output file (.csv)")
  d_frm <- ttklabelframe(v_frm, text = "Output file", padding = 10)
  d_e <- tkentry(d_frm, width = 50, textvariable = dest)
  tkpack(d_e, side = "left")
  tkpack(tkbutton(d_frm, text = "Select file", command = select.destination), side = "left")
  tkgrid(d_frm, column = 3, row = 2, sticky = 'nw', pady = 10)
  
  tkgrid(v_frm)
  
  tkpack(t_frm)
  
  # Buttons frame
  bfrm <- tkframe(tt)
  ab_frm <- tkbutton(bfrm, text = "Accept", command = accept, width = 15)
  tkpack(ab_frm, side = "left", padx = 10, pady = 10)
  cb_frm = tkbutton(bfrm, text = "Cancel", command = cancel, width = 15)
  tkpack(cb_frm, side = "right", padx = 10, pady = 10)
  tkpack(bfrm, side = "bottom")
  
  tkwait.window(tt)
  
  c(value_fus, value_CNVs, value_SNPs, value_InDels, value_MNPs, 
    tclvalue(var), tclvalue(dest), tclvalue(gene))
}

# Customizes the displayed final columns
# @param var_path Path of the variants .tsv file
# @param variants Data frame of the filtered variants
#
customize.cols <- function (var_path, variants) {
  
  # Assigns the column Sample_ID the sample specified by the filepath
  # @param path Path of the variants .tsv file
  # @param df   Data frame of the filtered variants
  #
  customize.ID <- function (path, df) {
    id_sample <- strsplit(strsplit(path, "/\\s*(?=[^/]+$)", # Split by the last "/"
                                   perl = "TRUE")[[1]][[2]], "_", fixed = "TRUE") # Split by "_"
    id_sample <- paste(id_sample[[1]][[1]],  id_sample[[1]][[2]])
    df$Sample_ID <- id_sample
    df
  }
  
  # Sets the destination format (numeric or character)
  # @param df   Data frame of the filtered variants
  #
  customize.format <- function (df) {
    df$INFO.A.LOD <- percent(as.numeric(as.character(df$INFO.A.LOD)), # Detection limit
                             accuracy = 0.00000001, decimal.mark = ",")
    df
  }
  
  # Assigns the MAPD the corresponding value extracted from the .tsv file
  # @param df   Data frame of the filtered variants
  #
  assign.MAPD <- function(df) {
    df$MAPD <- MAPD
    df
  }
  
  df <- as.data.frame(matrix(NA, 0, length(dest_cols_tsv)))
  names(df) <- dest_cols_tsv
  df <- merge(df, variants[, intersect(names(df), names(variants))], all = TRUE)[, dest_cols_tsv] # .tsv column names
  df <- customize.format(df)
  df <- assign.MAPD(df)
  names(df) <- dest_cols_names # Custom column names
  df <- customize.ID(var_path, df)
}

# Shows the results after filtering the variants
# @param fus    Fusions data frame
# @param CNVs   CNVs data frame
# @param SNPs   SNPs data frame
# @param InDels InDels data frame
# @param MNPs   MNPs data frame
#
show.results <- function (fus, CNVs, SNPs, InDels, MNPs) {
  
  # Accept button function on the results screen
  #
  accept <- function () {
    tkdestroy(tt)
    # quit() # Execution from the terminal
  }
  
  if ((is.null(fus) | identical(nrow(fus), integer(1))) & 
      (is.null(CNVs) | identical(nrow(CNVs), integer(1))) & 
      (is.null(SNPs) | identical(nrow(SNPs), integer(1))) &
      (is.null(InDels) | identical(nrow(InDels), integer(1))) &
      (is.null(MNPs) | identical(nrow(MNPs), integer(1)))) {
    msgBox <- tkmessageBox(title = "Filtering results", 
                           message = "No mutation has been filtered with the specified parameters", 
                           icon = "info", type = "ok")
  } else {
    tt <- tktoplevel()
    tkwm.title(tt, "Filtering results")
    
    t_frm <- tkframe(tt, padx = 15, pady = 10)
    v_frm <- ttklabelframe(t_frm, text = "Variants", padding = 10)
    
    # Displays the filtered variants row number 
    # @param  row Row to be printed
    #
    print.vcf_rn <- function (row) {
      tkpack(tklabel(v_frm, text = paste("Row Number (vcf.rownum): ", row[1], 
                                         "     Alternate Allele ID (ALT.idx): ", row[8])), 
             side='top', anchor = "nw", padx = 15)
    }
    
    # Fusions
    if (!is.null(fus) & !identical(nrow(fus), integer(1))) {
      if (nrow(fus) > 1) {
        tkpack(tklabel(v_frm, text = paste("Fusions: ", nrow(fus), " variants have been filtered")), 
               side='top', anchor = "nw", pady = 5)
      } else {
        tkpack(tklabel(v_frm, text = paste("Fusions: ", nrow(fus), " variant has been filtered")), 
               side='top', anchor = "nw", pady = 5)
      }
      
      apply(fus, 1, print.vcf_rn) 
    }
    
    # CNVs
    if (!is.null(CNVs) & !identical(nrow(CNVs), integer(1))) {
      if (nrow(CNVs) > 1) {
        tkpack(tklabel(v_frm, text = paste("CNVs: ", nrow(CNVs), " variants have been filtered")), 
               side='top', anchor = "nw", pady = 5)
      } else {
        tkpack(tklabel(v_frm, text = paste("CNVs: ", nrow(CNVs), " variant has been filtered")), 
               side='top', anchor = "nw", pady = 5)
      }
      
      apply(CNVs, 1, print.vcf_rn)
    }
      
    # SNPs
    if (!is.null(SNPs) & !identical(nrow(SNPs), integer(1))) {
      if (nrow(SNPs) > 1) {
        tkpack(tklabel(v_frm, text = paste("SNPs: ", nrow(SNPs), " variants have been filtered")), 
               side='top', anchor = "nw", pady = 5)
      } else {
        tkpack(tklabel(v_frm, text = paste("SNPs: ", nrow(SNPs), " variant has been filtered")), 
               side='top', anchor = "nw", pady = 5)
      }
      
      apply(SNPs, 1, print.vcf_rn)
    }
    
    # InDels
    if (!is.null(InDels) & !identical(nrow(InDels), integer(1))) {
      if (nrow(InDels) > 1) {
        tkpack(tklabel(v_frm, text = paste("InDels: ", nrow(InDels), " variants have been filtered")), 
               side='top', anchor = "nw", pady = 5)
      } else {
        tkpack(tklabel(v_frm, text = paste("InDels: ", nrow(InDels), " variant has been filtered")), 
               side='top', anchor = "nw", pady = 5)
      }
      
      apply(InDels, 1, print.vcf_rn)
    }
    
    # MNPs
    if (!is.null(MNPs) & !identical(nrow(MNPs), integer(1))) {
      if (nrow(MNPs) > 1) {
        tkpack(tklabel(v_frm, text = paste("MNPs: ", nrow(MNPs), " variants have been filtered")), 
               side='top', anchor = "nw", pady = 5)
      } else {
        tkpack(tklabel(v_frm, text = paste("MNPs: ", nrow(MNPs), " variant has been filtered")), 
               side='top', anchor = "nw", pady = 5)
      }
      
      apply(MNPs, 1, print.vcf_rn)
    }
      
    tkpack(v_frm, side = "top", anchor = "nw", pady = 10)
    
    tkpack(t_frm)
    
    # Buttons frame
    bfrm <- tkframe(tt)
    ab_frm <- tkbutton(bfrm, text = "Accept", command = accept, width = 15)
    tkpack(ab_frm, side = "top", pady = 10)
    tkpack(bfrm, side = "bottom")
    
    tkwait.window(tt)
  }
}
# Main --------------------------------------------------------------------

fusions <- FALSE
CNVs <- FALSE
SNPs <- FALSE
InDels <- FALSE
MNPs <- FALSE

# Collection of variables from the dialog
inputs <- get.variables()
fusions <- as.logical(inputs[1])
CNVs <- as.logical(inputs[2])
SNPs <- as.logical(inputs[3])
InDels <- as.logical(inputs[4])
MNPs <- as.logical(inputs[5])
var_path <<- inputs[6]
dest_path <<- inputs[7]
gene <<- inputs[8]

# Progress message
tt_p <- tktoplevel()
tkwm.title(tt_p, "Filtering results")
t_frm <- tkframe(tt_p, padx = 15, pady = 15)
tkpack(tklabel(t_frm, text = "Please wait for the filtering to finish..."), side = "top")
tkpack(t_frm)

# Creation and saving of the data frames
result_df <- NULL
fus_df <- NULL
CNVs_df <- NULL
SNPs_df <- NULL
InDels_df <- NULL
MNPs_df <- NULL

if (fusions | CNVs | SNPs | InDels | MNPs) {
  variants_df <- read.table(var_path, header = TRUE, dec = ",")
  
  if (fusions) {
    fus_df <<- get.fusions(variants_df, gene)
    result_df <<- merge(result_df, fus_df, all = TRUE)
  }

  if (CNVs) {
    CNVs_df <<- get.CNVs(variants_df, gene)
    result_df <<- merge(result_df, CNVs_df, all = TRUE)
  }
  
  if (SNPs)  {
    SNPs_df <<- get.SNPs(variants_df, gene)
    result_df <<- merge(result_df, SNPs_df, all = TRUE)
  }
  
  if (InDels)  {
    InDels_df <<- get.InDels(variants_df, gene)
    result_df <<- merge(result_df, InDels_df, all = TRUE)
  }
  
  if (MNPs)  {
    MNPs_df <<- get.MNPs(variants_df, gene)
    result_df <<- merge(result_df, MNPs_df, all = TRUE)
  }
  
  if (nrow(result_df) > 0) {
    save.mutations(dest_path, customize.cols(var_path, result_df))
  }
}

# Progress message
tkdestroy(tt_p)

# Shows the filtering results
show.results(fus_df, CNVs_df, SNPs_df, InDels_df, MNPs_df)