# query - vector of length 4 containing the GnPS "cluster.index", "ExactMass", "IonMode", "Adduct"
# all spectra - GnPS cluster spectra mgf file
# database - local database
# type - local or MetFrag Database type (PubChem, ExtendedPubChem, ChemSpider, KEGG, LocalInChI, LocalSDF, MySQL, PostgresSQL)  
# PrecursorIonMode (1 = [M+H]⁺, 0 = [M]⁺⁻, 18 = [M+NH4]⁺, 23 = [M+Na]⁺, 39 = [M+K]⁺, -1 = [M-H]⁻, 35 = [M+Cl]⁻, 45 = [M+HCOO]⁻, 59 = [M+CH3COO]⁻)  
# either ppm (if > 0) or abs difference

#[1] "M+H"        "[M+H]"      "M-H"        "M+NH4"      "M+2H"       "M+Na"       "M+K"        "M"          "M+H-H2O"    "M-H2O+H"    "[M+Na]"     "[M+2H]"     "[M+H]+"     "579.171"   
#[15] "[M]+*"      "[M]*+"      "[M]+"       "[M+Na]+"    "[M-H2O+H]+" "[2M+H]+"    "[M+2H]2+"   "[M+H-H2O]+" "[M+NH4]+"
#  c(1,1,-1, 18,1,23,39,0,1,1,23,1,1,1,0,0,0,23,1,1,1,1,18 )

get.metfrag <- function (query, allspectra, database = NULL, type = "local", path="/PATHTOOLSFOLDER/nap_ccms2/Snap/MetFrag2.3-CL.jar", 
    abs.diff = 0.1, ppm = 10, Mdiff = 1.00727, sl_param, IsPositiveIonMode = TRUE, 
    PrecursorIonMode = 1, idx, out) 
{
    if (type == "local") {
        if (is.null(database)) {
		return(write.table("", paste0(out, ".psv")))
	} else {
        	database2 <- apply(database, 1, function(x) paste0(x, collapse = "|"))
	}
        fn1 <- paste0("selHMDB_local_inchi_file", idx, ".txt")
        write.table(c("MonoisotopicMass|InChI|SMILES|Identifier|InChIKey2|InChIKey1|MolecularFormula", 
            database2), fn1, row.names = FALSE, 
            col.names = FALSE, quote = FALSE)
        sl_param[grep("MetFragDatabaseType", sl_param)] <- sub("TYPE", 
            "LocalPSV", sl_param[grep("MetFragDatabaseType", 
                sl_param)])
        sl_param[grep("LocalDatabasePath", sl_param)] <- sub("^#", 
            "", sl_param[grep("LocalDatabasePath", sl_param)])
        
        sl_param[grep("LocalDatabasePath", sl_param)] <- sub("PATH", 
            paste0("./", fn1), sl_param[grep("LocalDatabasePath", 
                sl_param)])

        sl_param[grep("NeutralPrecursorMass", sl_param)] <- sub("MASS_TO_REPLACE", 
            query[2], sl_param[grep("NeutralPrecursorMass", sl_param)])

        sl_param[grep("PeakListPath", sl_param)] <- sub("PATH_TO_REPLACE", 
            paste0(out, ".txt"), sl_param[grep("PeakListPath", 
                sl_param)])

        sl_param[grep("SampleName", sl_param)] <- sub("NAME_TO_REPLACE", 
            paste0(out), sl_param[grep("SampleName", 
                sl_param)])

        if (!IsPositiveIonMode) 
            sl_param[grep("IsPositiveIonMode", sl_param)] <- sub("True", 
                "False", sl_param[grep("IsPositiveIonMode", sl_param)])
        if (PrecursorIonMode != 1) 
            sl_param[grep("PrecursorIonMode", sl_param)] <- sub("1", PrecursorIonMode, 
                 sl_param[grep("PrecursorIonMode", sl_param)])
        npar <- paste0(Sys.Date(), idx, "par.txt")
        write.table(sl_param, npar, row.names = FALSE, col.names = FALSE, 
            quote = FALSE)
    }
    else {
        sl_param[grep("MetFragDatabaseType", sl_param)] <- sub("TYPE", 
            type, sl_param[grep("MetFragDatabaseType", sl_param)])

        sl_param[grep("DatabaseSearchRelativeMassDeviation", 
            sl_param)] <- sub("^#", "", sl_param[grep("DatabaseSearchRelativeMassDeviation", 
            sl_param)])

        sl_param[grep("DatabaseSearchRelativeMassDeviation", 
            sl_param)] <- sub("PPM", ppm, sl_param[grep("DatabaseSearchRelativeMassDeviation", 
            sl_param)])

        sl_param[grep("NeutralPrecursorMass", sl_param)] <- sub("MASS_TO_REPLACE", 
            query[2], sl_param[grep("NeutralPrecursorMass", sl_param)])

        sl_param[grep("PeakListPath", sl_param)] <- sub("PATH_TO_REPLACE", 
            paste0(out, ".txt"), sl_param[grep("PeakListPath", 
                sl_param)])

        sl_param[grep("SampleName", sl_param)] <- sub("NAME_TO_REPLACE", 
            paste0(out), sl_param[grep("SampleName", 
                sl_param)])

        if (!IsPositiveIonMode) 
            sl_param[grep("IsPositiveIonMode", sl_param)] <- sub("True", 
                "False", sl_param[grep("IsPositiveIonMode", sl_param)])
        if (PrecursorIonMode != 1) 
            sl_param[grep("PrecursorIonMode", sl_param)] <- sub("1", PrecursorIonMode, 
                 sl_param[grep("PrecursorIonMode", sl_param)])
        npar <- paste0(Sys.Date(), idx, "par.txt")
        write.table(sl_param, npar, row.names = FALSE, col.names = FALSE, 
            quote = FALSE)
    }
    begin <- which(paste0("SCANS=", as.matrix(query[1])) == allspectra)
    end <- grep("END IONS", allspectra[begin:length(allspectra)])[1]
    write.table(allspectra[(begin + 1):(begin + end - 2)], paste0(out, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

    # Run MetFrag
    log <- system(paste("java -jar", path, npar))
	
    # Delete parameter file, spectrum and local db
    system(paste("rm", npar))
    if (type == "local") system(paste("rm", fn1))
    system(paste("rm", paste0(out, ".txt")))
}
