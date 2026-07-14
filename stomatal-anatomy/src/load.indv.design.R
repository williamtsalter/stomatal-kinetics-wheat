#'# Wheat stomatal kinetics (heat & drought) — bundled design data loader
#'## Reconstructs `indv.dat` from the repository-bundled `indv.design.csv`,
#'## replacing the external `indv.v2.Rdata` produced by the water-use analysis.
#'## Column types and factor level orders are set explicitly so that models
#'## fitted downstream are identical to the original analysis.

indv.dat <- read.csv("../data/indv.design.csv", colClasses = "character",
                     check.names = FALSE, na.strings = c("", "NA"))

#'### Numeric columns
indv.dat[["WU.Total"]] <- as.numeric(indv.dat[["WU.Total"]])

#'### Character columns are left as-is (Cell.ID, Pot.ID, Treatment)

#'### Factor columns (explicit level order preserved from the original object)
indv.dat[["Room"]] <- factor(indv.dat[["Room"]], levels = c("CER7", "CER8"))
indv.dat[["Lane"]] <- factor(indv.dat[["Lane"]], levels = c("A", "B", "C", "D", "E", "F"))
indv.dat[["Position"]] <- factor(indv.dat[["Position"]], levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"))
indv.dat[["PCoord"]] <- factor(indv.dat[["PCoord"]], levels = c("1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5", "6", "6.5", "7", "8", "8.5", "9", "9.5", "10", "10.5", "11", "11.5", "12", "12.5", "13", "13.5", "14"))
indv.dat[["Side"]] <- factor(indv.dat[["Side"]], levels = c("N", "S"))
indv.dat[["Location"]] <- factor(indv.dat[["Location"]], levels = c("Inner", "Outer"))
indv.dat[["Block"]] <- factor(indv.dat[["Block"]], levels = c("1", "2"))
indv.dat[["SubUnit"]] <- factor(indv.dat[["SubUnit"]], levels = c("1", "2"))
indv.dat[["Cell"]] <- factor(indv.dat[["Cell"]], levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
indv.dat[["Replicate"]] <- factor(indv.dat[["Replicate"]], levels = c("1", "2", "3", "4"))
indv.dat[["Rm_Temperature"]] <- factor(indv.dat[["Rm_Temperature"]], levels = c("7,Heat", "8,Ambient"))
indv.dat[["Watering"]] <- factor(indv.dat[["Watering"]], levels = c("Limited", "Unlimited"))
indv.dat[["Line"]] <- factor(indv.dat[["Line"]], levels = c("SUNCHASER", "ROCKSTAR", "PBI17N006-0N-0N-89N", "PBI17N013-0C-0N-0N-010N-8N", "PBIC15034-0C-24N-010N-1N-0N", "PBI16C009-0C-0N-0N-010N-14N", "PBI19N007-0N-3N", "PBI17N007-0C-0N-0N-010N-8N", "PBI17N016-0C-0N-0N-010N-27N", "PBI17N002-0C-0N-0N-010N-16N"))
indv.dat[["Line.ID"]] <- factor(indv.dat[["Line.ID"]], levels = c("L01", "L02", "L03", "L04", "L05", "L06", "L07", "L08", "L09", "L10"))
indv.dat[["AWatering"]] <- factor(indv.dat[["AWatering"]], levels = c("rest", "Limited", "Unlimited"))
indv.dat[["HWatering"]] <- factor(indv.dat[["HWatering"]], levels = c("rest", "Limited", "Unlimited"))
indv.dat[["ALine"]] <- factor(indv.dat[["ALine"]], levels = c("rest", "SUNCHASER", "ROCKSTAR", "PBI17N006-0N-0N-89N", "PBI17N013-0C-0N-0N-010N-8N", "PBIC15034-0C-24N-010N-1N-0N", "PBI16C009-0C-0N-0N-010N-14N", "PBI19N007-0N-3N", "PBI17N007-0C-0N-0N-010N-8N", "PBI17N016-0C-0N-0N-010N-27N", "PBI17N002-0C-0N-0N-010N-16N"))
indv.dat[["HLine"]] <- factor(indv.dat[["HLine"]], levels = c("rest", "SUNCHASER", "ROCKSTAR", "PBI17N006-0N-0N-89N", "PBI17N013-0C-0N-0N-010N-8N", "PBIC15034-0C-24N-010N-1N-0N", "PBI16C009-0C-0N-0N-010N-14N", "PBI19N007-0N-3N", "PBI17N007-0C-0N-0N-010N-8N", "PBI17N016-0C-0N-0N-010N-27N", "PBI17N002-0C-0N-0N-010N-16N"))
indv.dat[["TW"]] <- factor(indv.dat[["TW"]], levels = c("Ambient-Limited", "Ambient-Unlimited", "Heat-Limited", "Heat-Unlimited"))
indv.dat[["WT"]] <- factor(indv.dat[["WT"]], levels = c("Limited-Ambient", "Limited-Heat", "Unlimited-Ambient", "Unlimited-Heat"))
indv.dat[["ASide"]] <- factor(indv.dat[["ASide"]], levels = c("N", "S"))
indv.dat[["ABlock"]] <- factor(indv.dat[["ABlock"]], levels = c("1", "2"))

#'### Sanity check
stopifnot(nrow(indv.dat) == 160L, !anyDuplicated(indv.dat[c("Watering","Replicate","Rm_Temperature","Line.ID")]))

