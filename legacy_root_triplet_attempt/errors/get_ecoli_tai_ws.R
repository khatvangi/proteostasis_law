#!/usr/bin/env Rscript

suppressMessages({
  library(tAI)
})

# Load example E. coli K-12 tRNA data
data("ecolik12", package = "tAI")

# Compute relative adaptiveness values (ws) for codons
ws <- get.ws(tRNA = ecolik12$trna, sking = 1)

# Codon order used by codonM / tAI (from the GitHub README)
codon61 <- c(
  # 1–13
  "TTT","TTC","TTA","TTG",
  "TCT","TCC","TCA","TCG",
  "TAT","TAC","TGT","TGC","TGG",
  # 14–29
  "CTT","CTC","CTA","CTG",
  "CCT","CCC","CCA","CCG",
  "CAT","CAC","CAA","CAG",
  "CGT","CGC","CGA","CGG",
  # 30–45
  "ATT","ATC","ATA","ATG",
  "ACT","ACC","ACA","ACG",
  "AAT","AAC","AAA","AAG",
  "AGT","AGC","AGA","AGG",
  # 46–61
  "GTT","GTC","GTA","GTG",
  "GCT","GCC","GCA","GCG",
  "GAT","GAC","GAA","GAG",
  "GGT","GGC","GGA","GGG"
)

# tAI ignores Met (ATG) and STOP codons; ecotik12 tutorial drops col 33 (ATG)
codon60 <- codon61[-33]

if (length(ws) != length(codon60)) {
  stop(sprintf("Length mismatch: ws=%d, codons=%d", length(ws), length(codon60)))
}

df <- data.frame(
  codon = codon60,
  w_tai = as.numeric(ws),
  stringsAsFactors = FALSE
)

out_path <- "ecoli_tai_ws.tsv"
write.table(df, file = out_path,
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("[INFO] Wrote", nrow(df), "codons with tAI weights to", out_path, "\n")

