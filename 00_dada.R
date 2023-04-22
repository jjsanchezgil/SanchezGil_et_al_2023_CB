require(dada2, quietly = TRUE)
require(stringr, quietly = TRUE)
require(assertthat, quietly = TRUE)
require(Biostrings, quietly = TRUE)
require(magrittr, quietly = TRUE)

##Define the sequences in the study
sq_SAL <- "TGGGGAATCTTGCACAATGGGGTCACCCCTGATGCAGCCATGCCGCGTGGAGGAAGACACCCCTATGGGGCGTAAACTCCTTTTCTGAATGAAGAAACCCCTGTAGCTTCAGGGCGCGACGGTAGTTCAGGAATAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGTGCAGGCGGGGCAGCAAGTCGGATGTGAAACCCCATGGCTTAACCATGGAGGTGCATTCGAAACTGTTGCTCTTGAGTCCCGGAGAGGCTGTCGGAATTCGTGGTGTAGCGGTGAAATGCGTAGATATCACGAGGAACACCAGAGGCGAAAGCGGACAGCTGGACGGGTACTGACGCTCAGGCACGAAAGCGTGGGGAGCAAACA"
sq_CLP <- "TGGGGAATTTTCCGCAATGGGCGAAAGCCTGACGGAGCAATGCCGCGTGGAGGTAGAAGGCCTACGGGTCCTGAACTTCTTTTCCCAGAGAAGAAGCAATGACGGTATCTGGGGAATAAGCATCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACCAAGCTTGAGTACGGTAGGGGCAGAGGGAATTTCCGGTGGAGCGGTGAAATGCGTAGAGATCGGAAAGAACACCAACGGCGAAAGCACTCTGCTGGGCCGACACTGACACTGAGAGACGAAAGCTAGGGGAGCGAATG"
sq_MIT <- "TGGGGAATCTTGGACAATGGGCGAAAGCCCGATCCAGCAATATCGCGTGAGTGAAGAAAGGCAATGCCGCTTGTAAAGCTCTTTCGTCGAGTGCGCGATCATGACAGGACTCGAGGAAGAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAAAACGGGGGGGGCAAGTGTTCTTCGGAATGACTGGGCGTAAAGGGCACGTAGGCGGTGAATCGGGTTGAAAGTGAAAGTCGCCAAAAAGTGGCGGAATGCTTTCGAAACCAATTCACTTGAGTGAGACAGAGGAGAGTGGAATTTCGTGTGGAGGGGTGAAATCTACAGATCTACGAAGGAACGCCAAAAGCGAAGGCAGCTCTCTGGGTCCCTACCGACGCTGGGGTGCGAAAGCATGGGGAGCGAACG"
sq_CHA <- "TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCACTTTAAGTTGGGAGGAAGGGCAGTTACCTAATACGTGATTGTTTTGACGTTACCGACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTGGCAAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACA"
sq_417 <- "TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCACTTTAAGTTGGGAGGAAGGGCAGTTACCTAATACGTGATTGTTTTGACGTTACCGACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACTGACTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTAATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACA"
sq_365 <- "TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCACTTTAAGTTGGGAGGAAGGGCAGTAAGCGAATACCTTGCTGTTTTGACGTTACCGACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGAATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTGGCAAGCTAGAGTACGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACA"
sq_358 <- "TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCACTTTAAGTTGGGAGGAAGGGCAGTAAGCGAATACCTTGCTGTTTTGACGTTACCGACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGAATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTGGCAAGCTAGAGTACGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACA"
sq_317 <- "TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCACTTTAAGTTGGGAGGAAGGGTTGTAGATTAATACTCTGCAATTTTGACGTTACCGACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACA"
sq_158 <- "TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCACTTTAAGTTGGGAGGAAGGGCAGTAAGTTAATACCTTGCTGTTTTGACGTTACCGACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTCGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTGGCGAGCTAGAGTACGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACA"
sq_134 <- "TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCACTTTAAGTTGGGAGGAAGGGTTGTAGATTAATACTCTGCAATTTTGACGTTACCGACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTCGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGTCGAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACA"
sq_V3V4 <- c("SAL" = sq_SAL, "CLP" = sq_CLP, "MIT" = sq_MIT,
             "358" = sq_358, "417" = sq_417, "CHA" = sq_CHA,
             "158" = sq_158, "365" = sq_365, "134" = sq_134,
             "317" = sq_317)

####truncSeqs: truncate all sequences in vector to a given length
#param seqVector: vector or list of sequences
#param trunc_len: desired length to truncate each sequence
#param side: side to truncate, from the start or the end of the sequence
#param rc: if TRUE, reverse complement the sequences (for reverse reads)
truncSeqs <- function(seqVector, trunc_len, side = "start", rc = F) {
  valid.sides <- c("start", "end")
  assertthat::assert_that(is(seqVector, "character") | is(seqVector, "list"), msg = "seqVector is not character or list")
  assertthat::assert_that(all(sapply(seqVector, is, "character")), msg = "Elements inside vector are not of character class")
  assertthat::assert_that(side %in% valid.sides, msg = "Side should be either 'start' or 'end")
  side <- if (side == "start") "right" else "left"
  trunc.vec <- sapply(seqVector, stringr::str_trunc, width = trunc_len, side = side, ellipsis = "")
  if (rc) trunc.vec <- dada2::rc(trunc.vec)
  return(trunc.vec)
}

MAIN <- "./"

#raw reads in reads/raw
fwd.reads <- list.files(path = file.path(MAIN, "reads", "raw"), pattern = "*_R1_001.fastq.gz", full.names = T, recursive = TRUE)
rev.reads <- list.files(path = file.path(MAIN, "reads", "raw"), pattern = "*_R2_001.fastq.gz", full.names = T, recursive = TRUE)

#cutadapt-filtered reads in reads/filtered
fwd.filt <- file.path(MAIN, "reads", "filtered", fwd.reads %>% sub("_000000000-J59L9_S.*_R1_001", "-R1", x = .) %>% basename)
rev.filt <- file.path(MAIN, "reads", "filtered", rev.reads %>% sub("_000000000-J59L9_S.*_R2_001", "-R2", x = .) %>% basename)

treatments <- fwd.filt %>% basename %>% strsplit("-") %>% sapply('[', 1)
samples <- fwd.filt %>% basename %>% strsplit("-R") %>% sapply('[', 1)

if (!identical(treatments, rev.filt %>% basename %>% strsplit("-") %>% sapply('[', 1))) stop()
if (!identical(samples, rev.filt %>% basename %>% strsplit("-R") %>% sapply('[', 1))) stop()

fwd.reads <- fwd.filt
rev.reads <- rev.filt

#final reads in reads/filterandtrim
fwd.filt <- file.path(MAIN, "reads", "filterandtrim", fwd.reads %>% basename)
rev.filt <- file.path(MAIN, "reads", "filterandtrim", rev.reads %>% basename)

#truncation length in cutadapt
trunc.len.f <- 265
trunc.len.r <- 190

filterAndTrim(fwd = fwd.reads,
              rev = rev.reads,
              filt = fwd.filt,
              filt.rev = rev.filt,
              truncLen = c(trunc.len.f, trunc.len.r),
              maxN = 0,
              maxEE = 4,
              truncQ = 2,
              multithread = 20,
              matchIDs = TRUE)

#get prior sequences
priors.f <- truncSeqs(sq_V3V4, trunc.len.f, "start", rc = F)
priors.r <- truncSeqs(sq_V3V4, trunc.len.r, "end", rc = T)

#control (uninoculated samples)
fwd.con <- fwd.filt %>% basename %>% startsWith("CTR") %>% which %>% '['(fwd.filt, .)
rev.con <- rev.filt %>% basename %>% startsWith("CTR") %>% which %>% '['(rev.filt, .)

iter.treat <- treatments[!treatments %in% c("CTR", "365")] %>% unique

#mergers list to contain results
mergers <- vector(mode = "list", length = length(iter.treat))
names(mergers) <- iter.treat

for (tr in seq_along(iter.treat)) {

  trt <- iter.treat[tr]

  #files for inoculated samples and controls
  fwd.trt <- c(fwd.filt %>% basename %>% startsWith(trt) %>% which %>% '['(fwd.filt, .),
               fwd.con)
  rev.trt <- c(rev.filt %>% basename %>% startsWith(trt) %>% which %>% '['(rev.filt, .),
               rev.con)

  #priors for Salinibacter, chloroplasts, mitochondria, and expected isolate V3V4
  priors.f.trt <- priors.f[c("SAL", "CLP", "MIT", trt)]
  priors.r.trt <- priors.r[c("SAL", "CLP", "MIT", trt)]

  out.files <- outer(c("fwd", "rev"), c(".err.RDS", ".drp.RDS", ".dd.RDS"), paste0) %>%
    as.vector %>%
    paste0(trt, ".", .) %>%
    file.path(MAIN, "dada", .)

  fwd.err <- learnErrors(fwd.trt, multithread = 30, verbose = 0)
  saveRDS(fwd.err, out.files[1])

  rev.err <- learnErrors(rev.trt, multithread = 30, verbose = 0)
  saveRDS(rev.err, out.files[2])

  fwd.drp <- derepFastq(fwd.trt)
  saveRDS(fwd.drp, out.files[3])

  rev.drp <- derepFastq(rev.trt)
  saveRDS(rev.drp, out.files[4])

  fwd.dd <- dada(fwd.drp, err = fwd.err, multithread = 30, pool = FALSE, priors = priors.f.trt)
  saveRDS(fwd.dd, out.files[5])

  rev.dd <- dada(rev.drp, err = rev.err, multithread = 30, pool = FALSE, priors = priors.r.trt)
  saveRDS(rev.dd, out.files[6])

  mergers[[tr]] <- mergePairs(fwd.dd, fwd.drp, rev.dd, rev.drp, trimOverhang = TRUE)

}

#Save mergers RDS
saveRDS(mergers, file.path(MAIN, "dada", "mergers.RDS"))

#create a sequence table for every treatment
mergers <- lapply(mergers, makeSequenceTable)
#modify the names for every table:
#uninoculated samples (controls) are the same per inoculation, and they all have
#the same name. Change that name to CTR_[isolate]
for (nm in names(mergers)) rownames(mergers[[nm]]) %<>% sub(x = ., "CTR", paste0("CTR_", nm))

#merge sequence tables
seqtab.raw <- mergeSequenceTables(tables = mergers, tryRC = T)
saveRDS(seqtab.raw, file.path(MAIN, "sequence.table.raw.RDS"))

#remove bimeras and save last
seqtab <- removeBimeraDenovo(seqtab.raw, method = "consensus", multithread = 30)
saveRDS(seqtab, file.path(MAIN, "sequence.table.clean.RDS"))



