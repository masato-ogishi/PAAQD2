#' PAAQD2: a code-refined version of PAAQD
#' @param peptideSet A set of peptide sequences.
#' @import randomForest
#' @export
#' @rdname PAAQD2
#' @name PAAQD2
PAAQD2 <- function(peptideSet){
  EPICPATH = system.file("extdata", package="PAAQD2")
  WA9 = paste(EPICPATH, "/ALLELE_WA9_rang1_100.txt", sep = "")
  WC = paste(EPICPATH, "/HLA_A_WC.txt", sep = "")
  HLAseq = paste(EPICPATH, "/ALLELE_seq.txt", sep = "")
  INSC = list.files(paste(EPICPATH, "/AAindex_norm/", sep = ""))
  INSC.dex = c(1:40)
  INSC = INSC[INSC.dex]
  WA9 = read.table(WA9, sep = "\t", header = F)
  WC = read.table(WC, sep = "\t", header = F)
  HLAseq = read.table(HLAseq, sep = "\t", header = F)
  HLAseq = substr(HLAseq[1, 1], 1, 300)
  HLAindex = c()
  indexer <- function(aachar){
    if (aachar == "A") {
      return(1)
    }
    if (aachar == "R") {
      return(2)
    }
    if (aachar == "N") {
      return(3)
    }
    if (aachar == "D") {
      return(4)
    }
    if (aachar == "C") {
      return(5)
    }
    if (aachar == "Q") {
      return(6)
    }
    if (aachar == "E") {
      return(7)
    }
    if (aachar == "G") {
      return(8)
    }
    if (aachar == "H") {
      return(9)
    }
    if (aachar == "I") {
      return(10)
    }
    if (aachar == "L") {
      return(11)
    }
    if (aachar == "K") {
      return(12)
    }
    if (aachar == "M") {
      return(13)
    }
    if (aachar == "F") {
      return(14)
    }
    if (aachar == "P") {
      return(15)
    }
    if (aachar == "S") {
      return(16)
    }
    if (aachar == "T") {
      return(17)
    }
    if (aachar == "W") {
      return(18)
    }
    if (aachar == "Y") {
      return(19)
    }
    if (aachar == "V") {
      return(20)
    }
  }
  for (i in 1:300) {
    SeqDex = indexer(substr(HLAseq, i, i))
    HLAindex = c(HLAindex, SeqDex)
  }
  if (c(".") %in% unlist(strsplit(peptideSet[1], ""))) {
    PepsTable = readLines(peptideSet)
  } else {
    PepsTable = peptideSet
  }
  TranPepIndex <- function(PepsTable) {
    Peps = PepsTable
    n_Peps = length(Peps)
    TranPeps = matrix(0, n_Peps, 9)
    for (i in 1:n_Peps) {
      for (j in 1:9) {
        TranPeps[i, j] = indexer(substr(Peps[i], j, j))
      }
    }
    return(TranPeps)
  }
  SCal <- function(TranPeps) {
    n_Peps = nrow(TranPeps)
    n_WC = nrow(WC)
    ScALL = c()
    for (i in 1:n_Peps) {
      ScRow = c()
      for (j in 1:9) {
        subSc = 0
        for (k in 1:n_WC) {
          subSc = subSc + InSc[TranPeps[i, j], HLAindex[WC[k, 10]]] * WC[k, j]
        }
        Score = WA9[TranPeps[i, j], j] * (subSc/sum(WC[, j] != 0))
        ScRow = c(ScRow, Score)
      }
      ScALL = rbind(ScALL, ScRow)
    }
    colnames(ScALL) = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")
    return(ScALL)
  }

  # QuantumFea
  TNonaPeps = TranPepIndex(PepsTable)
  Quandes = read.table(paste(EPICPATH, "/Quantums.txt", sep = ""), header = T, sep = "\t")
  QuantumRES = c()
  for (i in 1:nrow(TNonaPeps)) {
    Row.QuantumRES = c()
    for (j in 1:9) {
      Each.QuantumRES = Quandes[TNonaPeps[i, j], ]
      Row.QuantumRES = c(Row.QuantumRES, as.numeric(unlist(Each.QuantumRES)))
    }
    QuantumRES = rbind(QuantumRES, Row.QuantumRES)
  }
  Qnames = readLines(paste(EPICPATH, "/Qnames.txt", sep = ""))
  QuantumRES = as.matrix(QuantumRES)
  colnames(QuantumRES) = Qnames
  Quan.data <- QuantumRES

  # Transformer
  Peps_Tr = c()
  for (i in 1:length(INSC)) {
    InSc = read.table(paste(EPICPATH, "/AAindex_norm/",
                            INSC[i], sep = ""), sep = "\t", header = F)
    EACH_Peps_Tr = SCal(TranPepIndex(PepsTable))
    Peps_Tr = cbind(Peps_Tr, EACH_Peps_Tr)
  }
  NaNames = c()
  for (l in 1:length(INSC)) {
    for (k in 1:9) {
      EACH_name = paste("Pos", k, "AA", INSC.dex[l], sep = "")
      NaNames = c(NaNames, EACH_name)
    }
  }
  colnames(Peps_Tr) = NaNames
  Trs.data <- Peps_Tr

  # PAAQD analysis
  Trs.data = cbind(Trs.data, Quan.data)
  if (c(".") %in% unlist(strsplit(peptideSet[1], ""))) {
    RAW.seq = readLines(peptideSet)
  } else {
    RAW.seq = peptideSet
  }
  rownames(Trs.data) = c()
  pred = as.character(predict(model, Trs.data))
  pred[pred=="Y"] <- "Positive"
  pred[pred=="N"] <- "Negative"
  predPROB = predict(model, Trs.data, type="prob")[,"Y"]
  RES = data.frame("Peptide"=RAW.seq,
                   "PredictedImmunogenicity"=factor(pred, levels=c("Positive", "Negative")),
                   "Probability"=predPROB,
                   stringsAsFactors = F)
  return(RES)
}
