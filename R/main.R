#' One-hot encode DNA sequences into nucleotide position factors
#'
#' Converts a character vector of fixed-length DNA strings into a data frame
#' where each column corresponds to a nucleotide position (e.g. \code{nt_pos1}
#' to \code{nt_pos5}). Each position is stored as a factor with levels
#' \code{c("A", "T", "C", "G")}, which is suitable for use in machine learning
#' models such as random forests.
#'
#' @param dna_strings A character vector of DNA sequences of equal length.
#'   Each string should contain only the characters \code{"A"}, \code{"T"},
#'   \code{"C"}, and \code{"G"}.
#'
#' @return A data frame with one column per nucleotide position
#'   (\code{nt_pos1}, \code{nt_pos2}, ...). Each column is a factor with
#'   levels \code{c("A", "T", "C", "G")}. The number of rows equals
#'   \code{length(dna_strings)}.
#'
#' @examples
#' dna_encoding(c("ATCGA", "GGACA"))
#' @export
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A probability and status for multiple observations
#'
#' Uses a fitted machine learning model (typically a random forest) to predict
#' m6A site probabilities and binary status labels for a data frame of input
#' features. The function standardises factor levels for categorical variables,
#' encodes the \code{DNA_5mer} feature into nucleotide position factors, and
#' appends prediction results as new columns.
#'
#' @import randomForest
#' @importFrom stats predict
#' @param ml_fit A fitted classification model that supports
#'   \code{predict(ml_fit, newdata, type = "prob")}, e.g. an object returned by
#'   \code{randomForest::randomForest()} trained on the m6A data.
#' @param feature_df A data frame containing at least the following columns:
#'   \itemize{
#'     \item \code{gc_content} (numeric)
#'     \item \code{RNA_type} (character or factor; one of
#'       \code{"mRNA"}, \code{"lincRNA"}, \code{"lncRNA"}, \code{"pseudogene"})
#'     \item \code{RNA_region} (character or factor; one of
#'       \code{"CDS"}, \code{"intron"}, \code{"3'UTR"}, \code{"5'UTR"})
#'     \item \code{exon_length} (numeric)
#'     \item \code{distance_to_junction} (numeric)
#'     \item \code{evolutionary_conservation} (numeric)
#'     \item \code{DNA_5mer} (character; DNA string of length 5)
#'   }
#' @param positive_threshold A numeric value between 0 and 1 specifying the
#'   probability cutoff above which an observation is labelled
#'   \code{"Positive"}. Defaults to \code{0.5}.
#'
#' @return A data frame equal to \code{feature_df} with two additional columns:
#'   \describe{
#'     \item{\code{predicted_m6A_prob}}{Numeric vector of predicted probabilities
#'       for the positive class.}
#'     \item{\code{predicted_m6A_status}}{Character vector with values
#'       \code{"Positive"} or \code{"Negative"}, based on
#'       \code{positive_threshold}.}
#'   }
#'
#' @examples
#' rf_fit <- readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#' example_data <- read.csv(system.file("extdata", "m6A_input_example.csv", package="m6APrediction"))
#' preds <- prediction_multiple(rf_fit, example_data, positive_threshold = 0.5)
#' head(preds)
#'
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) #Check errors if incorrect column names of input data.frame

  # 设定因子水平顺序
  feature_df$RNA_type <- factor(
    feature_df$RNA_type,
    levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")
  )
  feature_df$RNA_region <- factor(
    feature_df$RNA_region,
    levels = c("CDS", "intron", "3'UTR", "5'UTR")
  )

  # 将 DNA_5mer one-hot 编码为 nt_pos1-5（因子 levels 固定为 A/T/C/G）
  encoded_df <- cbind(feature_df, dna_encoding(feature_df$DNA_5mer))

  # 预测 Positive 类别的概率
  prob_mat <- predict(ml_fit, newdata = encoded_df, type = "prob")
  pred_prob <- prob_mat[, "Positive"]

  # 添加预测列
  feature_df$predicted_m6A_prob <- pred_prob
  feature_df$predicted_m6A_status <- ifelse(
    pred_prob > positive_threshold,
    "Positive",
    "Negative"
  )

  return(feature_df) #return a data.frame with supplied columns of predicted m6A prob and predicted m6A status
}

#' Predict m6A probability and status for a single observation
#'
#' Convenience wrapper around \code{\link{prediction_multiple}} that takes the
#' feature values for a single site as separate arguments, constructs the input
#' data frame, and returns a concise prediction vector.
#'
#' @param ml_fit A fitted classification model compatible with
#'   \code{prediction_multiple()}, typically a random forest trained on the m6A
#'   dataset.
#' @param gc_content Numeric GC content of the site (between 0 and 1).
#' @param RNA_type Character string specifying the RNA type. Must be one of
#'   \code{"mRNA"}, \code{"lincRNA"}, \code{"lncRNA"}, or \code{"pseudogene"}.
#' @param RNA_region Character string specifying the RNA region. Must be one of
#'   \code{"CDS"}, \code{"intron"}, \code{"3'UTR"}, or \code{"5'UTR"}.
#' @param exon_length Numeric exon length (e.g. in nucleotides).
#' @param distance_to_junction Numeric distance to the nearest exon–intron
#'   junction.
#' @param evolutionary_conservation Numeric conservation score (typically
#'   between 0 and 1).
#' @param DNA_5mer Character DNA sequence of length 5, using only
#'   \code{"A"}, \code{"T"}, \code{"C"}, and \code{"G"}.
#' @param positive_threshold Numeric probability cutoff between 0 and 1 used to
#'   determine the \code{"Positive"} vs \code{"Negative"} label. Defaults to
#'   \code{0.5}.
#'
#' @return A named vector of length 2 containing:
#'   \describe{
#'     \item{\code{"predicted_m6A_prob"}}{Numeric predicted probability for the
#'       positive class.}
#'     \item{\code{"predicted_m6A_status"}}{Character label,
#'       \code{"Positive"} or \code{"Negative"}.}
#'   }
#'
#' @examples
#' rf_fit <- readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#' example_data <- read.csv(system.file("extdata", "m6A_input_example.csv", package="m6APrediction"))
#' prediction_single(
#'   ml_fit = rf_fit,
#'   gc_content = 0.5,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 10,
#'   distance_to_junction = 8,
#'   evolutionary_conservation = 0.5,
#'   DNA_5mer = "GGACA",
#'   positive_threshold = 0.5
#' )
#'
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  # 把单个观测值组装成 data.frame
  feature_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )

  # 使用上面写好的 prediction_multiple() 来做预测（内部会处理因子水平和编码）
  pred_df <- prediction_multiple(
    ml_fit = ml_fit,
    feature_df = feature_df,
    positive_threshold = positive_threshold
  )

  # 抽取概率和状态，组成长度为 2 的命名向量
  returned_vector <- c(
    predicted_m6A_prob   = pred_df$predicted_m6A_prob[1],
    predicted_m6A_status = as.character(pred_df$predicted_m6A_status[1])
  )

  return(returned_vector) #return a named vector with values for predicted m6A prob and predicted m6A status
}
