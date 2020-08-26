library(ggplot2)
library(unpakathon)
library(dplyr)
library(plyr)
library(readr)
snps = read.csv('./snps_AT1G52990.1.csv')
#head(snps)


#Open all csv files in snps_downloads folder
dir2 = "snps_downloads"
all_snps = list.files(path = dir2,
                      pattern = "*.csv",
                      full.names = TRUE)
all_snp_data = ldply(all_snps, read_csv) #Combine all csv files in low_mutation_snps folder

#Open all csv files in ins_downloads folder
dir5 = "ins_downloads"
all_ins = list.files(path = dir5,
                     pattern = "*.csv",
                     full.names = TRUE)
all_ins_data = ldply(all_ins, read_csv)

#Open all csv files in dels_downloads folder
dir6 = "dels_downloads"
all_dels = list.files(path = dir6,
                      pattern = "*.csv",
                      full.names = TRUE)
all_dels_data = ldply(all_dels, read_csv)

#number of SNP occurrences in AT numbers w/ high mutation level
SNPtotal_hi <-
  aggregate(strain ~ pos , hi_data, function(x)
    length(unique(x)))
nrow(SNPtotal_hi) # number of unique positions in total

#number of SNP occurrences in AT numbers w/ low mutation level
SNPtotal_low <-
  aggregate(strain ~ pos , all_snp_data, function(x)
    length(unique(x)))
nrow(SNPtotal_low) # number of unique positions in total

# Get merged data from snps and indels with unpakathon data
# get_data("all") to get all data in folder. get_data("AT Number Here") without .1 for specific AT Num.
get_data <- function(AT_num) {
  if (AT_num == "all_joined") {
    gene_data <- get_data("all")
    joined_data <-
      gene_data %>% left_join(
        phenolong %>% dplyr::filter(variable == "fruitnum") %>%
          dplyr::group_by(locus, variable) %>%
          dplyr::summarize(value = mean(value, na.rm = T)) %>% pivot_wider(names_from = "variable")
      )
    return(joined_data)
  }
  if (AT_num != "all") {
    getall_df <- get_data("all")
    a <- getall_df %>% filter(locus == AT_num)
    a
  } else{
    df_total = data.frame()
    for (i in 1:length(all_snps)) {
      all_data <- read_csv(all_snps[i])
      snpnum = 0
      tryCatch(
        expr = {
          all_lo_snps <-
            aggregate(strain ~ pos , all_data, function(x)
              length(unique(x)))
          snpnum = c(nrow(all_lo_snps))
        },
        error = function(e) {
          snpnum = 0
        },
        warning = function(w) {
          snpnum = 0
        }
      )
      df1 <- data.frame(AT_Number_path = c(all_snps[i]),
                        SNP_Number = c(snpnum))
      df_total <- rbind(df_total, df1)
      
    }
    tmp <- df_total
    shortened_AT <-
      gsub(".csv.*$", "", gsub("^.*snps_", "", tmp$AT_Number))
    dfreturn = tmp %>% mutate(locus = shortened_AT) %>% relocate(locus, .before = SNP_Number)
    # get indels at all_ins
    df_tmp_indels = data.frame()
    for (j in 1:length(dfreturn$locus)) {
      h <-
        grepl(paste("ins_", dfreturn$locus[j], ".csv", sep = ""),
              all_ins,
              fixed = TRUE)
      index_ins <- min(which(h == TRUE))#integer
      f <-
        grepl(paste("dels_", dfreturn$locus[j], ".csv", sep = ""),
              all_dels,
              fixed = TRUE)
      index_dels <- min(which(f == TRUE))#integer
      if (index_ins != Inf) {
        all_ins_df <- read_csv(all_ins[index_ins])
        insnum = 0
        tryCatch(
          expr = {
            instotal_lo <-
              aggregate(strain ~ pos , all_ins_df, function(x)
                length(unique(x)))
            insnum <- c(nrow(instotal_lo))
          },
          error = function(e) {
            insnum = 0
          },
          warning = function(w) {
            insnum = 0
          }
        )
      } else {
        insnum = 0
      }
      if (index_dels != Inf) {
        all_dels_df <- read_csv(all_dels[index_dels])
        delsnum = 0
        tryCatch(
          expr = {
            delstotal_lo <-
              aggregate(strain ~ pos , all_dels_df, function(x)
                length(unique(x)))
            delsnum <- c(nrow(delstotal_lo))
          },
          error = function(e) {
            delsnum = 0
          },
          warning = function(w) {
            delsnum = 0
          }
        )
      } else {
        delsnum = 0
      }
      df2 <- data.frame(Indels = c(insnum + delsnum))
      df_tmp_indels <- rbind(df_tmp_indels, df2)
    }
    final_return = dfreturn %>% mutate(Indels = df_tmp_indels$Indels)
    final_return$locus <- sub("\\.\\d+$", "", final_return$locus)
    return(final_return)
  }
}
get_data("all")

test_save <- as.data.frame(get_data("all"))
test_save
df_sv <- data.frame(
  locus = c(test_save$locus),
  SNP_Number = c(test_save$SNP_Number),
  Indels = c(test_save$Indels)
)
#write.csv(df_sv,"my_csv/SNPs_and_Indels.csv", row.names = TRUE)

sorted_SNPs <- df_sv[order(df_sv$SNP_Number), ]
bottom_20_SNPs <- head(sorted_SNPs, n = 117)
top_20_SNPs <- tail(sorted_SNPs, n = 117)

sorted_Indels <- df_sv[order(df_sv$Indels), ]
bottom_20_Indels <- head(sorted_Indels, n = 117)
top_20_Indels <- tail(sorted_Indels, n = 117)

common_top <-
  inner_join(top_20_Indels, top_20_SNPs) # AT numbers that were at the top 20% of Indels and SNPs
common_bottom <-
  inner_join(bottom_20_Indels, bottom_20_SNPs) # AT numbers that were at the bottom 20% of Indels and SNPs
common_top
common_bottom

common_top_joined <-
  common_top %>% left_join(
    phenolong %>% dplyr::filter(variable == "fruitnum") %>%
      dplyr::group_by(locus, variable) %>%
      dplyr::summarize(value = mean(value, na.rm = T)) %>% pivot_wider(names_from =
                                                                         "variable")
  )
common_bottom_joined <-
  common_bottom %>% left_join(
    phenolong %>% dplyr::filter(variable == "fruitnum") %>%
      dplyr::group_by(locus, variable) %>%
      dplyr::summarize(value =
                         mean(value, na.rm = T)) %>% pivot_wider(names_from = "variable")
  )

common_top_joined$quantile = "top"
common_bottom_joined$quantile = "bottom"
topbottom = rbind(common_top_joined, common_bottom_joined)
summary(aov(fruitnum ~ quantile, data = topbottom))
boxplot(fruitnum ~ quantile, data = topbottom)
joined_dataframe <- get_data("all_joined")
joined_dataframe

plot(joined_data$fruitnum ~ joined_data$SNP_Number)
summary(lm(joined_data$fruitnum ~ joined_data$Indels))
summary(lm(
  log(joined_dataframe$fruitnum + 1) ~ log(joined_dataframe$SNP_Number + 1)
))
abline(lm(
  log(joined_dataframe$fruitnum + 1) ~ log(joined_dataframe$SNP_Number + 1)
))

ggplot(data = get_data("all"), aes(x = locus, y = SNP_Number)) +
  geom_bar(stat = "identity")

plotsnps <- get_data("all")
ggplot(test_save, aes(x = Indels)) + geom_histogram(binwidth = 10, color = "black", fill = "darkblue")
ggplot(test_save, aes(x = SNP_Number)) + geom_histogram(binwidth = 100, color = "black", fill = "darkblue")
plotsnps

ggplot(plotsnps$Indels, aes(x = Indels)) + geom_histogram(binwidth = 5,
                                                          color = "black",
                                                          fill = "green")
ggplot(plotsnps, aes(x = SNP_Number)) + geom_histogram(binwidth = 100,
                                                       color = "black",
                                                       fill = "blue")
#get random AT Numbers from Unpakathon
rand_1000 <-
  data.frame(Rand_ATNum = c(sample(phenowide$locus, 1000)))
rand_1000[rand_1000 == ""] <- NA
rand_1000 <- na.omit(rand_1000)
rand_1000$Rand_ATNum[length(rand_1000$Rand_ATNum)]
length(rand_1000$Rand_ATNum)

twobytwo_df = data.frame(
  AT_Number = c(
    "AT1G52990",
    "AT3G59400",
    "AT2G20430",
    "AT1G62300",
    "AT2G43240",
    "AT2G32510",
    "AT1G24020",
    "AT1G64561",
    "AT4G12780"
    ,
    "AT5G18940",
    "AT5G09940",
    "AT3G55830",
    "AT2G20100",
    "AT3G20260",
    "AT1G17500",
    "AT2G45750",
    "AT4G06686",
    "AT2G37860",
    "AT3G47710",
    "AT3G32168",
    "AT5G32345"
  )
)
twobytwo_df$AT_Number

dl_AT_Num <- function(z) {
  for (i in 1:length(twobytwo_df$AT_Number)) {
    #download snps
    destfile_snps <- "twobytwo_snps/snps_ATXXXXXXX.1.csv"
    snps_url <-
      "https://tools.1001genomes.org/api/v1.1/variants.csv?gid=ATXXXXXXX.1;type=snps;accs=all"
    snps_url_cust <-
      gsub("ATXXXXXXX", twobytwo_df$AT_Number[i], snps_url)
    destfile_cust_snps <-
      gsub("ATXXXXXXX", twobytwo_df$AT_Number[i], destfile_snps)
    try(dl_snps <- download.file(snps_url_cust, destfile_cust_snps))
    dl_snps
    
    #download ins
    destfile_ins <- "twobytwo_ins/ins_ATXXXXXXX.1.csv"
    ins_url <-
      "https://tools.1001genomes.org/api/v1.1/variants.csv?gid=ATXXXXXXX.1;type=ins;accs=all"
    ins_url_cust <-
      gsub("ATXXXXXXX", twobytwo_df$AT_Number[i], ins_url)
    destfile_cust_ins <-
      gsub("ATXXXXXXX", twobytwo_df$AT_Number[i], destfile_ins)
    try(dl_ins <- download.file(ins_url_cust, destfile_cust_ins))
    dl_ins
    
    #download dels
    destfile_dels <- "twobytwo_dels/dels_ATXXXXXXX.1.csv"
    dels_url <-
      "https://tools.1001genomes.org/api/v1.1/variants.csv?gid=ATXXXXXXX.1;type=dels;accs=all"
    dels_url_cust <-
      gsub("ATXXXXXXX", twobytwo_df$AT_Number[i], dels_url)
    destfile_cust_dels <-
      gsub("ATXXXXXXX", twobytwo_df$AT_Number[i], destfile_dels)
    try(dl_dels <- download.file(dels_url_cust, destfile_cust_dels))
    dl_dels
    
  }
}
dl_AT_Num("")