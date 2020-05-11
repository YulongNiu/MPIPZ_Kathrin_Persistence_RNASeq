
## originally by Yulong Niu
## yulong.niu@hotmail.com

##############################merge/rename fq###########################
library('foreach')
library('doParallel')
library('tidyverse')

rawfqPath <- '/biodata/dep_psl/grp_rgo/ljsphere/atlj_cros_rnaseq'
resFolder <- '/netscratch/dep_psl/grp_rgo/yniu/KathrinPersistence/raw_data'
catPath <- '/bin/cat'

ncore <- 40

anno <- read_delim('/netscratch/dep_psl/grp_rgo/yniu/KathrinPersistence/results/mapping.txt', delim = '\t') %>%
  mutate(samplePrefix = str_replace(library_number, '\\.', '_'))

rawfq <- dir(rawfqPath,
             pattern = 'fastq.gz')

fqs <- rawfq %>%
  strsplit('_', fixed = TRUE) %>%
  lapply('[', c(1, 2, 7)) %>%
  sapply(paste, collapse = '_')

fqIdx <- split(seq_along(fqs), fqs)
fqPrefix <- names(fqIdx)

registerDoParallel(cores = ncore)
fqmd5 <- foreach (i = seq_along(fqIdx)) %dopar% {
  eachmd5 <- fqIdx[[i]] %>%
    {file.path(rawfqPath, rawfq[.])} %>%
    paste('md5sum ', .) %>%
    system(intern = TRUE) %>%
    strsplit(split = ' ', fixed = TRUE) %>%
    unlist %>%
    .[c(1, 3)]

  return(eachmd5)
} %>% do.call(rbind, .)

(fqmd5[, 1] %>% unique %>% length) == nrow(fqmd5)
stopImplicitCluster()

## merge
registerDoParallel(cores = ncore)
foreach (i = seq_along(fqIdx), .combine = c) %dopar% {
  ## input files
  fqin <- fqIdx[[i]] %>%
    {file.path(rawfqPath, rawfq[.])} %>%
    paste(collapse = ' ')

  fqout <- fqPrefix[i] %>%
    paste0('.fq.gz') %>%
    {file.path(resFolder, .)}

  mergeC <- paste(catPath,
                  fqin,
                  '>',
                  fqout)

  print(mergeC)

  system(mergeC)

  return(NULL)
}
stopImplicitCluster()


## change file names
for (i in seq_len(nrow(anno))) {

  ## merge different batch
  fqin <- anno[i, 9] %>%
    as.character %>%
    .[!is.na(.)] %>%
    file.path(resFolder, .) %>%
    paste0('_R1.fq.gz') %>%
    paste(collapse = ' ')

  fqout <- anno[i, 2] %>%
    file.path(resFolder, .) %>%
    paste0('_R1.fq.gz')

  mergeC <- paste(catPath,
                  fqin,
                  '>',
                  fqout)
  print(mergeC)

  system(mergeC)

  fqin <- anno[i, 9] %>%
    as.character %>%
    .[!is.na(.)] %>%
    file.path(resFolder, .) %>%
    paste0('_R2.fq.gz') %>%
    paste(collapse = ' ')

  fqout <- anno[i, 2] %>%
    file.path(resFolder, .) %>%
    paste0('_R2.fq.gz')

  mergeC <- paste(catPath,
                  fqin,
                  '>',
                  fqout)
  print(mergeC)

  system(mergeC)
}

system('rm `ls | grep "4276\\|4316"`')
########################################################################
