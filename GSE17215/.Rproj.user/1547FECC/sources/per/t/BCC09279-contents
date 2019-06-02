fname='GSE76275_series_matrix.txt.gz' 
AnnotGPL=FALSE
destdir=tempdir()
getGPL=TRUE
parseCharacteristics=TRUE

library(readr )                       
dat <- read_lines(fname)
## get the number of !Series and !Sample lines
series_header_row_count <- sum(grepl("^!Series_", dat)) 
sample_header_start <- grep("^!Sample_", dat)[1]
samples_header_row_count <- sum(grepl("^!Sample_", dat))
series_table_begin_line = grep("^!series_matrix_table_begin", dat) 
##  colClasses <- c('character',rep('numeric',nrow(sampledat)))
datamat <- read_tsv(fname,quote='"',
                    na=c('NA','null','NULL','Null'), skip = series_table_begin_line,
                    comment = '!series_matrix_table_end')
datamat[1:4,1:4]