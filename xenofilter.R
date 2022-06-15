library("XenofilteR")

args <- commandArgs(trailingOnly = TRUE)

samplelist <- args[1]
inputFolder <- args[2]
JOBS <- args[3]
align_mouse = args[4]
align_human = args[5]

# make a dataframe/table with these two lists
list.files(pattern='*.mm10.srt.bam$',path=align_mouse,full.names = TRUE, recursive=T)

list.files(pattern='*.hg38.srt.bam$',path=align_human,full.names = TRUE, recursive=T)

# figure out how to create a table with two lists

bp.param <- SnowParam(workers = JOBS, type = "SOCK")

sample.list <- read.csv(samplelist, header = FALSE)

output.names <- gsub('.human.srt.bam', '.processed.bam', sample.list$'V1')

output.names <- gsub(inputFolder, '', output.names)

XenofilteR(sample.list, destination.folder = "./", bp.param = bp.param, output.names,
           MM_threshold = 12, Unmapped_penalty = 4


