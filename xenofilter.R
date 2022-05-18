library("XenofilteR")

args <- commandArgs(trailingOnly = TRUE)

samplelist <- args[1]
inputFolder <- args[2]
JOBS <- args[3]

bp.param <- SnowParam(workers = JOBS, type = "SOCK")

sample.list <- read.csv(samplelist, header = FALSE)

output.names <- gsub('.human.srt.bam', '.processed.bam', sample.list$'V1')

output.names <- gsub(inputFolder, '', output.names)

XenofilteR(sample.list, destination.folder = "./", bp.param = bp.param, output.names,
           MM_threshold = 12, Unmapped_penalty = 4


