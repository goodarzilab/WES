library("XenofilteR")

args <- commandArgs(trailingOnly = TRUE)

align_mouse = args[1]
pattern_mouse = args[2] # *.mm10.srt

align_human = args[3]
pattern_human = args[4] # *.hg38.srt

xenofilterDir = args[5]
JOBS <- args[6]

bp.param <- SnowParam(workers = JOBS, type = "SOCK")

## make a dataframe/table with these two lists
# samplelist <- args[1]
# sample.list <- read.csv(samplelist, header = FALSE)

sample.list <- data.frame(
	human=list.files(pattern=paste0(pattern_human,'.bam$'),path=align_human,full.names = TRUE, recursive=T),
	mouse=list.files(pattern=paste0(pattern_mouse,'.bam$'),path=align_mouse,full.names = TRUE, recursive=T)
)

output.names <- gsub(paste0(pattern_human,'.bam$'), '.xenofilter.bam', sample.list$human)
output.names <- gsub(align_human, '/', output.names)

XenofilteR(
           sample.list, 
           destination.folder = xenofilterDir, 
           bp.param = bp.param, 
           output.names,
           MM_threshold = 12, 
           Unmapped_penalty = 4
)
