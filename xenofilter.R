suppressMessages(suppressWarnings(library (XenofilteR)))
message ('XenofilteR loaded!')

args <- commandArgs(trailingOnly = TRUE)

align_human = args[1]
pattern_human = args[2] # *.hg38.srt
align_mouse = args[3]
pattern_mouse = args[4] # *.mm10.srt

JOBS <- args[5]

#JOBS = 10
#align_human = "test"#"align_human"
#align_mouse = "test"#"align_mouse"
#pattern_human = ".hg38.srt"
#pattern_mouse = ".mm10.srt"

bp.param <- SnowParam(workers = JOBS, type = "SOCK")

## make a dataframe/table with these two lists
# samplelist <- args[1]
# sample.list <- read.csv(samplelist, header = FALSE)

sample.list <- data.frame(
	human=list.files(pattern=paste0(pattern_human,'.bam$'),path=align_human,full.names = TRUE, recursive=T),
	mouse=list.files(pattern=paste0(pattern_mouse,'.bam$'),path=align_mouse,full.names = TRUE, recursive=T)
)

output.names <- gsub(paste0(pattern_human,'.bam$'), '.bam', sample.list$human)
output.names <- gsub(align_human, '/', output.names)

XenofilteR(
           sample.list, 
           destination.folder = './', 
           bp.param = bp.param, 
           output.names,
           MM_threshold = 12, 
           Unmapped_penalty = 4
)
