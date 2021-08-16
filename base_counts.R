args <- commandArgs(trailingOnly = TRUE)


dat <- data.table::fread(paste0(args[1], "/", args[2]), 
                         header=FALSE) 

dat <- magrittr::set_colnames(dat, c("chr","chr.len"))


dat <- dplyr::mutate(dat, global.end = cumsum(as.numeric(chr.len)),
                     global.start = global.end - chr.len)

dat <- dplyr::relocate(dat, global.start, .before = global.end)

write.table(dat,
            file = paste0(args[1], "/", args[3]),
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")
