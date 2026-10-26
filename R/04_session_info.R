#!/usr/bin/env Rscript

sink(file.path(VALIDATION_DIR, "sessionInfo.txt"))
print(sessionInfo())
sink()
cat("Session information recorded.\n")

