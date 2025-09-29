slices <- c(hansson_normal_count, common_unique_count)
lbls <- paste(c("Non-Brain-Enriched Proteins", "Brain-Enriched Proteins"), "\n", slices, sep="")
pie(slices, labels = lbls,
    main="", col=c("#a65e5e","#6e1b1b"))


slices <- c(sum(fragment_length_by_enrichment[1] == "Non-Brain-Enriched"), sum(fragment_length_by_enrichment[1] == "Brain-Enriched"))
lbls <- paste(c("Non-Brain-Enriched Fragments", "Brain-Enriched Fragments"), "\n", slices, sep="")
pie(slices, labels = lbls,
    main="", col=c("#a65e5e","#6e1b1b"))
