# check if human and mosue symbols share the same homolog id
hm_homolog = readr::read_tsv("./metadata/HOM_MouseHumanSequence.rpt.txt")
hm_homolog_hg = hm_homolog[hm_homolog$`Common Organism Name` == "human", ]
hm_homolog_mm = hm_homolog[hm_homolog$`Common Organism Name` == "mouse, laboratory", ]

hm_homolog_sp = split(hm_homolog, hm_homolog$`HomoloGene ID`)
hm_sel = sapply(hm_homolog_sp, function(x) all(c("mouse, laboratory", "human") %in% x$`Common Organism Name`))
hm_homolog_sp = hm_homolog_sp[hm_sel]
# unique human symbols
human_symbols = sapply(hm_homolog_sp, function(x) {
        symbol = x[x$`Common Organism Name` == "human", ]$Symbol
        if (length(symbol) == 1) {
                return(symbol)
        } else {
                return(NA)
        }
})
# unique mouse symbols
mouse_symbols = sapply(hm_homolog_sp, function(x) {
        symbol = x[x$`Common Organism Name` == "mouse, laboratory", ]$Symbol
        if (length(symbol) == 1) {
                return(symbol)
        } else {
                return(NA)
        }
})
common_indices = which(!is.na(human_symbols) & !is.na(mouse_symbols))
homolog_id = names(human_symbols)[common_indices]
human_symbols = human_symbols[common_indices]
mouse_symbols = mouse_symbols[common_indices]
