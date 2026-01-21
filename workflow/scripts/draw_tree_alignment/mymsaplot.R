#
getCols <- function (n) {
    col <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
             "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
             "#ccebc5", "#ffed6f")
    col2 <- c("#1f78b4", "#ffff33", "#c2a5cf", "#ff7f00", "#810f7c",
              "#a6cee3", "#006d2c", "#4d4d4d", "#8c510a", "#d73027",
              "#78c679", "#7f0000", "#41b6c4", "#e7298a", "#54278f")
    col3 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
              "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
              "#ffff99", "#b15928")
    grDevices::colorRampPalette(col3)(n)
}

mymsaplot <- function(p, fasta, offset = 0, width = 1, color = NULL, window = NULL, bg_line = TRUE, height = 0.8) {
    if (missingArg(fasta)) {
        x <- NULL
    } else if (is(fasta, "DNAbin") || is(fasta, "AAbin")) {
        x <- fasta
    } else if (is(fasta, "character")) {
        x <- treeio::read.fasta(fasta)
    } else if (is(fasta, "BStringSet")) {
        if (requireNamespace("Biostrings", quietly = TRUE) == TRUE) {
            temp_fasta <- tempfile("temp_fasta", fileext = ".fasta")
            Biostrings::writeXStringSet(fasta, temp_fasta)

            x <- treeio::read.fasta(temp_fasta)
        } else {
            stop("object is of class 'BStringSet' but library 'Biostrings' is not installed...\n-> please install 'Biostrings' from https://bioconductor.org/packages/Biostrings for handling objects of type 'BStringSet'.")
        }
    } else if (is(fasta, "DNAStringSet")) {
        x <- ape::as.DNAbin(fasta)
    } else if (is(fasta, "AAStringSet")) {
        x <- ape::as.AAbin(fasta)
    } else {
        x <- NULL
    }


    if (is.null(x) && is(p, "treedata") && length(p@tip_seq)) {
        x <- p@tip_seq
        p <- ggtree(p) + geom_tiplab()
    }

    if (is.null(x)) {
        stop("multiple sequence alignment is not available...\n-> check the parameter 'fasta'...")
    }

    x <- as.matrix(x)

    if (!all(labels(x) %in% p$data$label)) {
        stop("taxa name in input sequences are not match with the ones on the tree, please check your input files...")
    }

    if (is.null(window)) {
        window <- c(1, ncol(x))
    }

    slice <- seq(window[1], window[2], by = 1)
    x <- x[, slice]

    seqs <- lapply(1:nrow(x), function(i) {
        seq <- as.vector(as.character(x[i, ]))
        seq[seq == "?"] <- "-"
        seq[seq == "*"] <- "-"
        seq[seq == " "] <- "-"
        return(seq)
    })

    names(seqs) <- labels(x)

    if (is.null(color)) {
        alphabet <- unlist(seqs) %>% unique()
        alphabet <- alphabet[alphabet != "-"]
        ## color <- rainbow_hcl(length(alphabet))
        color <- getCols(length(alphabet))
        names(color) <- alphabet
        color <- c(color, "-" = NA)
    }

    df <- p$data
    ## if (is.null(width)) {
    ##     width <- (df$x %>% range %>% diff)/500
    ## }

    ## convert width to width of each cell
    width <- width * (df$x %>% range() %>% diff()) / diff(window)

    df <- df[df$isTip, ]
    start <- max(df$x) * 1.02 + offset

    seqs <- seqs[df$label[order(df$y)]]
    ## seqs.df <- do.call("rbind", seqs)

    h <- ceiling(diff(range(df$y)) / length(df$y))
    xmax <- start + seq_along(slice) * width
    xmin <- xmax - width
    y <- sort(df$y)
    ymin <- y - height / 2 * h
    ymax <- y + height / 2 * h

    from <- to <- NULL

    lines.df <- data.frame(from = min(xmin), to = max(xmax), y = y)

    if (bg_line) {
        p <- p + geom_segment(
            data = lines.df, aes(x = from, xend = to, y = y, yend = y),
            size = h * .2, inherit.aes = FALSE, linetype = "dotted", color = "grey90"
        )
    }

    msa <- lapply(1:length(y), function(i) {
        data.frame(
            name = names(seqs)[i],
            xmin = xmin,
            xmax = xmax,
            ymin = ymin[i],
            ymax = ymax[i],
            seq = seqs[[i]]
        )
    })

    msa.df <- (
        do.call("rbind", msa)
        |> filter(seq != "-")
    )

    p <- p + geom_rect(
        aes(
            xmin = xmin, xmax = xmax,
            ymin = ymin, ymax = ymax,
            fill = seq
        ),
        data = msa.df, inherit.aes = FALSE
    ) +
        scale_fill_manual(values = color)

    breaks <- graphics::hist(seq_along(slice), breaks = 10, plot = FALSE)$breaks
    pos <- start + breaks * width
    data_axis <- data.frame(from = breaks + 1, to = pos)
    attr(p, "data_axis") <- data_axis

    return(p)
}
#gheatmap_colors_file <- "imports/maf_evolution/maf_cnc_colors.txt"
#family_colors_file <- "imports/maf_evolution/maf_cnc_family_colors.txt"
#group_colors_file <- "imports/maf_evolution/maf_cnc_group_colors.txt"
#class_colors_file <- "imports/maf_evolution/maf_cnc_class_colors.txt"

read_validate_fasta <- function(msa_fasta_file, tree_nwk) {
    x <- treeio::read.fasta(msa_fasta_file)
    if (is.null(x)) {
        stop("multiple sequence alignment is not available...\n-> check the parameter 'fasta'...")
    }
    print(head(x))

    # it seems treeio::read.fasta keep all header information as sequence name
    # we need to specifically remove everything after the first space character
    names(x) <- gsub("\\s.*", "", names(x))
    print(head(x))

    found <- names(x) %in% tree_nwk$tip.label

    print("MSA sequence names:")
    print(names(x))
    print("Tree TIPs:")
    print(tree_nwk$tip.label)

    if (!all(found)) {
        print(names(x)[!found])
        stop("Above taxa names in input sequences do not match with the ones on the tree, please check your input files...")
    }

    x
}

