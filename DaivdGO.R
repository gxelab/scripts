
#' parse David annotation clustering
#' @param x txt file path of annotCluster results
#' @return GO annotation clustering: term + cluster + cluster enrichment score
#' @example:
#'  x <- '~/tmp/t2t_B3AAF42FA92B1461018944337.txt'
#'  x` <- parseAnnotCluster(x, reduce = T)
#'
parseAnnotCluster <- function(x, short = T, reduce = F, ben = 0.01){
    # read input
    annot <- readLines(x)
    # split clusters
    annot.cl <- mapply(
        function(i, j){
            tc <- textConnection(paste(annot[(i + 1):(j - 1)], collapse = '\n'))
            res <- read.table(tc, sep = '\t', header = T, quote = '', stringsAsFactors = F)
            close(tc)
            colnames(res)[4] <- 'Percent'
            ac <- strsplit(annot[i], '\\s')[[1]]
            res <- cbind.data.frame(
                Cluster = ac[3], Enrichment.Score = ac[6], res, stringsAsFactors = F)
        },
        grep('^Annotation', annot), grep('^$', annot), SIMPLIFY = F
    )
    # modifying parameters
    if(reduce){
        annot.cl <- lapply(annot.cl, function(x){
            res <- x[x$Count == max(x$Count) | x$Benjamini == min(x$Benjamini), ]
        })
    }
    res <- Reduce(rbind, annot.cl)
    if(short){
        res <- res[, c('Cluster', 'Enrichment.Score', 'Category', 'Term', 'Count',
                       'Percent', 'PValue', 'Benjamini')]
    }
    # filter by adjusted P-values;
    res <- res[res$Benjamini < ben, ]
    return(res)
}

# x <- '~/tmp/t2t_B3AAF42FA92B1461018944337.txt'
# x2 <- parseAnnotCluster(x, reduce = T)


parseAnnotTable <- function(raw_file, short = T){
    # general parsing of david annotation table
    david <- read.table(raw_file, stringsAsFactors = F, sep='\t', quote='', header = T)
    colnames(david)[4] <- 'Percent'
    res <- david[, c(
            'Category', 'Term',  'Percent', 'Count', 'List.Total', 'Pop.Hits', 'Pop.Total',
            'Fold.Enrichment', 'PValue', 'Bonferroni', 'Benjamini', 'FDR')]
    if(short){
        res <- res[, c('Category', 'Term', 'Count', 'Percent', 'PValue', 'Benjamini')]
    }
    return(res)
}


parseAnnotTableGo <- function(raw_file, short = T){
    # if there are only terms with GO id, used this function to seperate id and description
    david <- read.table(
        raw_file,
        stringsAsFactors = F, sep='\t', quote='', header = T
    )

    david$Term.Id <- sapply(strsplit(david$Term, '~'), function(x) x[1])
    david$Term.Desp <- sapply(strsplit(david$Term, '~'), function(x) x[2])
    colnames(david)[4] <- 'Percent'

    res <- david[
        , c(
            'Category', 'Term.Id', 'Term.Desp',  'Percent', 'Count', 'List.Total', 'Pop.Hits',
            'Pop.Total', 'Fold.Enrichment', 'PValue', 'Bonferroni', 'Benjamini', 'FDR')]
    if(short){
        res <- res[, c('Category',  'Term.Id', 'Term.Desp', 'Count', 'Percent', 'PValue', 'Benjamini')]
    }
    return(res)
}

# A robust version based on `data.table` that parse all terms in each cluster and split Term into Term ID and Term
# description;
# specifically for David GO with only GO & KEGG terms;
ParseDavidCluster <- function(david.res.path){
    require(data.table)
    x <- strsplit(readLines(david.res.path), '\t')
    x <- x[sapply(x, length) == 13]
    h <- x[[1]]
    x <- unique(as.data.table(transpose(x)))
    x <- x[2:nrow(x), ]
    x[, names(x) := lapply(.SD, type.convert, as.is = TRUE)]
    setnames(x, names(x), h)
    x[, c('Term.ID', 'Term.desp') := tstrsplit(Term, split = '(?<=\\d)[:~]', perl = TRUE)]
    x[, Term.desp := paste0(toupper(substring(Term.desp, 1, 1)), substring(Term.desp, 2))]

    y <- setNames(
        c('GO Biological pathway', 'GO Molecular function', 'GO Cellular component', 'KEGG pathway'),
        c('GOTERM_BP_DIRECT', 'GOTERM_MF_DIRECT', 'GOTERM_CC_DIRECT', 'KEGG_PATHWAY'))
    x[, Category := y[Category]]

    res <- x[Benjamini < 0.05, .(
        Category, Term.ID, Term.desp, Count,
        Pvalue = ifelse(PValue < 0.01, sprintf('%.1e', PValue), sprintf('%.2f', PValue)),
        Benjamini = ifelse(Benjamini < 0.01, sprintf('%.1e', Benjamini), sprintf('%.2f', Benjamini))
    )]
    return(res)
}
