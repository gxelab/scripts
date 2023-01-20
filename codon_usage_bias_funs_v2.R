# TODO: how to deal with codon that does not appear in cf when calculate w for CAI

# last modified: 2021-03-31
# author: mt1022
# ------------------------------------------------
library(Biostrings)

# codon amino acid table: GENETIC_CODE

Base2Codon <- function(x){
    "convert CDS sequence to three letter codons"
    res <- substring(toString(x), seq(1, length(x), 3), seq(3, length(x), 3))
    return(res)
}

CodonFreq <- function(seqs){
    "tabulate codon frequency"
    cf <- oligonucleotideFrequency(seqs, width = 3, step = 3)
    # trinucleotideFrequency(seqs, step = 3)

    # return
    if(class(seqs) == 'DNAStringSet'){
        return(colSums(cf))
    }else{
        return(cf)
    }
}

Rscu <- function(seqs){
    "calculate RSCU"
    cf <- CodonFreq(seqs)

    # codon-aminoAcid-freq table
    codon.index <- as.data.frame(GENETIC_CODE)
    codon.index$codon <- rownames(codon.index)
    codon.index$codon.freq <- 0
    codon.index[names(cf), 'codon.freq'] <- cf

    # mean/max syn codon freq
    x <- aggregate(codon.index$codon.freq, list(codon.index$GENETIC_CODE), mean)
    cf.mean <- x[, 2]; names(cf.mean) <- x[, 1]
    codon.index <- merge(codon.index, cf.mean, by.x = 1, by.y = 0)
    colnames(codon.index)[4] <- 'mean.syn.codon.freq'

    y <- aggregate(codon.index$codon.freq, list(codon.index$GENETIC_CODE), max)
    cf.max <- y[, 2]; names(cf.max) <- y[, 1]
    codon.index <- merge(codon.index, cf.max, by.x = 1, by.y = 0)
    colnames(codon.index)[5] <- 'max.syn.codon.freq'

    # calc rscu (with pseudo-count for codon family with no occurence)
    codon.index$rscu <- codon.index$codon.freq / (codon.index$mean.syn.codon.freq + 1e-9)

    # calc CAI w (with pseudo-count for codon family with no occurence)
    codon.index$w <- codon.index$codon.freq / (codon.index$max.syn.codon.freq + 1e-9)

    return(codon.index)
}

Cai <- function(seq, codon.index){
    "calculate CAI"
    # seq is an BSstring object
    # codon.index is the output of Rscu
    seq.codon.freq <- CodonFreq(seq)
    logsum <- sum(seq.codon.freq[codon.index$codon] * log(codon.index$w))
    cai <- exp( logsum / sum(seq.codon.freq))
    return(cai)
}

CodonwCaiFile <- function(seqs_in, cai.file = ""){
    "generate cai ref file for codonW"
    seqs <- readDNAStringSet(seqs_in)
    codon.index <- Rscu(seqs)

    codon <- c()
    b <- c('T', 'C', 'A', 'G')
    for(i in b){for(j in b){for(k in b){codon <- c(codon, paste0(i, k, j))}}}

    codonw <- codon.index$w
    names(codonw) <- codon.index$codon
    codonw[codon.index$GENETIC_CODE == '*'] <- 1

    cat(as.character(codonw[codon]), sep = '\n', file = cai.file)
    return(NULL)
}

TrnaW <- function(x){
    # get tRNA W for calculating TAI
    # x: vector of copy number of tRNAs grouped by anti-codons, anti-codons are given in names attribute;

    # build data.frame
    cdn <- data.frame(
        GENETIC_CODE = GENETIC_CODE,
        codon = names(GENETIC_CODE),
        anti_codon = sapply(DNAStringSet(names(GENETIC_CODE)), function(x) as.character(reverseComplement(x))),
        stringsAsFactors = F
    )
    # process input
    rownames(cdn) <- cdn$anti_codon
    x <- x[names(x) %in% cdn$anti_codon]
    # add tRNA copy number
    cdn$cn <- 0
    cdn[names(x), 'cn'] <- x
    # omit Selenocysteine
    cdn['TCA', 'cn'] <- 0
    #' tRNA codon-anticodon constraints in Eukarya are from Table 4 of:
    #' Renana Sabi, Tamir Tuller, Modelling the Efficiency of Codon–tRNA Interactions Based on
    #' Codon Usage Bias, DNA Research, Volume 21, Issue 5, October 2014, Pages 511–526,
    #' https://doi.org/10.1093/dnares/dsu017
    s.iu <- 0; s.gc <- 0; s.cg <- 0; s.ua <- 0
    s.gu <- 0.7861; s.ic <- 0.4659; s.ia <- 0.9075; s.ug <- 0.6295
    # calculate W
    cdn$big_w <- 0
    for(i in seq(1, 64, 4)){
        cdn$big_w[i] <- (1 - s.iu) * cdn$cn[i] + (1 - s.gu) * cdn$cn[i + 1]
        cdn$big_w[i + 1] <- (1 - s.gc) * cdn$cn[i + 1] + (1 - s.ic) * cdn$cn[i]
        cdn$big_w[i + 2] <- (1 - s.ua) * cdn$cn[i + 2] + (1 - s.ia) * cdn$cn[i]
        cdn$big_w[i + 3] <- (1 - s.cg) * cdn$cn[i + 3] + (1 - s.ug) * cdn$cn[i + 2]
    }
    cdn$small_w <- cdn$big_w / max(cdn$big_w)
    cdn$small_w[cdn$small_w == 0] <- mean(cdn$small_w[cdn$small_w > 0])
    return(cdn)
}

Tai <- function(seq, codon.index){
    "calculate TAI"
    # seq is an BSstring object
    # codon.index is the output of TaiW
    seq.codon.freq <- CodonFreq(seq)
    logsum <- sum(seq.codon.freq[codon.index$codon] * log(codon.index$small_w))
    tai <- exp( logsum / sum(seq.codon.freq))
    return(tai)
}

# clean CDS sequence
# @param: dna --> DNAstringSet instance
# @return: DNAstringset inscance
# @bug: this function is extemely slow now.
CleanCDS <- function(seqs){
    # remove non 3k seq
    seqs <- seqs[sapply(seqs, length) %% 3 == 0]
    print(length(seqs))
    # check start and stop codon and internal codons
    k <- sapply(seqs, function(seq){
        seq.codon <- Base2Codon(seq)
        if(seq.codon[1] != 'ATG'){
            return(F)
        }else if(! seq.codon[length(seq.codon)] %in% c('TAA', 'TAG', 'TGA')){
            return(F)
        }else if(any(c('TAA', 'TAG', 'TGA') %in% seq.codon[2:(length(seq.codon) - 1)])){
            return(F)
        }else{
            return(T)
        }
    })
    print(sum(k))
    seqs <- seqs[k]
    return(seqs)
}

# CMH test for AT/CG bias in synonymous codon usage between two set of genes
RscuCMHtest <- function(d1, d2){
    # d1, d2: return of Rscu for gene set1 and set2;
    pp.table <- function(d){
        d <- d[d$GENETIC_CODE != '*', ]
        d$GENETIC_CODE <- droplevels(d$GENETIC_CODE)
        d.fam <- split(d, f = paste(d$GENETIC_CODE, substr(d$codon, 1, 2), sep = '_'))
        d.fam <- d.fam[sapply(d.fam, nrow) > 1]
        res <- sapply(d.fam, function(x){
            res <- c(
                AT = sum(x$codon.freq[substr(x$codon, 3, 3) %in% c('A', 'T')]),
                CG = sum(x$codon.freq[substr(x$codon, 3, 3) %in% c('C', 'G')])
            )
        })
    }
    x1 <- pp.table(d1)
    x2 <- pp.table(d2)
    chi.tab <- rbind(x1, x2)[c(1, 3, 2, 4), ]
    chi.tab <- array(
        chi.tab, dim = c(2, 2, ncol(x1)),
        dimnames = list(src = c('d1', 'd2'), syncodon = c('AT', 'CG'), fam = colnames(x1))
    )
    res <- mantelhaen.test(chi.tab)
}

# Preferred codons
#' Reference:
#'   Akashi H. (1995). Inferring weak selection from patterns of polymorphism and
#'   divergence at "silent" sites in Drosophila DNA. Genetics, 139(2), 1067–1076.
#' According to legend of Table 1 in Akashi 1995:
#'   Serine codons are divided into a fourfold and twofold degenerate family.

PreferredCodon <- function(seq_introns, seq_cds_rep){
    #' seq_introns: sequence of introns
    #' seq_cds_rep: CDS sequences of all protein-coding genes
    #'              (one representative transcript per gene)
    #' return     : data.frame, codons with rho > 0 (significantly) are preferred codons              
    # codon family
    synfam <- cbind.data.frame( 
        codon = names(GENETIC_CODE),aa = GENETIC_CODE, stringsAsFactors = F)
    rownames(synfam) <- NULL
    synfam$fam <- synfam$aa
    synfam$fam[synfam$codon %in% c('AGT', 'AGC')] <- 'S2'
    synfam$fam[synfam$codon %in% c('AGA', 'AGG')] <- 'R2'
    synfam$fam[synfam$codon %in% c('TTA', 'TTG')] <- 'L2'
    synfam <- synfam[!synfam$aa %in% c('*', 'M', 'W'), ]
    
    # proportion of A/T in introns (PAT)
    intron.base <- colSums(oligonucleotideFrequency(seq_introns, width = 1, step = 1))
    pat <- sum(intron.base[c('A', 'T')]) / sum(intron.base)
    
    # filter cds
    FilterSeq <- function(cds.seq){
        cds.seq <- cds.seq[width(cds.seq) %% 3 == 0]
        cds.seq <- cds.seq[subseq(cds.seq, 1, 3) == 'ATG']
        cds.seq <- cds.seq[
            subseq(cds.seq, width(cds.seq) - 2, width(cds.seq)) %in% c('TAA', 'TAG', 'TGA')]
    }
    seq_cds_rep <- FilterSeq(seq_cds_rep)
    
    # chi-square analysis
    rep.gene.codon <- trinucleotideFrequency(seq_cds_rep, step = 3)
    rep.gene.codon <- as.data.frame(rep.gene.codon)
    chisq <- lapply(split(synfam$codon, synfam$fam), function(x){
        x1 <- rowSums(rep.gene.codon[x[substr(x, 3, 3) %in% c('A', 'T')]])
        x2 <- rowSums(rep.gene.codon[x]) - x1
        atcg <- cbind.data.frame(at = x1, cg = x2)
        atcg$chisq <- apply(atcg, 1, function(y){
            if(any(y > 0)) chisq.test(y, p = c(pat, 1-pat))$statistic else 0
        })
        atcg
    })
    scaled.chisq <- sapply(names(chisq), function(x){
        other.chisq <- chisq[names(chisq) != x]
        tot.codons <- rowSums(sapply(other.chisq, function(d) rowSums(d[1:2])))
        tot.chisq <- rowSums(sapply(other.chisq, function(d) d$chi))
        tot.chisq / tot.codons
    })
    scaled.chisq <- as.data.frame(scaled.chisq)
    
    # determine preferred codons
    preferred.codon <- lapply(synfam$codon, function(x){
        fam <- synfam$fam[synfam$codon == x]
        codon.freq <- rep.gene.codon[[x]] /
            rowSums(rep.gene.codon[synfam$codon[synfam$fam == fam]])
        gene.schisq <- scaled.chisq[[fam]]
        t1 <- cor.test(codon.freq, gene.schisq, method = 'spearman')
        data.frame(codon = x, rho = t1$estimate, pval = t1$p.value, stringsAsFactors = F)
    })
    preferred.codon <- do.call(rbind, preferred.codon)
    preferred.codon$pval.adj <- p.adjust(preferred.codon$pval, method = 'BH')
    preferred.codon <- merge(synfam, preferred.codon, by = 'codon')
    return(preferred.codon)
}

make_codon_table <- function(pfcodons = NULL){
    #' make a codon table
    codon_table <- data.table(aa_code = GENETIC_CODE, codon = names(GENETIC_CODE))
    codon_table <- codon_table[!aa_code %in% c('*')]
    codon_table[, amino_acid := AMINO_ACID_CODE[aa_code]]
    codon_table[, `:=`(fam = aa_code, family = amino_acid)]
    codon_table[codon %in% c('TTA', 'TTG', 'AGA', 'AGG', 'AGC', 'AGT'),
                `:=`(fam = paste0(fam, '2'), family = paste0(family, '2'))]
    if(!is.null(pfcodons)){
        codon_table[pfcodons, pf_type := i.pf_type, on = .(codon)]
        codon_table[is.na(pf_type), pf_type := 'U']
    }
    return(codon_table)
}

get_genomic_bias <- function(seq_cds, pat=0.6){
    #' get genomic bias of each CDS sequence using Akashi's method
    #' @param seq_cds CDS sequence of each gene
    #' @param pat percentage of A/T content in introns
    codon_family <- make_codon_table()
    codon_family <- codon_family[!amino_acid %in% c('Met', 'Trp')]
    codon_family <- split(codon_family$codon, codon_family$family)
    tx_cf <- trinucleotideFrequency(seq_cds, step = 3)
    tx_cf <- tx_cf[, colnames(tx_cf) %in% unlist(codon_family)]
    tx_chisq <- sapply(codon_family, function(fam){
        m <- tx_cf[, fam]
        x <- substr(colnames(m), 3, 3)
        mat <- rowSums(m[, x %in% c('A', 'T'), drop = FALSE])
        mgc <- rowSums(m[, x %in% c('G', 'C'), drop = FALSE])
        chisq <- mapply(function(x1, x2){
            if(x1 + x2 > 0) chisq.test(c(x1, x2), p = c(pat, 1-pat))$statistic else 0
        }, mat, mgc)
    })
    data.table(tx_name = names(seq_cds), chi = rowSums(tx_chisq)/rowSums(tx_cf))
}
