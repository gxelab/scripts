library(xml2)
library(dplyr)
library(tibble)


parse_sra_xml <- function(xml_path){
    xml <- read_xml(xml_path)
    exps <- xml %>% xml_children()
    
    get_sample_attr <- function(y){
        tags <- y %>% xml_children() %>% xml_find_all('.//TAG') %>% xml_text()
        values <- y %>% xml_children() %>% xml_find_all('.//VALUE') %>% xml_text()
        names(values) <- tags
        return(values)
    }
    
    bind_rows(lapply(exps, function(node){
        exp_acc <- node %>% xml_child('EXPERIMENT') %>% xml_attr('accession')
        exp_title <- node %>% xml_child('EXPERIMENT') %>% xml_child('TITLE') %>% xml_text()
        
        lib <- node %>% xml_find_first('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR')
        lib_strategy <- lib %>% xml_child('LIBRARY_STRATEGY') %>% xml_text()
        lib_source <- lib %>% xml_child('LIBRARY_SOURCE') %>% xml_text()
        lib_selection <- lib %>% xml_child('LIBRARY_SELECTION') %>% xml_text()
        lib_layout <- lib %>% xml_child('LIBRARY_LAYOUT') %>% xml_child() %>% xml_name()
        
        platform <- node %>% xml_find_first('./EXPERIMENT/PLATFORM') %>%
            xml_child() %>% xml_name()
        model <- node %>% xml_find_first('./EXPERIMENT/PLATFORM') %>%
            xml_child() %>% xml_child('INSTRUMENT_MODEL') %>% xml_text()
        samp <- node %>% xml_child('SAMPLE')
        sample_acc <- samp %>% xml_attr('accession')
        sample_title <- samp %>% xml_child('TITLE') %>% xml_text()
        taxon <- samp %>% xml_child('SAMPLE_NAME') %>% xml_child('TAXON_ID') %>% xml_text()
        sci_name <- samp %>% xml_child('SAMPLE_NAME') %>% xml_child('SCIENTIFIC_NAME') %>% xml_text()
        species <- xml %>%
            xml_find_all('.//EXPERIMENT_PACKAGE/SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME') %>% xml_text()
        study <- node %>% xml_child('STUDY') %>% xml_attr('accession')
        geo <- node %>%
            xml_find_first('./STUDY/IDENTIFIERS/EXTERNAL_ID[@namespace="GEO"]') %>%
            xml_text()
        dt_meta <- tibble(exp_acc, exp_title, study, geo,
                          lib_strategy, lib_source, lib_selection, lib_layout,
                          platform, model, sample_acc, sample_title, taxon, sci_name)
        
        runs <- node %>% xml_find_all('./RUN_SET/RUN')
        run_acc <- runs %>% xml_attr('accession')
        run_url <- runs %>%
            xml_find_first('SRAFiles/SRAFile[@semantic_name="run"]/Alternatives[@org="NCBI"]') %>%
            xml_attr('url')
        dt_run <- tibble(run_acc, run_url)
        
        sample_attrs <- samp %>%
            xml_child('SAMPLE_ATTRIBUTES') %>%
            get_sample_attr() %>%
            as_tibble_row() %>%
            rename_with(~ gsub(' ', '_', .x))
        bind_cols(dt_run, dt_meta, sample_attrs)
    }))
}

parse_sra_xml('data/benchmark/SRP049168_Gao_2014_Nat_methods.xml')
parse_sra_xml('data/benchmark/2013_MBDC_RNA.xml')
