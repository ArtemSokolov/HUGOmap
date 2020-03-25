## Creates a map from ENSEMBL IDs to HUGO
##
## by Artem Sokolov

## Parse version number from the command line
argv <- commandArgs(TRUE)
if( length(argv) != 1 )
    stop( "Usage: Rscript genemap.R [build version]" )
v <- argv[1]

suppressMessages(library( tidyverse ))

## Location of the latest .gtf release
fnGTF <- str_c("ftp://ftp.ensembl.org/pub/release-", v,
               "/gtf/homo_sapiens/Homo_sapiens.GRCh38.", v,
               ".gtf.gz")

## Parse the raw .gtf file
cn <- c( "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute" )
ct <- cols( .default = col_character(), start = col_integer(), end = col_integer() )
X <- read_tsv( fnGTF, comment="#!", col_types=ct, col_names = cn )

## Isolate the gene entries and compute gene lengths
Y <- X %>% filter( feature == "gene" ) %>% mutate( gene_length = end-start+1 )

## Parse key-value pairs
G <- Y$attribute %>% str_sub(1,-2) %>% str_split("; ") %>% map(str_split, " ") %>%
    map( ~set_names(map(.x, pluck, 2), map(.x, pluck, 1)) ) %>%
    modify_depth( 2, str_sub, 2, -2 ) %>%
    map2( Y$gene_length, ~c(.x, list(gene_length=.y)) ) %>% bind_rows()
    
## Save a copy to a file
G %>% write_csv( str_c("GRCh38.",v,".genemap.csv") )
