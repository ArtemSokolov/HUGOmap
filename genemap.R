## Creates a map from ENSEMBL IDs to HUGO for protein-coding regions
##
## by Artem Sokolov

library( tidyverse )

## Location of the latest .gtf release
fnGTF <- "ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz"
fnLocal <- "GRCh38.97.gtf.gz"

main <- function()
{
    ## Download the raw .gtf file
    download.file( fnGTF, fnLocal )

    ## Parse the raw .gtf file
    cn <- c( "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute" )
    ct <- cols( .default = col_character(), start = col_integer(), end = col_integer() )
    X <- read_tsv( fnLocal, comment="#!", col_types=ct, col_names = cn )

    ## Isolate the gene entries and compute gene lengths
    Y <- X %>% filter( feature == "gene" ) %>% mutate( gene_length = end-start+1 ) %>%
        select( attribute, gene_length )

    ## Parse key-value pairs
    pkv <- function(v)
    {
        str_sub(v, 1,-2) %>% str_split("; ") %>% map(str_split, " ") %>%
            map( ~set_names(map(.x, pluck, 2), map(.x, pluck, 1)) ) %>%
            modify_depth( 2, str_sub, 2, -2 )
    }
    G <- Y %>% mutate_at("attribute", pkv) %>%
        mutate_at( "attribute", map2, .$gene_length, ~c(.x, list(gene_length=.y)) ) %>%
        pull( attribute ) %>% bind_rows()
    
    ## Save a copy to a file
    G %>% write_csv( "GRCh38.97.genemap.csv" )
}
