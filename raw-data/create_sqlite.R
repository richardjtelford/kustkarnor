####load packages####
library("tidyverse")
library("DBI")
library("RSQLite")
library("glue")

#### import function ####
import_table <- function(mdb, table, overwrite = TRUE){
  fix_names <- . %>%
    str_replace_all("/", "_")

  fix_data <- . %>%
    str_remove("\r") %>%
    str_remove("<br>")

  ### export data
  import_sql <- system(glue("mdb-export -H -I sqlite {mdb} {table}"), intern = TRUE) %>%
    fix_names() %>%
    fix_data() %>%
    #collapse and resplit on INSERT
    str_c(collapse = " ") %>%
    str_split("(?=INSERT)") %>%
    pluck(1)

  import_sql <- import_sql[import_sql != ""]#remove empty elements
  ## export schema
  schema <- system(glue("mdb-schema {mdb} -T {table}"), intern = TRUE) %>%
    paste(collapse = "\n") %>%
    fix_names()

  #drop table if it already exists and overwrite is TRUE
  if(isTRUE(overwrite) & table %in% dbListTables(con)){
    dbExecute(conn = con, glue("DROP TABLE {table}"))
  }

  #setup table
  dbExecute(statement = schema, conn = con)

  #import fchem data
  walk(import_sql, ~dbExecute(conn = con, statement = .x))
}



#### make sqlite3 database ####
#delete existing version
if(file.exists("data/define.sqlite")){
  file.remove("data/define.sqlite")
}

#make new empty database
con <- dbConnect(SQLite(), dbname = "data/define.sqlite")

#### import enviroment ####
import_table(mdb = "raw-data/define.mdb", table = "fchem")

#### import taxonomy ####
import_table(mdb = "raw-data/definetaxonomy.mdb", table = "cleanTaxa")
import_table(mdb = "raw-data/definetaxonomy.mdb", table = "Synonyms")
import_table(mdb = "raw-data/definetaxonomy.mdb", table = "Merges")
import_table(mdb = "raw-data/processCounts.mdb", table = "MergeCodes")

#counts
import_table(mdb = "raw-data/defineCounts.mdb", table = "AllCounts")
import_table(mdb = "raw-data/moltenSurfaceCounts2000.mdb", table = "MCodedCounts")


#excluded taxa
import_table(mdb = "raw-data/processCounts.mdb", table = "ExcludedTaxa")


#union count data
bind_rows(
  tbl(con, "MCodedCounts") %>%
    rename(siteID = siteId, taxonCode = MoltenCode) %>%
    collect(),
  tbl(con, "AllCounts") %>%
    rename(count = cnt) %>%
    mutate(count = as.numeric(count)) %>%
    collect()
) %>%
  dbWriteTable(conn = con, name = "counts")

#drop original counts
db_drop_table(con = con, "AllCounts")
db_drop_table(con = con, "MCodedCounts")



#check what tables have been added
dbListTables(con)

### cleanup ####
dbDisconnect(conn = con)
