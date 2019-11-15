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
    rename(taxonCode = MoltenCode) %>%
    collect(),
  tbl(con, "AllCounts") %>%
    rename(siteId = siteID, count = cnt) %>%
    mutate(count = as.numeric(count)) %>%
    collect()
) %>%
  dbWriteTable(conn = con, name = "counts")

#drop original counts
db_drop_table(con = con, "AllCounts")
db_drop_table(con = con, "MCodedCounts")

#fix cases in countryId
tbl(con, "Merges") %>%
  collect() %>%
  mutate(dataSet = str_to_title(dataSet)) %>%
  rename(countryId = dataSet) %>%
  dbWriteTable(conn = con, name = "Merges", overwrite = TRUE)

tbl(con, "fchem") %>%
  collect() %>%
  mutate(countryId = str_to_title(countryId)) %>%
  mutate(siteId = if_else(siteId == "Dk1", "DK1", siteId)) %>%
  dbWriteTable(conn = con, name = "fchem", overwrite = TRUE)


#fix case problems in spp names
case_fixes <- tribble(~correctCode,  ~synonymCode,
                      "RhiSpi", "RhiSPi",
                      "PlanHau", "Planhau")
dbWriteTable(conn = con, name = "Synonyms", value = case_fixes, append = TRUE)


#check what tables have been added
dbListTables(con)

### cleanup ####
dbDisconnect(conn = con)
