####load packages####
library("tidyverse")
library("DBI")
library("RSQLite")
library("glue")

#### import function ####
import_table <- function(mdb, table){
  ### export data
  import_sql <- system(glue("mdb-export -H -I sqlite {mdb} {table}"), intern = TRUE)

  ## export schema
  schema <- system(glue("mdb-schema {env_mdb} -T {table}"), intern = TRUE) %>%
    paste(collapse = "\n")

  #setup fchem table
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

# path to exported data
system('mkdir -p data/sql')

#raw data files
spp_mdb <- "raw-data/processCounts.mdb"

#### enviroment ####
import_table(mdb = "raw-data/define.mdb", table = "fchem")


#check fchem
tbl(con, "fchem") %>%
  collect()

#excluded taxa
#defineCounts.mdb -defineCounts
#molten/moltenSurfaceCounts2000.mdb
#definetaxonomy.mdb synonyms & merges

#check what tables have been added
dbListTables(con)

### cleanup ####
dbDisconnect(conn = con)




#extract raw counts  & merges from access files
#read sqlite database and process counts

