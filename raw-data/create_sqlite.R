####load packages####
library("tidyverse")
library("DBI")
library("RSQLite")
library("glue")

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
env_mdb <- "raw-data/define.mdb"
spp_mdb <- "raw-data/processCounts.mdb"

# ## export env data
fchem <- system(glue("mdb-export -H -I sqlite {env_mdb} fchem"), intern = TRUE)

## export env schema
env_setup <- system(glue("mdb-schema {env_mdb} -T fchem"), intern = TRUE) %>%
  paste(collapse = "\n")

#setup fchem table
dbExecute(statement = env_setup, conn = con)

#import fchem data
walk(fchem, ~dbExecute(conn = con, statement = .x))

#check fchem
tbl(con, "fchem") %>%
  collect()




#check what tables have been added
dbListTables(con)

### cleanup ####
dbDisconnect(conn = con)




#extract raw counts  & merges from access files
#read sqlite database and process counts

