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
fchem <- system(glue("mdb-export -H -I sqlite {env_mdb} fchem"), intern = TRUE) %>%
  paste(collapse = "\n")

## export env schema
env_setup <- system(glue("mdb-schema {env_mdb} -T fchem"), intern = TRUE) %>%
  paste(collapse = "\n")

#setup fchem table
dbExecute(statement = env_setup, conn = con)

#import fchem data
dbExecute(conn = con, statement = fchem)


#check what tables have been added
dbListTables(con)


### cleanup ####
dbDisconnect(conn = con)




#extract raw counts  & merges from access files
#create sqlite database
#read sqlite database and process counts

