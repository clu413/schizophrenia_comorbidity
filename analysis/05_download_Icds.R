library(SqlServerJtds)
library(SqlTools)
library(pryr)

my.username="XXXXX"
my.password="YYYYY"

cn = connect.sql.server(
  database="ZZZZZ",
  domain="XXXXX",
  server.address="YYYYYYY",
  user=my.username,
  password=my.password)

output_dir = 'XXX/schizo/data/retrieved-data/case_ctrl_icds'

nonSchizo = chunked.query(cn, table.name="XXX.dbo.all_ctrl_icds", order.by.clause="ORDER BY ctrl_id", max.n.rows = 2e+6)
write.csv(nonSchizo, paste(output_dir, "/all_ctrl_icds.txt", sep=''), row.names=FALSE)
