"""MySQL Querying of BBS to yield data for METE analyses"""

from get_data import get_data

get_data([
# Step 1. TAXONOMY_BIRDS TAXONOMY.AOU_IN linked to BBS counts.Aou
    # Left Join = ALL records from BBS counts are included 
    # and only values from TAXONOMY that match
# Step 2. Group By AOU IN and by TOO WHERE DIURNAL LANDBIRD = 1 - 
    # to yield AOU_TOO
          """CREATE TABLE queries.aous (aou INT(11)) 
             SELECT counts.Aou FROM BBS.counts 
             GROUP BY counts.Aou;
             """,
          """CREATE TABLE queries.aou_too_1 (aou INT(11)) 
             SELECT Aou AS AOU, TAXON_ORDER_OUT AS TOO FROM queries.aous 
             LEFT JOIN TAXONOMY_BIRDS.TAXONOMY 
             ON aous.Aou = TAXONOMY.AOU_IN 
             WHERE TAXONOMY.DIURNALLANDBIRD = 1;
             """,
          """CREATE TABLE queries.aou_too (aou INT(11)) 
             SELECT AOU, TOO FROM queries.aou_too_1 
             GROUP BY AOU, TOO;
             """,
# 3. To create table with SiteID - Year - RunType = 1:
    # from weather table
          """CREATE TABLE queries.weather_subquery
             SELECT (weather.statenum*1000+ weather.Route) AS SiteID, weather.Year, weather.RunType
             FROM BBS.weather
             WHERE weather.RunType = 1 AND weather.RPID = 101;
             """,
# 4. To create table with SiteID - Year - Sp - abund:
    # Link together AOU_TOO and BBS Counts by AOU - 
    # Group By SiteID = state*1000 + route - TOO - Year 
    # Sum SpeciesTotal
          """CREATE TABLE queries.counts_too
             SELECT (counts.statenum * 1000) + counts.Route AS SiteID,
             counts.Year, aou_too.TOO, 
             SUM(counts.SpeciesTotal) AS AB 
             FROM BBS.counts INNER JOIN queries.aou_too ON counts.Aou = aou_too.AOU
             GROUP BY SiteID, counts.Year, aou_too.TOO, counts.RPID
             HAVING (((counts.Year = 2009) AND (counts.RPID = 101)));
             """,
          """SELECT counts_too.SiteID, counts_too.Year, counts_too.TOO, counts_too.AB
             FROM queries.counts_too INNER JOIN queries.weather_subquery
             ON counts_too.SiteID = weather_subquery.SiteID 
             AND counts_too.Year = weather_subquery.Year
             INTO OUTFILE '/tmp/bbs_too_2009.csv'
             FIELDS TERMINATED BY ',' 
             LINES TERMINATED BY '\n';
             """,
])
