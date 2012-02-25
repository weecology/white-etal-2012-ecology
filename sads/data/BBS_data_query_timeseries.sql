CREATE TABLE queries.aous (aou INT(11)) 
SELECT counts.Aou FROM BBS.counts 
GROUP BY counts.Aou;

CREATE TABLE queries.aou_too_1 (aou INT(11)) 
SELECT Aou AS AOU, TAXON_ORDER_OUT AS TOO FROM queries.aous 
LEFT JOIN TAXONOMY_BIRDS.TAXONOMY 
ON aous.Aou = TAXONOMY.AOU_IN 
WHERE TAXONOMY.DIURNALLANDBIRD = 1;

CREATE TABLE queries.aou_too (aou INT(11)) 
SELECT AOU, TOO FROM queries.aou_too_1 
GROUP BY AOU, TOO;

CREATE TABLE queries.weather_subquery
SELECT (weather.statenum*1000+ weather.Route) AS SiteID, weather.Year, weather.RunType
FROM BBS.weather
WHERE weather.RunType = 1 AND weather.RPID = 101;

CREATE TABLE queries.bbs_site_yr
SELECT (counts.statenum * 1000) + counts.Route AS SiteID,
counts.Year FROM BBS.counts
GROUP BY SiteID, Year;

CREATE TABLE queries.bbs_num_yr
SELECT bbs_site_yr.SiteID, COUNT(bbs_site_yr.Year) AS num_yrs FROM queries.bbs_site_yr
INNER JOIN queries.weather_subquery
ON bbs_site_yr.SiteID = weather_subquery.SiteID AND bbs_site_yr.Year = weather_subquery.Year
GROUP BY SiteID HAVING ((num_yrs > 29));

CREATE TABLE queries.counts_too_ts
SELECT (counts.statenum * 1000) + counts.Route AS SiteID,
counts.Year, aou_too.TOO, counts.RPID,
SUM(counts.SpeciesTotal) AS AB 
FROM BBS.counts INNER JOIN queries.aou_too ON counts.Aou = aou_too.AOU
GROUP BY SiteID, counts.Year, aou_too.TOO, counts.RPID
HAVING (((counts.RPID = 101)));

CREATE TABLE queries.counts_ts_final
SELECT counts_too_ts.SiteID, counts_too_ts.Year, counts_too_ts.TOO, counts_too_ts.AB
FROM queries.counts_too_ts INNER JOIN queries.bbs_num_yr
ON counts_too_ts.SiteID = bbs_num_yr.SiteID;

SELECT * FROM queries.counts_ts_final
INTO OUTFILE '/tmp/bbs_too_ts.csv'
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n';