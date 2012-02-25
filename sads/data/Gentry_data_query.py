"""MySQL Querying of Gentry to yield data for METE analyses"""

from get_data import get_data

# Group By SITE and SPECIES_ID, COUNT stem_id
get_data(["""CREATE TABLE queries.gentry
             SELECT Stems.site, Stems.species_id, Stems.Morpho,
             Stems.id_status, sum(`Stems`.`N(Ind)`) AS ab
             FROM gentry_glenda.Stems
             GROUP BY Stems.site, Stems.species_id
             HAVING (((ab > 0) AND (Stems.Morpho = 1 OR Stems.id_status = "species")));
             """,
          """SELECT gentry.site, gentry.species_id, gentry.ab FROM queries.gentry
             INTO OUTFILE '/tmp/gentry_spab.csv'
             FIELDS TERMINATED BY ',' 
             LINES TERMINATED BY '\n';
             """,
])
