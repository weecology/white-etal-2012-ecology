"""MySQL Querying of Gentry to yield data for METE analyses"""

import MySQLdb as dbapi
import getpass
import shutil

# enter MySQL password
p=getpass.getpass()

connection = dbapi.connect(host='129.123.92.244', port=1977, user='kate', 
                       passwd=p)
cursor = connection.cursor()

cursor.execute("""DROP DATABASE IF EXISTS queries;""")
cursor.execute("""CREATE DATABASE queries;""")

# Group By SITE and SPECIES_ID, COUNT stem_id
cursor.execute("""
                CREATE TABLE queries.gentry
                SELECT Stems.site, Stems.species_id, Stems.Morpho,
                Stems.id_status, sum(`Stems`.`N(Ind)`) AS ab
                FROM gentry_glenda.Stems
                GROUP BY Stems.site, Stems.species_id
                HAVING (((ab > 0) AND (Stems.Morpho = 1 OR Stems.id_status = "species")));
                """)
                
cursor.execute("""
                SELECT gentry.site, gentry.species_id, gentry.ab FROM queries.gentry
                INTO OUTFILE '/tmp/gentry_spab.csv'
                FIELDS TERMINATED BY ',' 
                LINES TERMINATED BY '\n';
                """)
               
connection.commit()

shutil.copy('/tmp/gentry_spab.csv', '/home/kate/data/gentry_spab.csv')