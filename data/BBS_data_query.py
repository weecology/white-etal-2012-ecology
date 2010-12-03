"""MySQL Querying of BBS to yield data for METE analyses"""

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

# Step 1. TAXONOMY_BIRDS TAXONOMY.AOU_IN linked to BBS counts.Aou
    # Left Join = ALL records from BBS counts are included 
    # and only values from TAXONOMY that match
# Step 2. Group By AOU IN and by TOO WHERE DIURNAL LANDBIRD = 1 - 
    # to yield AOU_TOO

cursor.execute("""
                CREATE TABLE queries.aous (aou INT(11)) 
                SELECT counts.Aou FROM BBS.counts 
                GROUP BY counts.Aou;
                """) 

cursor.execute("""
                CREATE TABLE queries.aou_too_1 (aou INT(11)) 
                SELECT Aou AS AOU, TAXON_ORDER_OUT AS TOO FROM queries.aous 
                LEFT JOIN TAXONOMY_BIRDS.TAXONOMY 
                ON aous.Aou = TAXONOMY.AOU_IN 
                WHERE TAXONOMY.DIURNALLANDBIRD = 1;
                """) 

cursor.execute("""
                CREATE TABLE queries.aou_too (aou INT(11)) 
                SELECT AOU, TOO FROM queries.aou_too_1 
                GROUP BY AOU, TOO;
                """) 

# 3. To create table with Year - SiteID - Sp - abund:
    # Link together AOU_TOO and BBS Counts by AOU - 
    # Group By SiteID = state*1000 + route - TOO - Year 
    # Sum SpeciesTotal
    
cursor.execute("""
                CREATE TABLE queries.bbs_too
                SELECT (counts.statenum * 1000) + counts.Route AS SiteID,
                counts.Year, aou_too.TOO, 
                SUM(counts.SpeciesTotal) AS AB 
                FROM BBS.counts INNER JOIN queries.aou_too ON counts.Aou = aou_too.AOU
                GROUP BY (counts.statenum * 1000) + counts.Route, counts.Year, aou_too.TOO
                HAVING (((counts.Year) = 2009))
                INTO OUTFILE '/tmp/bbs_too_2009.csv';
                """)
                
connection.commit()

shutil.copy('/tmp/bbs_too_2009.csv', '/home/kate/data/bbs_too_2009.csv')