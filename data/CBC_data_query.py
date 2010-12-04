"""MySQL Querying of CBC to yield data for METE analyses"""

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

# Step 1. Group By SPECIES_CODE and by TOO WHERE DIURNAL LANDBIRD = 1 - 
    # to yield sp_too

cursor.execute("""
                CREATE TABLE queries.spcodes
                SELECT SPECIES.SPECIES_CODE FROM CBC.SPECIES 
                GROUP BY SPECIES.SPECIES_CODE;
                """) 

cursor.execute("""
                CREATE TABLE queries.sp_too_1
                SELECT SPECIES_CODE AS SPCODE, 
                TAXON_ORDER_OUT AS TOO FROM queries.spcodes 
                LEFT JOIN TAXONOMY_BIRDS.TAXONOMY 
                ON spcodes.SPECIES_CODE = TAXONOMY.CBCSPCODE
                WHERE TAXONOMY.DIURNALLANDBIRD = 1;
                """) 

cursor.execute("""
                CREATE TABLE queries.sp_too 
                SELECT SPCODE, TOO FROM queries.sp_too_1 
                GROUP BY SPCODE, TOO;
                """) 

# Step 2. LINK OBS and SUB_AUX by SUB_ID, adding COUNT_YR = 109 (2008-2009)
    # to yield OBSDATA_CTYR_STEP1 (SUB_ID, COUNT_YR, SPECIES_CODE, HOW_MANY)

cursor.execute("""
                CREATE TABLE queries.obs_1
                SELECT SUB_AUX.SUB_ID, SUB_AUX.COUNT_YR, OBS.SPECIES_CODE, OBS.HOW_MANY
                FROM CBC.SUB_AUX INNER JOIN CBC.OBS ON SUB_AUX.SUB_ID = OBS.SUB_ID
                WHERE SUB_AUX.COUNT_YR = 109;
                """)

# Step 3. LINK obs_1 to SUB by SUB_ID to add LOC_ID 

cursor.execute("""
                CREATE TABLE queries.obs_2
                SELECT SUB.LOC_ID, obs_1.* 
                FROM queries.obs_1 INNER JOIN CBC.SUB ON obs_1.SUB_ID = SUB.SUB_ID;
                """)

# Step 4. LINK obs_2 to LOC by LOC_ID to remove records WHERE SUBNATIONAL1_CODE 
    # IS NOT US-HI AND COUNTRY CODE IS 'CA', 'US', OR 'US-CA'

cursor.execute("""
                CREATE TABLE queries.obs_3
                SELECT obs_2.* 
                FROM queries.obs_2 INNER JOIN CBC.LOC ON obs_2.LOC_ID = LOC.LOC_ID
                WHERE LOC.SUBNATIONAL1_CODE != 'US-HI' AND 
                LOC.COUNTRY_CODE = 'CA' OR 'US' OR 'US-CA';
                """)

# Step 5. CREATE obs_4 by adding sp_too.TOO

cursor.execute("""
                CREATE TABLE queries.obs_4
                SELECT obs_3.*, sp_too.TOO
                FROM queries.obs_3 INNER JOIN queries.sp_too 
                ON obs_3.SPECIES_CODE = sp_too.SPCODE;
                """)

# Step 6. To create table with SiteID - Year - Sp - abund:
    # GROUPing BY LOC_ID and SPECIES_CODE and SUMmming over HOW_MANY FROM obs_4

cursor.execute("""
                SELECT obs_4.LOC_ID, obs_4.COUNT_YR AS YEAR, obs_4.TOO, 
                SUM(obs_4.HOW_MANY) AS AB
                FROM queries.obs_4 
                GROUP BY obs_4.LOC_ID, obs_4.COUNT_YR, obs_4.TOO
                INTO OUTFILE '/tmp/cbc_too_1093.csv'
                FIELDS TERMINATED BY ',' 
                LINES TERMINATED BY '\n';
                """)
               
connection.commit()

shutil.copy('/tmp/cbc_too_1093.csv', '/home/kate/data/cbc_too_109.csv')