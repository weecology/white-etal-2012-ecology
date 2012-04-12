"""MySQL Querying of CBC to yield data for METE analyses"""

from get_data import get_data

get_data([
# Step 1. Group By SPECIES_CODE and by TOO WHERE DIURNAL LANDBIRD = 1 - 
    # to yield sp_too
          """CREATE TABLE queries.spcodes
             SELECT SPECIES.SPECIES_CODE FROM CBC.SPECIES 
             GROUP BY SPECIES.SPECIES_CODE;
             """, 
          """CREATE TABLE queries.sp_too_1
             SELECT SPECIES_CODE AS SPCODE, 
             TAXON_ORDER_OUT AS TOO FROM queries.spcodes 
             LEFT JOIN TAXONOMY_BIRDS.TAXONOMY 
             ON spcodes.SPECIES_CODE = TAXONOMY.CBCSPCODE
             WHERE TAXONOMY.DIURNALLANDBIRD = 1;
             """, 
          """CREATE TABLE queries.sp_too 
             SELECT SPCODE, TOO FROM queries.sp_too_1 
             GROUP BY SPCODE, TOO;
             """, 
# Step 2. LINK OBS and SUB_AUX by SUB_ID, adding COUNT_YR = 109 (2008-2009)
    # to yield OBSDATA_CTYR_STEP1 (SUB_ID, COUNT_YR, SPECIES_CODE, HOW_MANY)
          """CREATE TABLE queries.obs_1
             SELECT SUB_AUX.SUB_ID, SUB_AUX.COUNT_YR, OBS.SPECIES_CODE, OBS.HOW_MANY
             FROM CBC.SUB_AUX INNER JOIN CBC.OBS ON SUB_AUX.SUB_ID = OBS.SUB_ID
             WHERE SUB_AUX.COUNT_YR = 109;
             """,
# Step 3. LINK obs_1 to SUB by SUB_ID to add LOC_ID 
          """CREATE TABLE queries.obs_2
             SELECT SUB.LOC_ID, obs_1.* 
             FROM queries.obs_1 INNER JOIN CBC.SUB ON obs_1.SUB_ID = SUB.SUB_ID;
             """,
# Step 4. LINK obs_2 to LOC by LOC_ID to remove records WHERE SUBNATIONAL1_CODE 
    # IS NOT US-HI AND COUNTRY CODE IS 'CA', 'US', OR 'US-CA'
          """CREATE TABLE queries.obs_3
             SELECT obs_2.* 
             FROM queries.obs_2 INNER JOIN CBC.LOC ON obs_2.LOC_ID = LOC.LOC_ID
             WHERE LOC.COUNTRY_CODE = 'CA' OR LOC.COUNTRY_CODE = 'US' 
             OR LOC.COUNTRY_CODE = 'US-CA' AND LOC.SUBNATIONAL1_CODE != 'US-HI';
             """,
# Step 5. CREATE obs_4 by adding sp_too.TOO
          """CREATE TABLE queries.obs_4
             SELECT obs_3.*, sp_too.TOO
             FROM queries.obs_3 INNER JOIN queries.sp_too 
             ON obs_3.SPECIES_CODE = sp_too.SPCODE;
             """,
# Step 6. To create table with SiteID - Year - Sp - abund:
    # GROUPing BY LOC_ID and SPECIES_CODE and SUMmming over HOW_MANY FROM obs_4
          """CREATE TABLE queries.obs_5                
             SELECT obs_4.LOC_ID, obs_4.COUNT_YR AS YEAR, obs_4.TOO, 
             SUM(obs_4.HOW_MANY) AS AB
             FROM queries.obs_4 
             GROUP BY obs_4.LOC_ID, obs_4.COUNT_YR, obs_4.TOO;
             """,
# Step 7. Remove zeroes from data, and save table to file
          """SELECT * FROM queries.obs_5 
             WHERE obs_5.AB > 0
             INTO OUTFILE '/home/kate/cbc_too_109.csv'
             FIELDS TERMINATED BY ',' 
             LINES TERMINATED BY '\n';
             """,
])               

