"""MySQL Querying of BBS, CBC, FIA, Gentry, MCDB,and NABC for METE analyses"""

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

# BBS
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

# 3. To create table with SiteID - Year - RunType = 1:
    # from weather table
    
cursor.execute("""
                CREATE TABLE queries.weather_subquery
                SELECT (weather.statenum*1000+ weather.Route) AS SiteID,
                weather.Year, weather.RunType
                FROM BBS.weather
                WHERE weather.RunType = 1 AND weather.RPID = 101;
                """)
                 
# 4. To create table with SiteID - Year - Sp - abund:
    # Link together AOU_TOO and BBS Counts by AOU - 
    # Group By SiteID = state*1000 + route - TOO - Year 
    # Sum SpeciesTotal
    
cursor.execute("""
                CREATE TABLE queries.counts_too
                SELECT (counts.statenum * 1000) + counts.Route AS SiteID,
                counts.Year, aou_too.TOO, counts.RPID,
                SUM(counts.SpeciesTotal) AS AB 
                FROM BBS.counts INNER JOIN queries.aou_too ON 
                counts.Aou = aou_too.AOU
                GROUP BY SiteID, counts.Year, aou_too.TOO, counts.RPID
                HAVING (((counts.Year = 2009) AND (counts.RPID = 101)));
                """)
                
cursor.execute("""
                SELECT counts_too.SiteID, counts_too.Year, counts_too.TOO, 
                counts_too.AB
                FROM queries.counts_too INNER JOIN queries.weather_subquery
                ON counts_too.SiteID = weather_subquery.SiteID 
                AND counts_too.Year = weather_subquery.Year
                INTO OUTFILE '/tmp/bbs_too_2009.csv'
                FIELDS TERMINATED BY ',' 
                LINES TERMINATED BY '\n';
                """)

# CBC
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
                SELECT SUB_AUX.SUB_ID, SUB_AUX.COUNT_YR, OBS.SPECIES_CODE,
                OBS.HOW_MANY
                FROM CBC.SUB_AUX INNER JOIN CBC.OBS ON 
                SUB_AUX.SUB_ID = OBS.SUB_ID
                WHERE SUB_AUX.COUNT_YR = 109;
                """)

# Step 3. LINK obs_1 to SUB by SUB_ID to add LOC_ID 

cursor.execute("""
                CREATE TABLE queries.obs_2
                SELECT SUB.LOC_ID, obs_1.* 
                FROM queries.obs_1 INNER JOIN CBC.SUB ON 
                obs_1.SUB_ID = SUB.SUB_ID;
                """)

# Step 4. LINK obs_2 to LOC by LOC_ID to remove records WHERE SUBNATIONAL1_CODE 
    # IS NOT US-HI AND COUNTRY CODE IS 'CA', 'US', OR 'US-CA'

cursor.execute("""
                CREATE TABLE queries.obs_3
                SELECT obs_2.* 
                FROM queries.obs_2 INNER JOIN CBC.LOC ON obs_2.LOC_ID = LOC.LOC_ID
                WHERE LOC.COUNTRY_CODE = 'CA' OR LOC.COUNTRY_CODE = 'US' 
                OR LOC.COUNTRY_CODE = 'US-CA' AND 
                LOC.SUBNATIONAL1_CODE != 'US-HI';
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
                CREATE TABLE queries.obs_5                
                SELECT obs_4.LOC_ID, obs_4.COUNT_YR AS YEAR, obs_4.TOO, 
                SUM(obs_4.HOW_MANY) AS AB
                FROM queries.obs_4 
                GROUP BY obs_4.LOC_ID, obs_4.COUNT_YR, obs_4.TOO;
                """)

# Step 7. Remove zeroes from data, and save table to file

cursor.execute("""             
                SELECT * FROM queries.obs_5 
                WHERE obs_5.AB > 0
                INTO OUTFILE '/tmp/cbc_too_109.csv'
                FIELDS TERMINATED BY ',' 
                LINES TERMINATED BY '\n';
                """)

# FIA
cursor.execute("""
                CREATE TABLE queries.fia_survey1
                SELECT SURVEY.CN FROM FIA.SURVEY
                WHERE SURVEY.ANN_INVENTORY = '"Y"';
                """)

cursor.execute("""
                CREATE TABLE queries.fia_plot1
                SELECT PLOT.CN, PLOT.SRV_CN, PLOT.CTY_CN, PLOT.INVYR, 
                PLOT.STATECD, PLOT.UNITCD, PLOT.COUNTYCD, PLOT.PLOT, PLOT.LAT, 
                PLOT.LON, PLOT.ELEV FROM FIA.PLOT 
                WHERE PLOT.PLOT_STATUS_CD = 1 AND PLOT.KINDCD > 0 AND 
                PLOT.KINDCD < 4 AND PLOT.DESIGNCD = 1 OR PLOT.DESIGNCD = 311 
                OR PLOT.DESIGNCD = 312 OR PLOT.DESIGNCD = 313 OR 
                PLOT.DESIGNCD = 314 OR PLOT.DESIGNCD = 328 OR 
                PLOT.DESIGNCD = 220 OR PLOT.DESIGNCD = 240 AND 
                PLOT.MANUAL >= 1 AND PLOT.QA_STATUS = 1 AND 
                PLOT.SAMP_METHOD_CD = 1 AND PLOT.INVYR < 3000; 
                """)

cursor.execute("""
                CREATE TABLE queries.fia_plot2
                SELECT fia_plot1.CN, fia_plot1.CTY_CN, fia_plot1.INVYR,
                fia_plot1.STATECD, fia_plot1.UNITCD, fia_plot1.COUNTYCD, 
                fia_plot1.PLOT, fia_plot1.LAT, fia_plot1.LON, fia_plot1.ELEV
                FROM queries.fia_plot1 INNER JOIN queries.fia_survey1 
                ON fia_plot1.SRV_CN = fia_survey1.CN 
                ORDER BY fia_plot1.STATECD, fia_plot1.UNITCD, 
                fia_plot1.COUNTYCD, fia_plot1.PLOT, fia_plot1.INVYR DESC;
                """)

cursor.execute("""
                CREATE TABLE queries.fia_cond1
                SELECT COND.PLT_CN FROM FIA.COND 
                WHERE (((COND.STDORGCD)="0" OR (COND.STDORGCD)='') 
                AND ((COND.TRTCD1) = '' OR (COND.TRTCD1)="0") 
                AND ((COND.TRTCD2) = '' OR (COND.TRTCD2)="0") 
                AND ((COND.TRTCD3) = '' OR (COND.TRTCD3)="0"))
                GROUP BY COND.PLT_CN;
                """)

cursor.execute("""
                CREATE TABLE queries.fia_plot3
                SELECT fia_plot2.CN, fia_plot2.STATECD, fia_plot2.UNITCD,
                fia_plot2.COUNTYCD, fia_plot2.PLOT, 
                MAX(fia_plot2.INVYR) AS INVYR, AVG(fia_plot2.LAT) AS LAT,
                AVG(fia_plot2.LON) AS LON, AVG(fia_plot2.ELEV) AS ELEV
                FROM queries.fia_plot2
                GROUP BY fia_plot2.STATECD, fia_plot2.UNITCD, 
                fia_plot2.COUNTYCD, fia_plot2.PLOT;
                """)

cursor.execute("""
                CREATE TABLE queries.fia_plot4
                SELECT fia_plot3.* FROM queries.fia_plot3
                INNER JOIN queries.fia_cond1 ON fia_plot3.CN = fia_cond1.PLT_CN; 
                """)

cursor.execute("""
                CREATE TABLE queries.fia_tree1
                SELECT TREE.PLT_CN, TREE.STATECD, TREE.UNITCD, TREE.COUNTYCD, 
                TREE.PLOT, TREE.SPCD FROM FIA.TREE
                WHERE TREE.STATUSCD = 1;
                """)

cursor.execute("""
                CREATE TABLE queries.fia_tree2
                SELECT ((fia_tree1.STATECD*10000000000) + 
                (fia_tree1.UNITCD*1000000000) + (fia_tree1.COUNTYCD*1000000) + 
                fia_tree1.PLOT) AS PlotID, fia_tree1.PLT_CN, fia_tree1.SPCD
                FROM queries.fia_tree1
                INNER JOIN queries.fia_plot4 ON fia_tree1.PLT_CN = fia_plot4.CN 
                WHERE fia_plot4.INVYR < 3000;
                """)

cursor.execute("""
                CREATE TABLE queries.fia_tree3
                SELECT fia_tree2.PlotID, fia_tree2.PLT_CN, fia_tree2.SPCD,
                COUNT(fia_tree2.SPCD) AS AB
                FROM queries.fia_tree2 
                GROUP BY fia_tree2.PlotID, fia_tree2.PLT_CN, fia_tree2.SPCD;
                """)

cursor.execute("""
                SELECT * FROM queries.fia_tree3
                INTO OUTFILE '/tmp/fia_spab.csv'
                FIELDS TERMINATED BY ',' 
                LINES TERMINATED BY '\n';
                """)

# Gentry
# Group By SITE and SPECIES_ID, COUNT stem_id
cursor.execute("""
                CREATE TABLE queries.gentry
                SELECT Stems.site, Stems.species_id, Stems.Morpho, 
                Stems.id_status, sum(`Stems`.`N(Ind)`) AS ab 
                FROM gentry_glenda.Stems 
                GROUP BY Stems.site, Stems.species_id 
                HAVING (((ab > 0) AND (Stems.Morpho = 1 OR 
                Stems.id_status = "species")));
                """)
                
cursor.execute("""
                SELECT gentry.site, gentry.species_id, gentry.ab 
                FROM queries.gentry
                INTO OUTFILE '/tmp/gentry_spab.csv'
                FIELDS TERMINATED BY ',' 
                LINES TERMINATED BY '\n';
                """)

# MCDB
# Select out raw abundance data
cursor.execute("""
                CREATE TABLE queries.mcdb1
                SELECT communities.Site_ID AS site,
                communities.Initial_year AS year,
                communities.Species_ID AS sp, communities.Abundance AS ab,
                sites.Abundance_data_format AS format 
                FROM MCDB.communities
                INNER JOIN MCDB.sites USING (Site_ID)
                HAVING (((ab > 0) AND (format = 'raw')));
                """)

# Create intermediate table that has the value of the earliest year of 
# sampling for each site 
cursor.execute("""
                CREATE TABLE queries.mcdb2
                SELECT mcdb1.site, Min(mcdb1.year) AS year
                FROM queries.mcdb1
                GROUP BY mcdb1.site;
                """)

# Use intermediate table to select out only one year of data per site
cursor.execute("""
                CREATE TABLE queries.mcdb3
                SELECT mcdb1.site, mcdb1.year, mcdb1.sp, mcdb1.ab 
                FROM queries.mcdb1
                INNER JOIN queries.mcdb2 USING (site, year);
                """)

# Dump into csv file
cursor.execute("""
                SELECT mcdb3.* FROM queries.mcdb3
                INTO OUTFILE '/tmp/mcdb_spab.csv'
                FIELDS TERMINATED BY ',' 
                LINES TERMINATED BY '\n';
                """)

# NABC
cursor.execute("""
                CREATE TABLE queries.nabc_sample_dates 
                SELECT NABA_2009.Count_State, NABA_2009.Count_Name, NABA_2009.Date
                FROM NABC.NABA_2009
                GROUP BY NABA_2009.Count_State, NABA_2009.Count_Name,
                NABA_2009.Date;
                """)

cursor.execute("""
                CREATE TABLE queries.nabc_min_sample_date
                SELECT nabc_sample_dates.Count_State,
                nabc_sample_dates.Count_Name, 
                MIN(nabc_sample_dates.Date) AS MinDate
                FROM queries.nabc_sample_dates
                GROUP BY nabc_sample_dates.Count_State, 
                nabc_sample_dates.Count_Name;
                """)

cursor.execute("""
                CREATE TABLE queries.nabc_sp_ab_2009a
                SELECT NABA_2009.Count_State, NABA_2009.Count_Name, 
                CONCAT(NABA_2009.Count_State,"_",Left(NABA_2009.Count_Name,5),
                "_",Right(NABA_2009.Count_Name,5)) AS SiteID,
                YEAR(NABA_2009.Date) AS Year, NABA_2009.Scientific_Name, 
                Sum(NABA_2009.Number_Butterflies) AS AB
                FROM NABC.NABA_2009 INNER JOIN queries.nabc_min_sample_date ON 
                (nabc_min_sample_date.Count_State = NABA_2009.Count_State) AND 
                (nabc_min_sample_date.Count_Name = NABA_2009.Count_Name) AND 
                (nabc_min_sample_date.MinDate = NABA_2009.Date)
                GROUP BY NABA_2009.Count_State, NABA_2009.Count_Name,
                NABA_2009.Date, NABA_2009.Scientific_Name;
                """)

cursor.execute("""
                CREATE TABLE queries.nabc_sp_ab_2009
                SELECT nabc_sp_ab_2009a.SiteID, nabc_sp_ab_2009a.Year,
                CONCAT(NABA_species.Genus, "_", NABA_species.Species) AS SpID, 
                SUM(nabc_sp_ab_2009a.AB) AS Abund
                FROM queries.nabc_sp_ab_2009a INNER JOIN NABC.NABA_species ON 
                NABA_species.Scientific_Name = nabc_sp_ab_2009a.Scientific_Name 
                GROUP BY nabc_sp_ab_2009a.SiteID, NABA_species.Genus, 
                NABA_species.Species;
                """)

# Dump into csv file, removing two sites that include significant 
# outliers in abundance, one being the NABA Butterfly park
cursor.execute("""
                SELECT nabc_sp_ab_2009.* FROM queries.nabc_sp_ab_2009
                WHERE nabc_sp_ab_2009.SiteID != "TX_NABA _ Park" AND 
                nabc_sp_ab_2009.SiteID != "MN_Bear _ction"
                INTO OUTFILE '/tmp/nabc_spab.csv'
                FIELDS TERMINATED BY ',' 
                LINES TERMINATED BY '\n';
                """)

connection.commit()

shutil.copy('/tmp/bbs_too_2009.csv', '/home/kate/data/bbs_too_2009.csv')
shutil.copy('/tmp/cbc_too_109.csv', '/home/kate/data/cbc_too_109.csv')
shutil.copy('/tmp/fia_spab.csv', '/home/kate/data/fia_spab.csv')
shutil.copy('/tmp/gentry_spab.csv', '/home/kate/data/gentry_spab.csv')
shutil.copy('/tmp/mcdb_spab.csv', '/home/kate/data/mcdb_spab.csv')
shutil.copy('/tmp/nabc_spab.csv', '/home/kate/data/nabc_spab.csv')