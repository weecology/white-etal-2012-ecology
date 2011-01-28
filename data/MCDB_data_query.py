"""MySQL Querying of MCDB to yield data for METE analyses"""

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

# Select out raw abundance data
cursor.execute("""
                CREATE TABLE queries.mcdb1
                SELECT community_data.Site_ID AS site, community_data.Initial_year AS year,
                community_data.Species_ID AS sp, community_data.Abundance AS ab,
                sites.Abundance_data_format AS format 
                FROM MCDB.community_data
                INNER JOIN MCDB.sites USING (Site_ID)
                HAVING (((ab > 0) AND (format = 'raw')));
                """)

# Create intermediate table that has the value of the earliest year of sampling for each site 
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
               
connection.commit()

shutil.copy('/tmp/mcdb_spab.csv', '/home/kate/data/mcdb_spab.csv')