import MySQLdb as dbapi
import getpass
import shutil

def get_data(queries):
    """Function to connect to database, create query tables, and output to CSV.
    """

    # enter MySQL password
    p=getpass.getpass()

    connection = dbapi.connect(host='129.123.92.244', port=1977, user='kate', 
                           passwd=p)
    cursor = connection.cursor()

    cursor.execute("""DROP DATABASE IF EXISTS queries;""")
    cursor.execute("""CREATE DATABASE queries;""")
    
    for query in queries:
        cursor.execute(query)

    connection.commit()
