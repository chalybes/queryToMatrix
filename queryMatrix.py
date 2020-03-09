# Need to download sqlalchemy-redshift, and psycopg2 (PostgreSQL)
# if you haven't done so in order for this notebook to work
import pandas as pd
import numpy
import scipy.io
import scipy.sparse
# import psycopg2
import pyreadr
import csv
from pathlib import Path
from sqlalchemy import create_engine
import sensitive

# Redshift Postgres connection string goes below
connstr = redshiftlink

engine = create_engine(connstr)

'''Currently this script takes in an RObject that's a list of sample names with ARID. The script will then
    query Redshift for the transcript counts of these samples and in turn create a sparse matrix from the
    count matrix.'''


def getQuery(rdafile, tablename):

    getQuery.tablename = tablename

    if rdafile:
        loadRDA = str(rdafile)
        pyrfile = pyreadr.read_r(loadRDA)
        getQuery.rdalist = list(pyrfile.values())
        rdadf = pd.DataFrame(data=getQuery.rdalist[0])
        getQuery.rdalist = rdadf['sampled.cells']

        if len(getQuery.rdalist) > 100000:
            getQuery.query_chunks = numpy.array_split(getQuery.rdalist, (len(getQuery.rdalist)/100000))
            print(type(getQuery.query_chunks))

    return getQuery


executeQ = getQuery(rdafile=allen_internal, tablename="arid_counts_table_10x")


def makeconnection(gQ):

    # dfcollection = []
    with engine.connect() as conn, conn.begin():
        if not gQ.query_chunks:
            placeholders = ", ".join(["%s" for _ in getQuery.rdalist])
            filteredQuery = "SELECT * FROM " + gQ.tablename + " WHERE sample_id IN ({});".format(placeholders)
            df = pd.read_sql(filteredQuery, conn, params=gQ.rdalist)

            # Reshape the format from long to wide
            df = df.reset_index().pivot(index='sample_id', columns='gene', values='counts')
        else:
            for part in gQ.query_chunks:
                placeholders = ", ".join(["%s" for _ in part])
                filteredQuery = "SELECT * FROM " + gQ.tablename + " WHERE sample_id IN ({});".format(placeholders)
                dfpart = pd.read_sql(filteredQuery, conn, params=part)
                # dfpart = dfpart.reset_index().pivot(index='sample_id', columns='gene', values='counts')
                df = pd.DataFrame()
                df = df.append(dfpart)

        # df = pd.concat(dfcollection, axis=0, sort=False, ignore_index=True)
        # del dfcollection
        df = df.pivot(index='sample_id', columns='gene', values='counts')
    return df


DFtoSave = makeconnection(executeQ)


# dbob is a function call of makeconnection()
def outofile(dfob, outname, fiformat):
    # if fiformat == "csv":
    #     # Output to CSV
    #     outname = str(outname) + ".csv"
    #     # currently separating by commas, can use sep='\t'
    #     dfob.to_csv(outname, sep=',', index=True)

    print("Reached file output stage!")
    if fiformat == "mtx":
        # The columns and rows are actually flipped
        rowNamesArr = list(dfob.columns.values)
        colNamesArr = list(dfob.index.values)

        with open('batch_row_index.tsv', 'w', newline='\n') as row_output:
            tsv_col = csv.writer(row_output)
            tsv_col.writerow(rowNamesArr)
        with open('batch_col_index.tsv', 'w', newline='\n') as col_output:
            tsv_row = csv.writer(col_output)
            tsv_row.writerow(colNamesArr)
        print("Reached making sparse matrix!")
        # Output to matrix format
        outname = str(outname) + ".mtx"
        asparsemat = scipy.sparse.csr_matrix(dfob.values)
        scipy.io.mmwrite(outname, asparsemat)


outofile(DFtoSave, 'testOutput', 'mtx')
