import os
import hashlib
import logging
import pandas as pd
import importlib.resources as resources
import sqlite3 as sl
import tempfile

__version__ = "1.0.0"

SQL = os.path.join(os.path.dirname(__file__),"pcalf_db_scheme.sql")

class Pcalf_DB_Error(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return 'Pcalf_DB_Error, {0} '.format(self.message)
        else:
            return 'Pcalf_DB_Error has been raised'


class PcalfDB:

    def __init__(self,dbfile):        
        self.dbfile = self._parse_file(dbfile)
        self.md5 = self._digest_scheme()
        self.vmd5 = self._current_md5()
        if not self.is_schema_valid():
            logging.warning("{} version is different from current DB schema.".format(dbfile))

    def _parse_file(self, file):
        if os.path.exists(file):
            return self._load(file)    
        return self._create(file)
    
    def _current_md5(self):
        tmp = tempfile.NamedTemporaryFile(mode="w")          
        f = self._create(tmp.name)
        tmp.flush()
        con = sl.connect(f)
        cur = con.cursor()
        ex  = cur.execute("SELECT * FROM sqlite_master")
        exstr = str(ex.fetchall())
        return hashlib.md5(exstr.encode()).hexdigest()         


    def _load(self,file):
        return file
        
    def _create(self,file):     
        if os.path.exists(file) and os.stat(file).st_size != 0:
            raise Pcalf_DB_Error("{} file exits and doesn't seem empty".format(file))
        con,cur = self.connect(file)
        with open(SQL) as fp: 
            cur.executescript(fp.read()) # or con.executescript
        return file

    def is_schema_valid(self):
        return self.md5 == self.vmd5


    def connect(self,file=None):
        if file is None:
            file = self.dbfile
        con = sl.connect(file)
        cur = con.cursor()
        return con,cur

    def is_same_schema(self, db):
        assert isinstance(db,PcalfDB)
        #assert os.path.exists(self.dbfile) and os.path.exists(db)        
        if self.md5 != db.md5:
            return False
        return True
    
    def _digest_scheme(self):
        con,cur  = self.connect()        
        ex  = cur.execute("SELECT * FROM sqlite_master")
        exstr = str(ex.fetchall())
        return hashlib.md5(exstr.encode()).hexdigest()
        
    def is_empty(self,table):
        con,cur = self.connect()
        if self.check_if_table_exist(cur,table):    
            cur.execute('''SELECT COUNT(*) from {} '''.format(table))
            result=cur.fetchall()
            return result[0][0]==0
        raise  Pcalf_DB_Error("{} doesn't exist in {}".format(table,self.dbfile))
        
    def list_table(self):
        con,cur = self.connect()
        cur.execute(''' SELECT name FROM sqlite_master WHERE type='table' ''')
        return [r[0] for r in cur.fetchall()]

    def list_columns(self,table):
        con, cur = self.connect()
        if self.check_if_table_exist(cur,table):    
            cur.execute('''PRAGMA table_info({});'''.format(
                table
            ))        
            return [r[1] for r in cur.fetchall()]
        raise  Pcalf_DB_Error("{} doesn't exist in {}".format(table,self.dbfile))

    
    def to_df(self,table):
       con,_ = self.connect()
       return pd.read_sql_query("SELECT * from {}".format(table), con)
    
    def to_fasta(self, table ):
        if table in ["gly1", "gly2" , "gly3" , "glyx3"]:
            con, cur = self.connect()
            if self.check_if_table_exist(cur,table):
                cur.execute(''' SELECT `sequence_id`, `sequence` FROM {} '''.format(table))
                #result = cur.fetchall()
                return {r[0]:r[1] for r in cur.fetchall()}
            raise  Pcalf_DB_Error("{} doesn't exist in {}".format(table,self.dbfile))
        raise Pcalf_DB_Error("{} can't be converted into fasta dictionnary.".format(table))
    
    def to_msa(self,table):
        msa_dict = self.to_fasta(table)
        if msa_dict:
            length = len(list(msa_dict.values())[0])
            for seq in msa_dict.values():
                if len(seq) != length:
                    raise Pcalf_DB_Error("Sequence in a MSA should have the same length [{} != {}].".format(length,len(seq)))
            return msa_dict
        raise Pcalf_DB_Error("{} table from {} seems empty and can't be converted into fasta dictionnary.".format(table,self.dbfile))
    
    def generate_nter_db(self):
        rows = []
        con, cur = self.connect()
        if self.check_if_table_exist(cur,"summary") and self.check_if_table_exist(cur,"features"):            
            cur.execute('''SELECT f.sequence_id, 
                    f.feature_src,
                    f.feature_seq
                FROM features as f JOIN summary as s 
                ON f.sequence_id = s.sequence_accession 
                WHERE f.feature_id = "N-ter" AND s.flag = "Calcyanin with known N-ter"''')
        
            for r in cur.fetchall():
                nter,_ = r[1].split("||")
                rows.append((r[0],nter,r[-1]))
        return rows

    def get_col_values(self,table,key="pk"):    
        con,cur = self.connect()
        if self.check_if_table_exist(cur,table):
            with con:
                cur.execute('''
                    SELECT `{col}` FROM {table}
                '''.format(
                    table=table,
                    col=key,            
                ))
                rows = cur.fetchall()
                vals = [row[0] for row in rows]               
            return vals
        else:
            return []

    def check_if_table_exist(self,cur,table):
        #get the count of tables with the name        
        cur.execute('''
            SELECT count(name) FROM sqlite_master WHERE type='table' AND name='{}' 
        '''.format(table))
        #if the count is 1, then table exists
        if cur.fetchone()[0]==1 : 
            return True
        return False
    
    def create_from_df(self,dataframe,table,pk):
        con,cur = self.connect()
        logging.info("Create {} with {} as primary key".format(table,pk))
        return dataframe.to_sql(name=table,con=con,  if_exists = 'fail' , index_label = pk )

    def remove_duplicate(self,dataframe,con,cur,table,pk):
        with con:
            cur.execute('''
                SELECT `{pk}` FROM {table}
            '''.format(pk=pk,table=table))
            rows = cur.fetchall()        
        indb = [row[0] for row in rows]        
        return dataframe[~dataframe[pk].isin(indb)]

    def update_table(self,dataframe,table,pk):
        con,cur = self.connect()       
        
        filtered_df = self.remove_duplicate(dataframe,con,cur,table,pk)
        logging.debug("dataframe shape before and after fitlering : {} -> {}".format(
            dataframe.shape , filtered_df.shape))
        logging.info("Update {} with {} as primary key".format(table,pk))
        if not filtered_df.empty :
            filtered_df.set_index(pk,inplace=True)
            return filtered_df.to_sql(name=table, con=con,  if_exists = 'append')#, index=True , index_label = pk )
        else:
            logging.warning("Nothing to add after redundancy removal.")
            return None

    def feed_db(self,df,table,pk):
        con,cur = self.connect()

        sqlite_cols = self.list_columns(table)
        df_cols = df.columns
        if (sorted(sqlite_cols) != sorted(df_cols)):
            raise Pcalf_DB_Error("dataframe columns doesn't match sqlite table columns. Get {} [Expected : {}]".format(
                ", ".join(sorted(df_cols)) , ", ".join(sorted(sqlite_cols)) 
            ))

        if self.check_if_table_exist(cur,table):            
            self.update_table(df,table,pk)
        else:            
            self.create_from_df(df,table,pk)


    def _fasta_to_tsv(self,fasta):
        rows = []
        row = []
        with open(fasta) as fh:            
            for line in fh.readlines():
                if line.startswith(">"):
                    if row:
                        rows.append(row)
                    row = [ line.split()[0].replace(">","") , "" ]
                else:                
                    if row:
                        row[-1]+=line.strip()
                    else:
                        continue
        print(len(rows[0][-1]))
        rows.append(row)
        return pd.DataFrame(rows,columns=["sequence_id","sequence"])
        

    def feed_msa_db(self,msa,table,pk):
        if os.path.exists(msa):
            msa_table = self._fasta_to_tsv(msa)
        else:
            raise Pcalf_DB_Error("MSA [{}] doesn't exists ... ".format(msa))
        self.feed_db(
            msa_table,
            table,
            pk
        )


