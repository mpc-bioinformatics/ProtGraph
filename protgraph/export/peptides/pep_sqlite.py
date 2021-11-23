

import sqlite3
import apsw

import os
import zlib
from protgraph.export.peptides.abstract_peptide_exporter import \
    APeptideExporter

from protgraph.export.peptides.pep_fasta import PepFasta


class PepSQLite(APeptideExporter):
    """
    A PostGreSQL - Exporter to export PEPTIDES
    into the peptides table

    Those tables will contain all output generated by
    each of the processes. Keep in mind that this table can
    be extremely large, depending on the parmeters set in this tool.

    NOTE: Maybe even exceeding trillion of results for one protein!
    """

    def start_up(self, **kwargs):

        # Traversal parameters:
        self._set_up_taversal(
            kwargs["pep_sqlite_skip_x"],
            kwargs["pep_sqlite_min_pep_length"],
            kwargs["pep_sqlite_miscleavages"],
            kwargs["pep_sqlite_use_igraph"],
            kwargs["pep_sqlite_hops"],
            kwargs["pep_sqlite_batch_size"]
        )

        self.compact_repr = not kwargs["pep_sqlite_non_compact"]

        
        # Create database file
        os.makedirs(kwargs["pep_sqlite_output_dir"], exist_ok=True)
        # self.conn = sqlite3.connect(kwargs["pep_sqlite_database"], timeout=10, isolation_level=None)
        self.conn = apsw.Connection(os.path.join(kwargs["pep_sqlite_output_dir"], kwargs["pep_sqlite_database"]))
        self.conn.setbusytimeout(10000)

        self.pepfasta = PepFasta()


        # Set sqlite specific parameters
        self.use_blob = kwargs["pep_sqlite_use_blobs"]
        self.compression_level = kwargs["pep_sqlite_compression_level"]

        # Create tables if they not exist
        try:
            self._create_tables(**kwargs)
        except Exception as e:
            raise Exception("Could not set up sqlite.", e)

    def _create_tables(self, **kwargs):
        """ Set Pragmas """
        try:
            # Create the large peptides table containing most information
            cur = self.conn.cursor()
            cur.execute("""
            PRAGMA journal_mode=OFF
            """
            )
        except Exception as e:
            print("Warning: Pragma journal_mode was not set, inserting may be slower (Reason: {})".format(str(e)))
        finally:
            # self.conn.commit()
            cur.close()

        try:
            # Create the large peptides table containing most information
            cur = self.conn.cursor()
            cur.execute("""
            PRAGMA synchronous=OFF
            """
            )
        except Exception as e:
            print("Warning: pragma synchronous was not set, inserting may be slower (Reason: {})".format(str(e)))
        finally:
            # self.conn.commit()
            cur.close()

        try:
            # Create the large peptides table containing most information
            cur = self.conn.cursor()
            cur.execute("""
            PRAGMA temp_store = MEMORY;
            """
            )
        except Exception as e:
            print("Warning: pragma synchronous was not set, inserting may be slower (Reason: {})".format(str(e)))
        finally:
            # self.conn.commit()
            cur.close()

        try:
            # Create the large peptides table containing most information
            cur = self.conn.cursor()
            cur.execute("""
            PRAGMA mmap_size=268435456;
            """
            )
        except Exception as e:
            print("Warning: pragma synchronous was not set, inserting may be slower (Reason: {})".format(str(e)))
        finally:
            # self.conn.commit()
            cur.close()


        try:
            # Create the large peptides table containing most information
            cur = self.conn.cursor()
            cur.execute("""
            PRAGMA cache_size = 1000000;
            """
            )
        except Exception as e:
            print("Warning: pragma synchronous was not set, inserting may be slower (Reason: {})".format(str(e)))
        finally:
            # self.conn.commit()
            cur.close()


        
        try:
            # Create the peptides meta information (can also be extremely large), larger than the peptides tables
            cur = self.conn.cursor()
            if self.compact_repr:
                cur.execute("""
                CREATE TABLE if not exists peptide_meta (
                    peptide {B} PRIMARY KEY,
                    meta {B});""".format(B="BLOB" if self.use_blob else "TEXT")
                )
            else: 
                cur.execute("""
                CREATE TABLE if not exists peptide_meta (
                    peptide {B},
                    meta {B});""".format(B="BLOB" if self.use_blob else "TEXT")
                )
        except Exception as e:
            print("Warning: Failed creating table 'peptides_meta' (Reason: {})".format(str(e)))
        finally:
            # self.conn.commit()
            cur.close()

        if self.use_blob:
            self.set_upsert_query = \
                """
                INSERT INTO peptide_meta (peptide, meta)
                VALUES (?, ?)
                ON CONFLICT(peptide) DO UPDATE SET meta = meta || ? ;
                """
        else:
            self.set_upsert_query = \
                """
                INSERT INTO peptide_meta (peptide, meta)
                VALUES (?, ?)
                ON CONFLICT(peptide) DO UPDATE SET meta = meta || ',' || ? ;
                """

    def export(self, prot_graph, queue):
        # We continue with the export function
        super().export(prot_graph, queue)

        # and commit everything in the conenction for a protein
        # self.conn.commit()

    def export_peptides(self, prot_graph, l_path_nodes, l_path_edges, l_peptide, l_miscleavages, _):
        # Retrieve Meta Infos
        accs = [self.pepfasta._get_accession_or_isoform(prot_graph.vs[nodes[1]]) for nodes in l_path_nodes]
        start_poses = [str(self.pepfasta._get_position_or_isoform_position(prot_graph.vs[nodes[1]])) for nodes in l_path_nodes]
        end_poses = [str(self.pepfasta._get_position_or_isoform_position(prot_graph.vs[nodes[-2]], end=True)) for nodes in l_path_nodes]
        qualifiers =  [",".join(self.pepfasta._get_qualifiers(prot_graph, edges)) for edges in l_path_edges]
        qualifiers = ["," + quali if quali else "" for quali in qualifiers ]
        metas = [
            "".join(
                [
                    acc, "(", str(start_pos), ":", str(end_pos), ",",
                    "mssclvg:", str(misses), quali_entries,  ")"
                ])
            for acc, start_pos, end_pos, quali_entries, misses in zip(accs, start_poses, end_poses, qualifiers, l_miscleavages)
        ]

        if self.use_blob:
            compressed_metas = [zlib.compress(x.encode(), level=self.compression_level) for x in metas]
            compressed_peptides = [zlib.compress(x.encode(), level=self.compression_level) for x in l_peptide]
        else:
            compressed_metas = metas
            compressed_peptides = l_peptide


        # Insert into database
        cur = self.conn.cursor()

        # executemany
        cur.execute('BEGIN TRANSACTION')
        if self.compact_repr:
            self._retry_query_many(
                cur, 
                self.set_upsert_query,
                zip(compressed_peptides, compressed_metas, compressed_metas)
            )
        else:
            self._retry_query_many(
                cur, 
                """
                INSERT INTO peptides_meta (peptide, meta)
                VALUES (?, ?);
                """,
                zip(compressed_peptides, compressed_metas)
            )
        cur.execute('COMMIT')




        # Getting blob entries:
        # cur.execute("""
        # Select ROWID from peptides_meta where peptide = x'78da0b020000530053';
        # """
        # )
        # print("=========================================================")
        # a = cur.fetchone()

        # bl = self.conn.blobopen("main", "peptides_meta", "meta", a[0], 1)
        # binar = bl.read()
        # bl.close()

        # stream = binar
        # while stream:
        #     dco = zlib.decompressobj()
        #     print(dco.decompress(stream))
        #     stream = dco.unused_data


        # self.conn.commit()
        cur.close()

    def _retry_query_many(self, cursor, statement, entries):
        # Execute statement. Retry if failed.
        try:
            cursor.executemany(statement, entries)
        except Exception:
            self._retry_query_many(cursor, statement, entries)

    def tear_down(self):
        # Close the connection to sqlite
        try:
            self.conn.close()  # Close connection
        except Exception as e:
            print("Connection to Sqlite  could not be closed. (Reason: {})".format(str(e)))
