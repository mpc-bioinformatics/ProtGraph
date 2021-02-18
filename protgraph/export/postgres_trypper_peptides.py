import networkx
import psycopg2

from protgraph.export.abstract_exporter import AExporter


class PostgresTrypperPeptides(AExporter):
    """
    A PostGreSQL Trypper - Exporter to export PEPTIDES
    into the peptides table

    Those tables will contain all output generated by
    each of the processes. Keep in mind that this table can
    be extremely large, depending on the parmeters set in this tool.

    NOTE: Maybe even exceeding trillion of results for one protein!
    """

    def start_up(self, **kwargs):
        # Here we generate a connection to postgres
        # and generate the corresponding tables
        self.host = kwargs["postgres_trypper_host"]  # Host
        self.port = kwargs["postgres_trypper_port"]  # Port
        self.user = kwargs["postgres_trypper_user"]  # User
        self.password = kwargs["postgres_trypper_password"]  # Password
        self.database = kwargs["postgres_trypper_database"]  # Database

        # OFFSET +1 since we have dedicated start and end nodes!
        # The number of hops is also equal to the number of results generated by
        # the dynamic programming approach, counting the hops!
        self.peptide_length = kwargs["postgres_trypper_hops"] + 1  # Number of hops. E.G. 2: s -> h_1 -> h_2 -> e
        self.miscleavages = kwargs["postgres_trypper_miscleavages"]  # A filter criterion how many miscleavages?
        self.peptide_min_length = kwargs["postgres_trypper_min_pep_length"]  # Peptide minimum length

        # Initialize connection
        try:
            self.conn = psycopg2.connect(
                host=self.host,
                port=self.port,
                user=self.user,
                password=self.password,
                dbname=self.database
            )
            # Set a cursor
            self.cursor = self.conn.cursor()
        except Exception as e:
            raise Exception("Could not establish a connection to Postgres (Trypper).", e)

        # Create tables if they not exist
        try:
            self._create_tables(**kwargs)
        except Exception as e:
            raise Exception("Could not create tables in Postgres (Trypper).", e)

    def _create_tables(self, **kwargs):
        """ Create the accessions and peptides tables """
        try:
            # create accessions, so that we only save numbers in the large table!
            cur = self.conn.cursor()
            cur.execute("""
                create table if not exists accessions (
                    id SERIAl PRIMARY KEY,
                    accession VARCHAR(15) NOT NULL
                );""")
        except Exception as e:
            print("Error createing accessions table. Continuing... (Reason: {})".format(str(e)))
        finally:
            self.conn.commit()
            cur.close()

        try:
            # Create the large peptides table containing most information
            cur = self.conn.cursor()
            cur.execute("""
            CREATE TABLE  if not exists peptides (
                id BIGSERIAL UNIQUE,
                weight {0} NOT NULL,
                a_count SMALLINT NOT NULL,
                b_count SMALLINT NOT NULL,
                c_count SMALLINT NOT NULL,
                d_count SMALLINT NOT NULL,
                e_count SMALLINT NOT NULL,
                f_count SMALLINT NOT NULL,
                g_count SMALLINT NOT NULL,
                h_count SMALLINT NOT NULL,
                i_count SMALLINT NOT NULL,
                j_count SMALLINT NOT NULL,
                k_count SMALLINT NOT NULL,
                l_count SMALLINT NOT NULL,
                m_count SMALLINT NOT NULL,
                n_count SMALLINT NOT NULL,
                o_count SMALLINT NOT NULL,
                p_count SMALLINT NOT NULL,
                q_count SMALLINT NOT NULL,
                r_count SMALLINT NOT NULL,
                s_count SMALLINT NOT NULL,
                t_count SMALLINT NOT NULL,
                u_count SMALLINT NOT NULL,
                v_count SMALLINT NOT NULL,
                w_count SMALLINT NOT NULL,
                x_count SMALLINT NOT NULL,  -- NOT SKIPPED
                y_count SMALLINT NOT NULL,
                z_count SMALLINT NOT NULL,
                n_terminus character(1) NOT NULL,
                c_terminus character(1) NOT NULL,
                PRIMARY KEY (weight, a_count, b_count, c_count, 
                    d_count, e_count, f_count, g_count, h_count, 
                    i_count, j_count, k_count, l_count, m_count, 
                    n_count, o_count, p_count, q_count, r_count, 
                    s_count, t_count, u_count, v_count, w_count, 
                    x_count, y_count, z_count, n_terminus, c_terminus
                )
            );""".format("BIGINT" if kwargs["mass_dict_type"] is int else "DOUBLE PRECISION"))
        except Exception as e:
            print("Error createing peptides table. Continuing... (Reason: {})".format(str(e)))
        finally:
            self.conn.commit()
            cur.close()
            self.peptides_keys = [
                "weight",
                "a_count", "b_count", "c_count", "d_count", "e_count", "f_count", "g_count", "h_count", 
                "i_count", "j_count", "k_count", "l_count", "m_count", "n_count", "o_count", "p_count", 
                "q_count", "r_count", "s_count", "t_count", "u_count", "v_count", "w_count", "x_count", 
                "y_count", "z_count", "n_terminus", "c_terminus"
            ]
        try:
            # Create the peptides meta information (can also be extremely large)
            cur = self.conn.cursor()
            cur.execute("""
            CREATE TABLE  if not exists peptides_meta (
                id BIGSERIAL,
                peptides_id BIGINT references peptides(id),
                accession_id INT references accessions(id),
                path INT[] NOT NULL,
                miscleavages INT NOT NULL,
                PRIMARY KEY (id)
            );""")
        except Exception as e:
            print("Error createing peptides_meta table. Continuing... (Reason: {})".format(str(e)))
        finally:
            self.conn.commit()
            cur.close()
            self.peptides_meta_keys = [
                "peptides_id",
                "accession_id",
                "path",
                "miscleavages"
            ]

    def export(self, prot_graph):
        # Export the protein
        self._export(prot_graph)

    def tear_down(self):
        # Close the connection to postgres
        try:
            self.cursor.close()  # Close cursor
            self.conn.close()  # Close connection
        except Exception as e:
            print("Connection to PostgreSQL (Trypper) could not be closed. (Reason: {})".format(str(e)))

    def _export(self, prot_graph):
        # First insert accession into accession table and retrieve its id:
        accession = prot_graph.vs[0]["accession"]
        self.cursor.execute(
            "INSERT INTO accessions(accession) VALUES (%s) RETURNING id;",
            (accession,)
        )
        accession_id = self.cursor.fetchone()[0]

        # Get start and end node
        [__start_node__] = prot_graph.vs.select(aminoacid="__start__")
        [__stop_node__] = prot_graph.vs.select(aminoacid="__end__")

        # Set insert statement for peptides
        statement_peptides = " LOCK TABLE peptides IN SHARE ROW EXCLUSIVE MODE; INSERT INTO peptides (" \
            + ",".join(self.peptides_keys) \
            + ") VALUES (" \
            + ",".join(["%s"]*len(self.peptides_keys)) \
            + ") ON CONFLICT (" \
            + ",".join(self.peptides_keys) \
            + ") do update set {0} = EXCLUDED.{0} ".format(self.peptides_keys[0]) \
            + "RETURNING id;" 
        
        statement_meta_peptides = "INSERT INTO peptides_meta (" \
            + ",".join(self.peptides_meta_keys) \
            + ") VALUES (" \
            + ",".join(["%s"]*len(self.peptides_meta_keys)) \
            + ")"

        # Iterate via neworkx over all peptides (robust to parallel edges)
        netx = prot_graph.to_networkx()
        for pep in networkx.algorithms.simple_paths.all_simple_paths(
            netx,
            __start_node__.index,
            __stop_node__.index,
            cutoff=self.peptide_length
        ):
            # Get the actual Peptide (aas)
            aas = "".join(prot_graph.vs[pep[1:-1]]["aminoacid"])

            # TrypperDB specific: Filter out peptides which contain 'X'
            # if "X" in aas:
            #     continue  # TODO should we make this as a option?

            # Filter by Peptide Length
            # if len(aas) < self.peptide_min_length:
            #     continue  # Skip this entry

            # Get the weight
            edge_ids = prot_graph.get_eids(path=pep)
            if "mono_weight" in prot_graph.es[edge_ids[0]].attributes():
                weight = sum(prot_graph.es[edge_ids]["mono_weight"])
            else:
                weight = -1

            if "cleaved" in prot_graph.es[edge_ids[0]].attributes():
                cleaved = sum(filter(None, prot_graph.es[edge_ids]["cleaved"]))
            else:
                cleaved = -1

            # Filter by Miscleavages
            # if self.miscleavages != -1:
            #     if cleaved > self.miscleavages:
            #         continue  # Skip this entry

            
            peptides_tup = (
                weight,
                # Counts of Aminoacids
                aas.count("A"), aas.count("B"), aas.count("C"), aas.count("D"), aas.count("E"), aas.count("F"),
                aas.count("G"), aas.count("H"), aas.count("I"), aas.count("J"), aas.count("K"), aas.count("L"),
                aas.count("M"), aas.count("N"), aas.count("O"), aas.count("P"), aas.count("Q"), aas.count("R"),
                aas.count("S"), aas.count("T"), aas.count("U"), aas.count("V"), aas.count("W"), aas.count("X"),
                aas.count("Y"), aas.count("Z"),
                # N and C Terminus
                aas[0], aas[-1]
            )


            # Insert new entry into database:
            self.cursor.execute(statement_peptides, peptides_tup)
            peptides_id_fetched = self.cursor.fetchone()
            peptides_id = peptides_id_fetched[0]

            # Inster meta data
            peptides_meta_tup = (
                peptides_id,
                accession_id,
                pep,
                cleaved
            )
            self.cursor.execute(statement_meta_peptides, peptides_meta_tup)

        # Commit conenction
        self.conn.commit()
