import igraph
import sqlalchemy
import json

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from sqlalchemy import Column, BigInteger, String, JSON
from sqlalchemy.ext.declarative import declarative_base
Base = declarative_base()


import pyorient 


import redis
from redisgraph import Node, Edge, Graph, Path


import psycopg2


# todo check if graph is dag and only has one starting and end pooint (via degree)

# class Node(Base):
#     __tablename__ = "nodes"
#     id = Column(BigInteger, primary_key=True)
#     accession = Column(String)
#     aminoacid = Column(String)
#     attributes = Column(JSON)

#     def __repr__(self):
#        return "<Node(accession='{}', aminoacid='{}', attributes='{}')>".format(self.accession, self.aminoacid, self.attributes)

# class Edge(Base):
#     __tablename__ = "edges"
#     id = Column(BigInteger, primary_key=True)
#     accession = Column(String)
#     aminoacid = Column(String)
#     attributes = Column(JSON)

#     def __repr__(self):
#        return "<Node(accession='{}', aminoacid='{}', attributes='{}')>".format(self.accession, self.aminoacid, self.attributes)




def _get_qualifiers(list_of_qs):
    if list_of_qs is None:
        return None
    elif len(list_of_qs) == 0:
        return None
    else:
        l = []
        for k in list_of_qs:
            d = k.qualifiers
            d["type"] = k.type
            d["id"] = k.id
            l.append(d)
        return l


def insert_to_database(digested_graphs_queue, output_query):




    # Version for Redis Graph Local installation

    # r = redis.Redis(host="localhost", port=6379)

    # keys1 = ["accession", "aminoacid", "avrg_end_weight", "mono_end_weight"]
    # keys2 = ["avrg_end_weight", "avrg_weight", "mono_end_weight", "mono_weight"]


    # while True:
    #     try: 
    #         graph_entry = digested_graphs_queue.get(timeout=180)
    #     except Exception:
    #         continue

    #     redis_graph = Graph("proteins", r)


 

    #     # redis_nodes = [
    #     #     Node( label = "node", properties = x.attributes() ) for x in graph_entry.vs[:]
    #     # ]
    #     redis_nodes = [ # Excluding most attributes due to having None inside of them (causing errors)
    #         Node( label = "node", properties = {key:value for key, value in x.attributes().items() if key in keys1} ) for x in graph_entry.vs[:]
    #     ]
    #     for x in redis_nodes:
    #         redis_graph.add_node(x)
        
    #     # redis_edges = [
    #     #     Edge( redis_nodes[x.source], "edge", redis_nodes[x.target], properties = x.attributes() )  for x in graph_entry.es[:]
    #     # ]
    #     redis_edges = [
    #         Edge( redis_nodes[x.source], "edge", redis_nodes[x.target], properties = {key:value for key, value in x.attributes().items() if key in keys2} )  for x in graph_entry.es[:]
    #     ]
    #     for x in redis_edges:
    #         redis_graph.add_edge(x)

    #     p = redis_graph.commit()

    #     output_query.put(1)






    ### BATCH VERSION OF  Orient DB
    # while True:
    #     try: 
    #         graph_entry = digested_graphs_queue.get(timeout=180)
    #     except Exception:
    #         continue


    #     # Batch Version
    #     client = pyorient.OrientDB("localhost", 2424)
    #     client.set_session_token(True)
    #     session_id = client.connect("root", "root")
    #     client.db_open("graph", "root", "root")

        
    #     batch_cmd = ["begin;\n"]
    #     batch_main_nodes = ["let e{} = create vertex V content {};\n".format(x, json.dumps(y.attributes())) for x, y in enumerate(graph_entry.vs[:])]
        

    #     nodes_ids = ["$e" + str(x) for x in range(graph_entry.vcount())]
    #     sources_nodes = [nodes_ids[x.source] for x in graph_entry.es[:]]
    #     target_nodes  = [nodes_ids[x.target] for x in graph_entry.es[:]]
    #     content = [x.attributes() for x in graph_entry.es[:]]
    #     content_updated_qualifiers = [_get_qualifiers(c["qualifiers"]) for c in content]
    #     for x, y in zip(content, content_updated_qualifiers):
    #         x["qualifiers"] = y
        
    #     batch_main_edges = ["create edge E from {} to {} content {};".format(x, y, json.dumps(z)) for x, y, z in zip(sources_nodes, target_nodes, content)]
    #     batch_end = ["commit retry 100;\n"]

    #     results = client.batch( "".join(batch_cmd + batch_main_nodes + batch_main_edges + batch_end ))


    #     output_query.put(1)
        # To get Information of how many edges and Nodes are in the database
        #
        # Select ( SELECT COUNT(*) FROM V ) AS count1, ( SELECT COUNT(*) FROM E ) AS count2




        ## Non Batch Version
        # client = pyorient.OrientDB("localhost", 2424)
        # client.set_session_token(True)
        # session_id = client.connect("root", "root")
        # client.db_open("graph", "root", "root")

        
        # node_objects = [client.command("insert into V content " +  json.dumps(x.attributes())) for x in graph_entry.vs[:]]



        # sources_nodes = [node_objects[x.source][0]._rid for x in graph_entry.es[:]]
        # target_nodes  = [node_objects[x.target][0]._rid for x in graph_entry.es[:]]
        # content = [x.attributes() for x in graph_entry.es[:]]
        # content_updated_qualifiers = [_get_qualifiers(c["qualifiers"]) for c in content]
        # for x, y in zip(content, content_updated_qualifiers):
        #     x["qualifiers"] = y

        # edge_objects = [client.command("create edge E from {} to {} content {}".format(x, y, json.dumps(z))) for x, y, z in zip(sources_nodes, target_nodes, content)]







    ### psycopg2 
    conn = psycopg2.connect(host="127.0.0.1", port=5433, user="postgres", password="developer", dbname="node_edges")
    cur = conn.cursor()


    # create nodes
    cur.execute("""
        create table if not exists nodes (
            id BIGSERIAL PRIMARY KEY,
            accession TEXT NOT NULL,
            aminoacid VARCHAR(10) NOT NULL,
            attributes JSON NOT NULL
        );""")
    conn.commit()


    # Create edges
    cur.execute("""
        create table if not exists edges (
            id BIGSERIAL PRIMARY KEY,
            source BIGINT references nodes(id),
            target BIGINT references nodes(id),
            mono_weight BIGINT NOT NULL,  
            mono_end_weight BIGINT NOT NULL, 
            avrg_weight BIGINT NOT NULL,  
            avrg_end_weight BIGINT NOT NULL, 
            attributes JSON NOT NULL
        );""")
    conn.commit()


    while True:
        try: 
            graph_entry, num_of_paths = digested_graphs_queue.get(timeout=180)
        except Exception:
            continue

        # Generate Key list for json
        json_keys = [x for x in graph_entry.vs[0].attribute_names() if not( x == "accession" or x == "aminoacid" )]

        # Generate Node Object
        db_nodes = [(x["accession"], x["aminoacid"], json.dumps({y: x[y] for y in json_keys})) for x in graph_entry.vs[:]]

        # Insert all nodes (bulk)
        
        
        
        
        
        
        statement = "INSERT INTO nodes(accession, aminoacid, attributes) VALUES " + ",".join(["(%s,%s,%s)"]*len(db_nodes)) + " RETURNING id"
        db_nodes_falttened = [y for x in db_nodes for y in x]
        cur.execute(cur.mogrify(statement, db_nodes_falttened))



        node_ids_bulk = cur.fetchall()


        # parse ids and map nodes to these
        node_ids = [x[0] for x in node_ids_bulk]




        sources = [node_ids[x.source] for x in graph_entry.es[:]]
        targets = [node_ids[x.target] for x in graph_entry.es[:]]
        mono_weights = graph_entry.es[:]["mono_weight"]
        end_mono_weights = graph_entry.es[:]["mono_end_weight"]
        avrg_weights = graph_entry.es[:]["avrg_weight"]
        end_avrg_weights = graph_entry.es[:]["avrg_end_weight"]
        attributes = [json.dumps(_get_qualifiers(x["qualifiers"])) for x in graph_entry.es[:]]


        cur.executemany("""INSERT INTO edges(source, target, mono_weight, mono_end_weight, avrg_weight, avrg_end_weight, attributes) VALUES 
        (%s, %s, %s, %s, %s, %s, %s)
        ;""", zip(sources, targets, mono_weights, end_mono_weights, avrg_weights, end_avrg_weights, attributes))



        conn.commit()


        output_query.put( (graph_entry.vs[1]["accession"], num_of_paths) )





    ### SQL ALchemy
    # engine = create_engine('postgresql://postgres:developer@127.0.0.1:5433/node_edges', echo=False)
    # Session = sessionmaker(bind=engine)
    # session = Session()



    # while True:
    #     try: 
    #         graph_entry = digested_graphs_queue.get(timeout=180)
    #     except Exception:
    #         continue

    #     # If table is not existing!
    #     if not engine.dialect.has_table(engine, "nodes"):
    #         Node.__table__.create(bind=engine)
    #     if not engine.dialect.has_table(engine, "edges"):
    #         # Node.__table__.create(bind=engine)
    #         # TODO
    #         pass

    #     # Generate Key list for json
    #     json_keys = [x for x in graph_entry.vs[0].attribute_names() if not( x == "accession" or x == "aminoacid" )]

    #     # Generate Node Object
    #     db_nodes = [Node(accession=x["accession"], aminoacid=x["aminoacid"], attributes=json.dumps({y: x[y] for y in json_keys})) for x in graph_entry.vs[:]]

    #     # Add to session and commit
    #     session.add_all(db_nodes)
    #     session.commit()


        # our_user = session.query(Node).filter_by(aminoacid='M').first() 







