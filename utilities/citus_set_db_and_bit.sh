#!/bin/bash
# The password for the postgres user
CLUSTER_PWD="<Password>" 
# All postgresql machines with citus
CLUSTER_MACHINES=("192.168.0.3", "192.168.0.4", "192.168.0.5")
# node names (shards)
CLUSTER_NODES=(worker-01, worker-02)
CLUSTER_CONTROLLER="192.168.0.3"

# Set here the database name
DB_NAME=proteins

create_database(){
    # Create database on controller and nodes (if not exists) and install citus extension
    for cluster_machine in "${CLUSTER_MACHINES[@]}"
    do
        PGPASSWORD=${CLUSTER_PWD} psql -h $cluster_machine -p 5432 -U postgres <<-EOSQL
            SELECT 'CREATE DATABASE ${DB_NAME}' WHERE NOT EXISTS (SELECT FROM pg_database WHERE datname = '${DB_NAME}')\\gexec
EOSQL
        PGPASSWORD=${CLUSTER_PWD} psql -h $cluster_machine -p 5432 -U postgres -d ${DB_NAME} -c "CREATE EXTENSION IF NOT EXISTS citus;"
    done

    # Add nodes to controller
    for cluster_node in "${CLUSTER_NODES[@]}"
    do
        PGPASSWORD=${CLUSTER_PWD} psql -h $CLUSTER_CONTROLLER -p 5432 -U postgres -d ${DB_NAME} -c "SELECT * from master_add_node('${cluster_node}', 5432);"
    done

    # Print nodes again
    PGPASSWORD=${CLUSTER_PWD} psql -h $CLUSTER_CONTROLLER -p 5432 -U postgres -d ${DB_NAME} -c "SELECT * FROM master_get_active_worker_nodes();"
}


create_bits(){
    # Create bit hash function for distributing
    for cluster_machine in "${CLUSTER_MACHINES[@]}"
    do
        PGPASSWORD=${CLUSTER_PWD} psql -h $cluster_machine -p 5432 -U postgres -d ${DB_NAME} <<-EOSQL
            CREATE FUNCTION equal_test_bit_function(bit(452), bit(452)) RETURNS boolean
            AS 'select \$1 = \$2;'
            LANGUAGE SQL
            IMMUTABLE
            RETURNS NULL ON NULL INPUT;
EOSQL
        PGPASSWORD=${CLUSTER_PWD} psql -h $cluster_machine -p 5432 -U postgres -d ${DB_NAME} <<-EOSQL
            CREATE OPERATOR = (
                LEFTARG = bit(452),
                RIGHTARG = bit(452),
                PROCEDURE = equal_test_bit_function,
                HASHES
            );
EOSQL
        PGPASSWORD=${CLUSTER_PWD} psql -h $cluster_machine -p 5432 -U postgres -d ${DB_NAME} <<-EOSQL
            CREATE FUNCTION bit_hash(bit(452)) RETURNS int
            AS 'SELECT hashtext(\$1::text);'   
            LANGUAGE SQL
            IMMUTABLE
            RETURNS NULL ON NULL INPUT;
EOSQL
        PGPASSWORD=${CLUSTER_PWD} psql -h $cluster_machine -p 5432 -U postgres -d ${DB_NAME} <<-EOSQL
            CREATE OPERATOR CLASS new_op_fam_hash_class
            DEFAULT FOR TYPE bit(452) USING HASH AS
            OPERATOR 1 = (bit(452), bit(452)),
            FUNCTION 1 bit_hash(bit(452));
EOSQL
        PGPASSWORD=${CLUSTER_PWD} psql -h $cluster_machine -p 5432 -U postgres -d ${DB_NAME} <<-EOSQL
            CREATE OPERATOR CLASS new_op_fam_btree_class
            DEFAULT FOR TYPE bit(452) USING BTREE AS
            OPERATOR 3 = (bit(452), bit(452));
EOSQL
    done
}

drop_database (){
    # Drop database if exists on controller and nodes
    # Also remove all functions which may have been generated
    for cluster_machine in "${CLUSTER_MACHINES[@]}"
    do
        PGPASSWORD=${CLUSTER_PWD} psql -h $cluster_machine -p 5432 -U postgres -c "DROP DATABASE IF EXISTS ${DB_NAME};"
        PGPASSWORD=${CLUSTER_PWD} psql -h $cluster_machine -p 5432 -U postgres -c "DROP operator class new_op_fam_hash_class using hash;"
        PGPASSWORD=${CLUSTER_PWD} psql -h $cluster_machine -p 5432 -U postgres -c "DROP function bit_hash;"
        PGPASSWORD=${CLUSTER_PWD} psql -h $cluster_machine -p 5432 -U postgres -c "drop function equal_test_bit_function cascade;"
    done
}

echo "Dropping database"
drop_database
# Give the database some time
echo "Sleeping before creating database and functions"
sleep 5
create_database 
create_bits
echo "Finished!"
