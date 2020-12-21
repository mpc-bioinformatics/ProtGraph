

# Digestion of Complete uniprot yields roughly:
23 132 890 796 Edges and
8 675 257 837 Nodes


# Old Content
Example queries to retrieve paths from this datastructure (in postgres)




## Starting databae on Ubuntu
> /usr/lib/postgresql/10/bin/postgres -D $(pwd)/.db/pgsql/data -h 127.0.0.1 -p 5433 -k $(pwd)/.db/pgsql/run


## Get One Recursion Depth of a Path from ALL Proteins
```sql
	With cte_bfs(source, target, weight, paths) as
	(
	Select source, target, mono_weight, array[source, target]
	from edges
	where source = ANY (select id from nodes where nodes.aminoacid = '__start__')	
	and mono_weight > 500 -- Which range of weight do we want!?
	and mono_end_weight + edges.mono_weight < 1000 -- Usinng mono_end_weight do eliminate even more rows
	 )
	
	Select cte_bfs.source, cte_bfs.target, cte_bfs.weight + edges.mono_weight as weight, cte_bfs.paths || edges.target as paths
	from edges
	join cte_bfs
	on edges.source = cte_bfs.target
	where cte_bfs.weight + edges.mono_weight > 500 -- Which range of weight do we want!?
	and cte_bfs.weight + edges.mono_weight + edges.mono_end_weight< 1000
	
```


 This query below can retrieve them specific weight

```sql
with recursive cte_bfs(source, target, weight, paths ) as (
	Select source, target, mono_weight, array[source, target]
	from edges
	where source = ANY (select id from nodes where nodes.aminoacid = '__start__')	
	--and mono_weight > 5000 -- Which range of weight do we want!?
	and mono_end_weight + edges.mono_weight < 6000 -- Usinng mono_end_weight do eliminate even more rows
	
	Union 
	
	Select cte_bfs.source, cte_bfs.target, 
		cte_bfs.weight + edges.mono_weight as weight, cte_bfs.paths || edges.target as paths
	from edges
	join cte_bfs
	on edges.source = cte_bfs.target
	--where cte_bfs.weight + edges.mono_weight > 5000 -- Which range of weight do we want!?
	where cte_bfs.weight + edges.mono_weight + edges.mono_end_weight< 6000 -- Same here to reduce the number of rows
)

Select * from cte_bfs
where cte_bfs.weight > 5000
and array_length(cte_bfs.paths, 1) = 4
limit 1000
```

FINAL "OPTIMIZED" VERSION
```sql
explain analyze

with recursive cte_bfs(source, target, weight, paths ) as (
	Select source, target, mono_weight, array[source, target]
	from edges
	where source = ANY (select id from nodes where nodes.aminoacid = '__start__')	
	--and mono_weight > 5000 -- Which range of weight do we want!?
	and mono_end_weight + edges.mono_weight < 6000 -- Usinng mono_end_weight do eliminate even more rows
	
	Union 
	
	Select cte_bfs.source, cte_bfs.target, 
		cte_bfs.weight + edges.mono_weight as weight, cte_bfs.paths || edges.target as paths
	from edges
	join cte_bfs
	on edges.source = cte_bfs.target
	--where cte_bfs.weight + edges.mono_weight > 5000 -- Which range of weight do we want!?
	where cte_bfs.weight + edges.mono_weight + edges.mono_end_weight< 6000 -- Same here to reduce the number of rows
)

Select * from cte_bfs
where cte_bfs.weight > 5000
and array_length(cte_bfs.paths, 1) = 4
```

---


Get a path from to:
```sql
with recursive params as (
  select 12009068::bigint as fromnode,
         12009147::bigint as tonode
),
paths as (
    select ARRAY[fromnode] pathnodes,
           fromnode lastnode,
           0::double precision sumdistance
    from   params
  union all
    select     pathnodes || e.target,
               e.target,
               sumdistance + e.mono_weight
    from       paths
    join       edges e on e.source = lastnode
    cross join params p
    where      e.source <> p.tonode
    and        e.target <> all(pathnodes)
)
select   pathnodes, sumdistance
from     paths, params
where    lastnode = tonode
order by sumdistance desc
limit 10
```

All Paths for Protein: P54296
```sql
with recursive params as (
  select 12009068::bigint as fromnode,
         12009147::bigint as tonode
),
paths as (
    select ARRAY[fromnode] pathnodes,
           fromnode lastnode,
           0::double precision sumdistance
    from   params
  union all
    select     pathnodes || e.target,
               e.target,
               sumdistance + e.mono_weight
    from       paths
    join       edges e on e.source = lastnode
    cross join params p
    where      e.source <> p.tonode
    and        e.target <> all(pathnodes)
)
select   pathnodes, sumdistance
from     paths, params
where    lastnode = tonode
order by sumdistance desc
```

# Complete Recursion do find arbitrary lenghts of paths 

```sql
with recursive cte_bfs(source, target, weight, paths ) as (
	Select source, target, mono_weight, array[source, target]
	from edges
	where source = ANY (select id from nodes where nodes.aminoacid = '__start__')	
	and mono_weight > 500 -- Which range of weight do we want!?
	and mono_end_weight + edges.mono_weight < 1000 -- Usinng mono_end_weight do eliminate even more rows
	
	Union 
	
	Select cte_bfs.source, cte_bfs.target, 
		cte_bfs.weight + edges.mono_weight as weight, cte_bfs.paths || edges.target as paths
	from edges
	join cte_bfs
	on edges.source = cte_bfs.target
	where cte_bfs.weight + edges.mono_weight > 500 -- Which range of weight do we want!?
	and cte_bfs.weight + edges.mono_weight + edges.mono_end_weight< 1000 -- Same here to reduce the number of rows
)

Select * from cte_bfs

```


## Optimization Tricks

This will create an index for all proteins, where they start (for fast look up)
> create index on nodes (id) where aminoacid = '__start__'
> create index on nodes (id) where aminoacid = '__end__'





## Path with limiting max weight
```sql
with recursive params as (
  select 12009068::bigint as fromnode, -- start node
         12009147::bigint as tonode, -- end node
		 80000.00::double precision as max_weight -- max weight of path ( here mono_weight )
	
),
paths as (
    select ARRAY[fromnode] pathnodes,
           fromnode lastnode,
           0::double precision sumdistance
    from   params
  union all
    select     pathnodes || e.target,
               e.target,
               sumdistance + e.mono_weight
    from       paths
    join       edges e on e.source = lastnode
    cross join params p
    where      e.source <> p.tonode
    and        e.target <> all(pathnodes) -- maybe can be removed, since we have here strict DAGs !?
	and 	   sumdistance + e.mono_weight + e.mono_end_weight < p.max_weight -- Condition to only check for paths up to a weight
)
select   pathnodes, sumdistance
from     paths, params
where    lastnode = tonode
order by sumdistance desc
```


## Traverse a Graph from end to start

```sql
WITH RECURSIVE traverse(id, depth) AS (
        SELECT source, 124515142 FROM edges
        WHERE target = ANY (select id from nodes where nodes.aminoacid = '__end__')
    UNION ALL
        SELECT edges.source, traverse.depth + 1 FROM edges
        JOIN traverse
        ON edges.source = traverse.id
)
SELECT id, depth FROM traverse
ORDER BY depth DESC;
```





## UNTESTED 
# All paths with max weight
```sql
with recursive params as ( -- GET ALL PATHS from ALL PROTEINS!
	select n1.id as fromnode, n2.id as tonode, 300.00::double precision as max_weight
	from nodes n1
	join nodes n2 on n1.accession = n2.accession
	where n1.aminoacid = '__start__'
	and n2.aminoacid = '__end__'
),
paths as (
    select ARRAY[fromnode] pathnodes,
           fromnode lastnode,
           0::double precision sumdistance
    from   params
  union all
    select     pathnodes || e.target,
               e.target,
               sumdistance + e.mono_weight
    from       paths
    join       edges e on e.source = lastnode
    cross join params p
    where      e.source <> p.tonode
    and        e.target <> all(pathnodes)
	and 	   sumdistance + e.mono_weight + e.mono_end_weight < p.max_weight
)
select   pathnodes, sumdistance
from     paths, params
where    lastnode = tonode
order by sumdistance desc
```






















## Cypher and RedisGraph


Get all paths (limit 3) where the mono weight is less than 300.0
```cypher
MATCH p = (n : node {aminoacid : '__start__'})-[e:edge*]->(t : node {aminoacid : '__end__'}) 
WITH n, t, sum([edge IN relationships(p) | edge.mono_weight]) as s 
WHERE s < 300.0 RETURN * limit 3
```