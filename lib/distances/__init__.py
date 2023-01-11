"""
This file contains functions to use the get distance database calculations. 
First written for snapperdb v3, updated for version 3-5.

original author: ulf.schaefer@phe.gov.uk
update author: jonathan.jenkins3@wales.nhs.uk (
    get_distances_para, get_all_pw_dists_new,
    independent_dist_calc, chunk_remaining_list)
"""

import logging
from time import time
import json
import gzip
from operator import itemgetter
from multiprocessing import Pool
from contextlib import closing
import psycopg2
from psycopg2.extras import DictCursor
# --------------------------------------------------------------------------------------------------

def independent_dist_calc(input_dic):
    """
    Sets up an independent database connection to run get_sample_distances_by_id.
    get_sample_distances_by_id calculates the distance of one sample (test_sample_id)
    to a list of other samples (small_chunk). 

    Parameters
    ----------
    input_dic: dictionary 
        dictionary of arguments:
            small_chunk,    : list of int
                list of sample ids to compare distance to test_sample_id
                [comparison_sample_id_1, comparison_sample_id_2, ...]
            test_sample_id, : int
            contig_ids,     : int
            chunk_count,    : int
            pool_size,      : int
            db : dictionary 
                user : str
                host : str
                password : str
                dbname : str
    Returns
    -------
    dist_dic: dictionary
        A dictionary of distance results from the test_sample_id to other 
        samples (e.g. comparison_sample_id_1) found in small_chunk list.
        distance_1 is distance between test_sample_id and comparison_sample_id_1
        { comparison_sample_id_1: distance_1,
          comparison_sample_id_2: distance_2, ... }
        { 953 : 52, 
          1000: 1523, ....}
    """
    db_dic=input_dic['db']
    # make new independent database connection, need for parallel running
    conn1 = psycopg2.connect(user=db_dic['user'],host=db_dic['host'],password=db_dic['password'],dbname=db_dic['dbname'])
    cur1 = conn1.cursor(cursor_factory=psycopg2.extras.DictCursor)
    dist_dic={}
    for cid in input_dic['contig_ids']:
        # send calculation to database 
        t0 = time()
        cur1.callproc("get_sample_distances_by_id", [input_dic['test_sample_id'], cid, input_dic['small_chunk']])
        result = cur1.fetchall()
        t1 = time()
        logging.info("Calculated %i (C%i) distances on contig %i with 'get_sample_distances_by_id' in %.3f seconds",
            len(result), input_dic['chunk_count'], cid, t1 - t0)
        # sum up if there is more than one contig
        for res in result:
            if res[2] == None:
                res[2] = 0
            try:
                dist_dic[res[0]] += res[2]
            except KeyError:
                dist_dic[res[0]] = res[2]
    cur1.close()
    conn1.close()
    return dist_dic

def chunk_remaining_list(remaining_list,pool_size=4,min_chunk_size=8,max_chunk_size=750):
    """
    Split a list in to chunks to optimise the use of parallel pools and place limits the 
    size of the distance queries.

    Parameters
    ----------
    remaining_list: list of int
        other samples to calculate the distance to
    pool_size: int
        multiprocessing pool size to be passed on to any distance calculation
    min_chunk_size: int
        Added to limit time wasted on small quick queries to the database 
        generally min_chunk_size should be 5-20
    max_chunk_size: int
        Set greater than double the min_chunk_size.
        Looking at previous calculation times there is a trend for the 
        distance calculation take longer per samples when number of samples
        inputted was greater than about 750.

    Returns
    -------
    chunk_list: list of lists
        The original list of integers chunked in to smaller lists 
        [[1,2,3],[4,5,6],[7,8,9],[10,11,12]]
    out_pool_size: int 
        only changes if chunk_size < min_chunk_size
    """
    # Setting initial chunk size
    chunk_size=int(len(remaining_list)/pool_size)
    pool_runs=1
    if chunk_size > max_chunk_size:
        # calculate how many chucks would be need at the max chunk_size
        # adding one to the chunks_at_max if there is a remainder
        chunks_at_max=int(len(remaining_list)/max_chunk_size)+(len(remaining_list) % max_chunk_size > 0)
        # calculate number of full pool submission need 
        pool_runs=int(chunks_at_max/pool_size)+(chunks_at_max % pool_size > 0)
        # optimise the chunk size for the the total number of submitted queries
        chunk_size=int(len(remaining_list)/(pool_size*pool_runs))
    elif chunk_size < min_chunk_size:
        # reduce the pool size if possible 
        reduced_pool_size=int(len(remaining_list)/min_chunk_size)
        if pool_size > 1 and reduced_pool_size > 1:
            # optimise the chunk size for reduced_pool_size
            chunk_size=int(len(remaining_list)/(reduced_pool_size*pool_runs))
            pool_size=reduced_pool_size
        else:   
            # too small for more than one pool
            chunk_size=len(remaining_list)
            pool_size=1
    # calculate how many list items would be remaining at with current chunk size
    num_items_remaining=len(remaining_list)-(pool_size*pool_runs*chunk_size)
    # make a list with desired chunk size total list length =(pool_size*pool_runs)
    chunk_size_list=[chunk_size] * ((pool_size*pool_runs)-num_items_remaining)
    # spread the remainder evenly over multiple chunks 
    chunk_size_list.extend([chunk_size+1] * num_items_remaining)
    # chunk the input list into chunks sizes in found chunk_size_list
    chunk_list=[]
    index_c=0
    for test_chunk_size in chunk_size_list:
        chunk_list.append(remaining_list[index_c:index_c+test_chunk_size])
        index_c=index_c+test_chunk_size
    out_pool_size=pool_size
    return chunk_list, out_pool_size

def get_distances_para(conn, cur, test_sample_id, comparison_sample_ids_list, pool_size=4):
    """
    Get the distances of this sample to the other samples from the database.
    Calculations with large numbers of distances to calculate are split
    into chunks that can run in parallel. The calculated distances are then stored
    in the distances table. An UPSERT (INSERT ... ON CONFLICT UPDATE) is used to 
    allow distances to be safely calculated by two different processes in the future.

    Parameters
    ----------
    conn: obj
        database connection
    cur: obj
        database cursor
    test_sample_id: int
        sample pk_id
    comparison_sample_ids_list: list of int
        list of other sample pk_id to calculate the distance to test_sample_id.
    pool_size: int
        multiprocessing pool size to be passed on to any distance calculation

    Returns
    -------
    d: list of tuples
        sorted list of tuples with (comparison_sample_id_1, distance) with closest sample first.
        The distance is between test_sample_id and other comparison_sample_id_1
        e.g. [(298, 0), (37, 3), (55, 4)]      
    """
    # Check for pre-calculated distances
    comp_sample_ids_set=set(comparison_sample_ids_list)
    sql = "SELECT * FROM distances WHERE id_1={0} OR id_2={0}".format(test_sample_id)
    cur.execute(sql)
    rows = cur.fetchall()
    sample_ids_with_calculated_distances=[]
    sample_dist_dic={}
    for row in rows:
        if int(row[0]) == test_sample_id:
            result_key=int(row[1])
        elif int(row[1]) == test_sample_id:
            result_key=int(row[0])
        else:
            raise ValueError(
                'ERROR: test_sample_id ({}) must match one id in db query return ({} or {})'.format(
                    test_sample_id,row[0],row[1]))
        # filter to only samples were looking for
        if result_key in comp_sample_ids_set:
            sample_dist_dic[result_key]=int(row[2])
    sample_ids_with_calculated_distances=set(sample_dist_dic.keys()) # 

    # remove samples that have a distance pre-calculated to not calculate again
    remaining_set=comp_sample_ids_set.difference(sample_ids_with_calculated_distances)
    if len(remaining_set) < 1:
        d = sorted(sample_dist_dic.items(), key=itemgetter(1), reverse=False)
        return d

    # split list in to optimised chunks
    remaining_list=list(remaining_set)
    chunk_list, out_pool_size=chunk_remaining_list(remaining_list,pool_size=pool_size)

    # get list of contig ids from database
    sql = "SELECT pk_id FROM contigs"
    cur.execute(sql)
    rows = cur.fetchall()
    contig_ids = [r['pk_id'] for r in rows]
   
    # calculate the missing distances
    dist_missing_dic = {}
    chunk_count=0
    pool_arg_list=[]
    for small_chunk in chunk_list:
        chunk_count=chunk_count+1
        input_dic = {
            'small_chunk'    : small_chunk,
            'test_sample_id' : test_sample_id,
            'contig_ids'     : contig_ids,
            'chunk_count'    : chunk_count,
            'pool_size'      : out_pool_size,
            'db' : {
                'user'    : conn.info.user,
                'host'    : conn.info.host,
                'password': conn.info.password,
                'dbname'  : conn.info.dbname
                }
        }
        pool_arg_list.append(input_dic)
    with closing(Pool(processes=out_pool_size)) as p:
        result_list = p.map(independent_dist_calc, pool_arg_list)
    for result in result_list:
        dist_missing_dic.update(result)

    # add missing distances to distances table and dist_missing_dic
    # distances (small_id, large_id, distance)
    ordered_list=[]
    for key, item in dist_missing_dic.items():
        comp_id=int(key)
        dist_test=int(item)
        sample_dist_dic[comp_id]=dist_test
        if test_sample_id < comp_id:
            o_list=[test_sample_id,comp_id,dist_test]
        else:
            o_list=[comp_id,test_sample_id,dist_test]
        ordered_list.append(tuple(o_list))
    ordered_tup=tuple(ordered_list)
    args_str = str(','.join(str(x) for x in ordered_tup))
    sql = "INSERT INTO distances VALUES {} ON CONFLICT (id_1, id_2) DO UPDATE SET distance = excluded.distance".format(
        args_str)
    cur.execute(sql)
    conn.commit()
    
    # convert to nested tuple list output
    d = sorted(sample_dist_dic.items(), key=itemgetter(1), reverse=False)
    return d

def get_all_pw_dists_new(conn, cur, comparison_sample_ids_list, pool_size=4):
    """
    Get all pairwise distances between the samples in the input list. 

    Parameters
    ----------
    conn: obj
        database connection
    cur: obj
        database cursor
    comparison_sample_ids_list: list of int
        list of sample pk_id
    pool_size: int
        multiprocessing pool size to be passed on to any distance calculation
    Returns
    -------
    dists: lists of ints
        lists with distances
    None if there is a problem
    """

    comp_sample_ids_set=set(comparison_sample_ids_list)
    done = set()
    dists_list = []

    for sample_id in comp_sample_ids_set:
        # note the ones that are already done, so nothing is calculated twice.
        done.add(sample_id)
        remaining_set = comp_sample_ids_set.difference(done)
        if len(remaining_set) < 1:
            break
        result = get_distances_para(conn, cur, sample_id, list(remaining_set), pool_size)
        for pair in result:
            # get distance from list of tuples produced by get_distances_para
            dists_list.append(pair[1])
    
    assert len(dists_list) == (len(comp_sample_ids_set) * (len(comp_sample_ids_set)-1))/2

    return dists_list


# --------------------------------------------------------------------------------------------------

def get_all_pw_dists(cur, samids):
    """
    Get all pairwise distances between the samples in the input list.


    Parameters
    ----------
    cur: obj
        database cursor
    samids: list of int
        sample ids

    Returns
    -------
    dists: lists of ints
        lists with distances
    None if there is a problem
    """

    # get list of contig ids from database
    sql = "SELECT pk_id FROM contigs"
    cur.execute(sql)
    rows = cur.fetchall()
    contig_ids = [r['pk_id'] for r in rows]

    samids = set(samids)
    done = set()
    dists = []

    for s in samids:
        # note the ones that are already done, so nothing is calculated twice.
        done.add(s)
        oths = samids.difference(done)

        d = {}
        for cid in contig_ids:

            cur.callproc("get_sample_distances_by_id", [s, cid, list(oths)])
            result = cur.fetchall()

            for res in result:
                if res[2] == None:
                    res[2] = 0
                try:
                    d[res[0]] += res[2]
                except KeyError:
                    d[res[0]] = res[2]

        dists += d.values()

    assert len(dists) == (len(samids) * (len(samids)-1))/2

    return dists

# --------------------------------------------------------------------------------------------------

def get_relevant_distances(conn, cur, sample_id, pool_size=4):
    """
    Get the distances to this sample from the database.

    Parameters
    ----------
    conn: obj
        database connection
    cur: obj
        database cursor
    sample_id: int
        sample pk_id
    pool_size: int
        multiprocessing pool size to be passed on to any distance calculation

    Returns
    -------
    d: list of tuples
        sorted list of tuples with (sample_id, distance) with closes sample first
        e.g. [(298, 0), (37, 3), (55, 4)]
        None if fail
    """

    # get the relevant samples from the database, these are the ones that have been clustered and are not ignored
    sql = "SELECT c.fk_sample_id FROM sample_clusters c, samples s WHERE s.pk_id=c.fk_sample_id AND s.ignore_sample IS FALSE"
    cur.execute(sql)
    rows = cur.fetchall()
    relv_samples = [r['fk_sample_id'] for r in rows]

    d = get_distances_para(conn, cur, sample_id, relv_samples, pool_size)

    return d

# --------------------------------------------------------------------------------------------------

def get_missing_distances(cur, sample_id, haves):
    """
    Get the distances to this sample that we don't already have from the database.

    Parameters
    ----------
    cur: obj
        database cursor
    sample_id: int
        sample pk_id
    haves: set
        set of sample ids tuples for which we already have the distance.

    Returns
    -------
    d: list of tuples
        sorted list of tuples with (sample_id, distance) with closes sample first
        e.g. [(298, 0), (37, 3), (55, 4)]
        None if fail
    """

    # get the relevant samples from the database, these are the ones that have been clustered and are not ignored
    sql = "SELECT c.fk_sample_id FROM sample_clusters c, samples s WHERE s.pk_id=c.fk_sample_id AND s.ignore_sample IS FALSE"
    cur.execute(sql)
    rows = cur.fetchall()
    all_samples = set([r['fk_sample_id'] for r in rows])

    relv_samples = list(all_samples.difference(haves))

    d = get_distances(cur, sample_id, relv_samples)

    return d


# --------------------------------------------------------------------------------------------------

def get_distances(cur, samid, others):
    """
    Get the distances of this sample to the other samples from the database.

    Parameters
    ----------
    cur: obj
        database cursor
    sample_id: int
        sample pk_id
    others: list of int
        other samples to calculate the distance to

    Returns
    -------
    d: list of tuples
        sorted list of tuples with (sample_id, distance) with closes sample first
        e.g. [(298, 0), (37, 3), (55, 4)]
        None if fail
    """

    # get list of contig ids from database
    sql = "SELECT pk_id FROM contigs"
    cur.execute(sql)
    rows = cur.fetchall()
    contig_ids = [r['pk_id'] for r in rows]

    d = {}
    for cid in contig_ids:
        t0 = time()
        cur.callproc("get_sample_distances_by_id", [samid, cid, others])
        result = cur.fetchall()
        t1 = time()
        logging.info("Calculated %i distances on contig %i with 'get_sample_distances_by_id' in %.3f seconds", len(result), cid, t1 - t0)

        # sum up if there are more than one contigs
        for res in result:
            if res[2] == None:
                res[2] = 0
            try:
                d[res[0]] += res[2]
            except KeyError:
                d[res[0]] = res[2]

    d = sorted(d.items(), key=itemgetter(1), reverse=False)

    return d

# --------------------------------------------------------------------------------------------------

def get_distance_matrix(cur, samids):
    """
    Get a distance matrix for the given samples.

    Parameters
    ----------
    cur: obj
        database cursor
    samids: int
        list of sample pk_ids

    Returns
    -------
    dist: dict
        complete matrix
        dist[s1][s2] = d
        dist[s2][s1] = d
    """

    # get list of contig ids from database
    sql = "SELECT pk_id FROM contigs"
    cur.execute(sql)
    rows = cur.fetchall()
    contig_ids = [r['pk_id'] for r in rows]

    samids = set(samids)
    done = set()
    dists = {}

    for s in samids:

        try:
            dists[s][s] = 0
        except KeyError:
            dists[s] = {s: 0}

        # note the ones that are already done, so nothing is calculated twice.
        done.add(s)
        oths = samids.difference(done)

        d = {}
        for cid in contig_ids:

            cur.callproc("get_sample_distances_by_id", [s, cid, list(oths)])
            result = cur.fetchall()

            for res in result:
                if res[2] == None:
                    res[2] = 0
                try:
                    d[res[0]] += res[2]
                except KeyError:
                    d[res[0]] = res[2]
        for osam, snpdi in d.items():
            dists[s][osam] = snpdi
            try:
                dists[osam][s] = snpdi
            except KeyError:
                dists[osam] = {s: snpdi}

    return dists

# --------------------------------------------------------------------------------------------------

def get_distances_precalc(cur, sam_id, sample_name, json_file_name):
    """
    Check the precalculated data against the database, get the missing distances,
    and add the precalculated ones.

    Parameters
    ----------
    cur: obj
        database cursor
    sample_id: int
        sample pk_id
    sample_name: str
        samples name
    json_file_name: str
        name of a json file containing

    Returns
    -------
    d: list of tuples
        sorted list of tuples with (sample_id, distance) with closes sample first
        e.g. [(298, 0), (37, 3), (55, 4)]
        None if fail
    """

    open_func = gzip.open if json_file_name.endswith('.gz') == True else open
    precalc_data = None
    with open_func(json_file_name) as f:
        precalc_data = json.load(f)

    # chck the sample name in the json against the current sample
    if precalc_data['sample_name'] != sample_name:
        logging.error("Sample name does not match precalculated data!")
        return None

    # check data consistency between the db and the precalculated distances, do samid and samname match?
    sql = "SELECT c.fk_sample_id AS sid, s.sample_name AS name FROM sample_clusters c, samples s WHERE s.pk_id=c.fk_sample_id AND s.ignore_sample IS FALSE"
    cur.execute(sql)
    rows = cur.fetchall()
    tbl_values = {r['sid']: r['name'] for r in rows}
    for (pre_id, pre_name, dis) in precalc_data['distances']:
        if tbl_values[pre_id] != pre_name:
            logging.error("Precalculated data does not match database. Precalculated samples name for id %i was %s, but in db it's %s",
                          pre_id, pre_name, tbl_values[pre_id])
            return None

    # get missing distances, add the precalculated ones and re-sort
    d = get_missing_distances(cur, sam_id, set([x[0] for x in precalc_data['distances']]))
    d += [(x[0], x[2]) for x in precalc_data['distances']]
    d.sort(key=lambda x: x[1])

    return d

# --------------------------------------------------------------------------------------------------

