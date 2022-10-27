"""
File contains some to get some distance calculations done for snapperdb v3.

original author: ulf.schaefer@phe.gov.uk
new author: jonathan.jenkins3@wales.nhs.uk (get_distances_new, get_all_pw_dists_new)

"""

import logging
from time import time
import json
import gzip
from operator import itemgetter

# --------------------------------------------------------------------------------------------------

def get_distances_new(conn, cur, test_sample_id, comp_id_list):
    """
    Get the distances of this sample to the other samples from the database.

    Parameters
    ----------
    cur: obj
        database cursor
    test_sample_id: int
        sample pk_id
    comp_id_list: list of int
        other samples to calculate the distance to

    Returns
    -------
    d: list of tuples
        sorted list of tuples with (sample_id, distance) with closes sample first
        e.g. [(298, 0), (37, 3), (55, 4)]
        None if fail
    """
    # Check for pre-calculated distances
    comp_id_set=set(comp_id_list)
    sql = "SELECT * FROM distances WHERE id_1={0} OR id_2={0}".format(test_sample_id)
    cur.execute(sql)
    rows = cur.fetchall()
    done=[]
    sample_dist_dic={}
    for row in rows:
        if int(row[0]) == test_sample_id:
            result_key=int(row[1])
        else:
            result_key=int(row[0])
        # filter to only samples were looking for
        if result_key in comp_id_set:
            sample_dist_dic[result_key]=int(row[2])
    done=set(sample_dist_dic.keys())

    # check stil samples reaminaing
    remaining_set=comp_id_set.difference(done)
    if len(remaining_set) < 1:
        d = sorted(sample_dist_dic.items(), key=itemgetter(1), reverse=False)
        return d

    # chunk
    remaining_list=list(remaining_set)
    chunk_size=750
    chunk_list=[]
    for i in range(0, len(remaining_list), chunk_size):
        chunk_list.append(remaining_list[i:i + chunk_size])

    # get list of contig ids from database
    sql = "SELECT pk_id FROM contigs"
    cur.execute(sql)
    rows = cur.fetchall()
    contig_ids = [r['pk_id'] for r in rows]

    # calculate the missing distances
    dist_missing_dic = {}
    chunk_count=0
    for small_chunk in chunk_list:
        chunk_count=chunk_count+1
        for cid in contig_ids:
            t0 = time()
            cur.callproc("get_sample_distances_by_id", [test_sample_id, cid, small_chunk])
            result = cur.fetchall()
            t1 = time()
            logging.info("Calculated %i (C%i) distances on contig %i with 'get_sample_distances_by_id' in %.3f seconds", len(result), chunk_count, cid, t1 - t0)

            # sum up if there are more than one contigs
            for res in result:
                if res[2] == None:
                    res[2] = 0
                try:
                    dist_missing_dic[res[0]] += res[2]
                except KeyError:
                    dist_missing_dic[res[0]] = res[2]

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
    args_str = b','.join(cur.mogrify('(%s,%s,%s)', x) for x in ordered_tup)
    # sql = "INSERT INTO distances VALUES " + args_str.decode("utf-8") 
    # cur.execute(sql)
    # conn.commit()
    
    # convert to nested tuple list output
    d = sorted(sample_dist_dic.items(), key=itemgetter(1), reverse=False)
    return d

def get_all_pw_dists_new(conn, cur, comp_id_list):
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

    comp_id_set=set(comp_id_list)
    done = set()
    dists_list = []

    for sample_id in comp_id_set:
        # note the ones that are already done, so nothing is calculated twice.
        done.add(sample_id)
        remaining_set = comp_id_set.difference(done)
        if len(remaining_set) < 1:
            break
        result = get_distances_new(conn, cur, sample_id, list(remaining_set))
        for pair in result:
            dists_list.append(pair[1])
    
    assert len(dists_list) == (len(comp_id_set) * (len(comp_id_set)-1))/2

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

def get_relevant_distances(conn, cur, sample_id):
    """
    Get the distances to this sample from the database.

    Parameters
    ----------
    cur: obj
        database cursor
    sample_id: int
        sample pk_id

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

    d = get_distances_new(conn, cur, sample_id, relv_samples)

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
